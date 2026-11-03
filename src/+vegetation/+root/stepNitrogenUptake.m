function result = stepNitrogenUptake(model, demandContext)
%STEPNITROGENUPTAKE  Compute independent current-step root NH4/NO3 uptake.
%
%   result = vegetation.root.stepNitrogenUptake(model, demandContext) computes the
%   soil mineral-N root sink without assembling or mutating PlantCN state.
%   lRsuef is the sole master switch. PlantCN may optionally provide demand
%   caps through demandContext.root_ctx. Without PlantCN demand, the module
%   falls back to the original passive max(Srt,0)*C mass-flow formulation.
%   Thus soil-N transport, PlantCN inventory, and passive root N removal are
%   separate processes: passive uptake is an untracked vegetation sink when
%   PlantCN is disabled, not a hidden PlantCN state update.
%
%   Inputs:
%     model          [in] struct  current staged model, including scope Srt
%     demandContext  [in] struct  optional output of preparePlantNitrogenDemand
%
%   Returns result:
%     enabled       [-]             lRsuef state
%     SinkS         [ug/cm3/s, NLx2x2] NH4/NO3 soil sinks
%     SPUptake      [ug/cm2/s, 1x2] depth-integrated NH4/NO3 uptake
%     U_N_total     [g N/m2/s]      total root N inflow
%     uptake_diag   struct           DAISY/Nacry layer diagnostics
%     demand_nh4/no3 [g N/m2/s]     optional PlantCN demand caps
%     scheme        char             off, passive_mass_flow, or daisy_nacry_demand
%     accounting_mode char           plant_cn_inventory or passive_sink_untracked
%
%   Physical closure:
%     SinkS is a pending soil-N flux. The accepted-step transaction commits it
%     once; rejected steps discard it. PlantCN can consume U_N_total but cannot
%     change SinkS or SPUptake.
%
%   Reference: soil.nitrogen.root_uptake.nitrogenUptake and
%   vegetation.root.applyRootUptake (DAISY single-root plus Nacry membrane coupling).
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
% Change: 2026-07-11 Jizixing Zhu -- use original water-linked passive uptake when PlantCN demand is absent.
% Change: 2026-07-20 Jizixing Zhu -- expose whether uptake enters PlantCN inventory or is the independent passive sink.

    if nargin < 2 || isempty(demandContext), demandContext = struct(); end
    built = model.built;
    NL = built.domain.NL;
    result = zeroResult(NL);
    if getFlag(built.control.switches, 'lRsuef', 0) ~= 1 ...
            || getFlag(built.control.switches, 'lNefc', 0) ~= 1
        return;
    end
    assert(isfield(built.workspace.scope_step, 'Srt') ...
        && ~isempty(built.workspace.scope_step.Srt), ...
        'stepNitrogenUptake:missingCurrentSrt', ...
        'Current scope_step.Srt is required before root nitrogen uptake.');

    if ~isfield(demandContext, 'root_ctx') || isempty(demandContext.root_ctx)
        [sinkS, spUptake, uptakeDiag] = vegetation.root.passiveNitrogenUptake( ...
            built.workspace.scope_step.Srt, built.domain.dz, ...
            built.Y.nh4, built.Y.no3);
        demandNh4 = 0;
        demandNo3 = 0;
        scheme = 'passive_mass_flow';
    else
        [daisyParams, nacryParams, demandNh4, demandNo3] = ...
            resolveDemandParameters(demandContext);
        soilState = init.assembleSoilStateN(built);
        rld = rootLengthDensity(built.p.plantHydro, built.domain);
        [~, sinkS, spUptake, uptakeDiag] = vegetation.root.applyRootUptake( ...
            built.workspace.scope_step.Srt, NL, built.domain.dz, ...
            built.domain.tot_depth, built.p.plantHydro.Lmax, 1, ...
            built.Y.nh4, built.Y.no3, soilState, rld, daisyParams, nacryParams);
        scheme = 'daisy_nacry_demand';
    end

    result.enabled = true;
    result.SinkS = sinkS;
    result.SPUptake = spUptake;
    result.U_N_total = sum(spUptake) * 1e-2;
    result.uptake_diag = uptakeDiag;
    result.demand_nh4 = demandNh4;
    result.demand_no3 = demandNo3;
    result.scheme = scheme;
    if strcmp(scheme, 'daisy_nacry_demand')
        result.accounting_mode = 'plant_cn_inventory';
    else
        result.accounting_mode = 'passive_sink_untracked';
    end
end

% -------------------------------------------------------------------------
function result = zeroResult(NL)
% Return the fixed result schema when root N uptake is disabled.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    result = struct('enabled', false, 'SinkS', zeros(NL, 2, 2), ...
        'SPUptake', zeros(1, 2), 'U_N_total', 0, ...
        'uptake_diag', struct(), 'demand_nh4', 0, 'demand_no3', 0, ...
        'scheme', 'off', 'accounting_mode', 'disabled');
end

% -------------------------------------------------------------------------
function [daisy, nacry, demandNh4, demandNo3] = resolveDemandParameters(context)
% Read the PlantCN demand and DAISY/Nacry parameter context.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    rootCtx = context.root_ctx;
    daisy = rootCtx.daisyParams;
    nacry = rootCtx.nacryParams;
    demandNh4 = rootCtx.demand_NH4_g_m2_s;
    demandNo3 = rootCtx.demand_NO3_g_m2_s;
end

% -------------------------------------------------------------------------
function rld = rootLengthDensity(plantHydro, domain)
% Resolve measured RLD or derive it from RAI without PlantCN state.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    if isfield(plantHydro, 'rootProfile') ...
            && isfield(plantHydro.rootProfile, 'root_length_density_m_m3') ...
            && ~isempty(plantHydro.rootProfile.root_length_density_m_m3)
        rld = plantHydro.rootProfile.root_length_density_m_m3(:);
    else
        dz_m = domain.dz(:) * 1e-2;
        rai = plantHydro.RAI(:);
        rai = rai(1:domain.NL);
        rld = rai ./ (2 * pi * plantHydro.r_root_m * dz_m);
    end
end

% -------------------------------------------------------------------------
function value = getFlag(s, name, fallback)
% Read one normalized switch.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    value = fallback;
    if isfield(s, name) && ~isempty(s.(name)), value = s.(name); end
end
