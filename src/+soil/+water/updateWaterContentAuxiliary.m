function built = updateWaterContentAuxiliary(built, dt_s, L_f)
%UPDATEWATERCONTENTAUXILIARY  Refresh soil-water auxiliary cache from state.
%
%   built = soil.water.updateWaterContentAuxiliary(built, dt_s, L_f) prepares
%   layer-endpoint inputs, calls soil.water.computeEndpointHydraulics, and
%   writes the returned auxiliary arrays into built.workspace.derived.
%
%   Inputs:
%     built  [inout] struct canonical model container with built.Y,
%                    built.parameters.soil, built.control, and workspace
%     dt_s   [in]    [s] outer update interval. Use zero/omit during
%                    initialization; Dall'Amico phase change only runs when
%                    dt_s>0, L_f~=0, and thermal fields are populated.
%     L_f    [in]    [J/kg] fusion latent heat coefficient used by the unified
%                    hydraulic kernel. If omitted, the site freeze switch
%                    supplies the default coefficient for compatibility.
%
%   Output:
%     built  [out] struct with refreshed workspace.derived:
%                    COR/CORh/Gamma_hh, Theta_LL/UU/II/U/I/V/g, Se, KL_h,
%                    KfL_h, DTheta_LLh, DTheta_UUh and Ratio_ice
%
%   Role:
%     This is a thin scheduler. The physical water-content and hydraulic
%     formulas live in soil.water.computeEndpointHydraulics.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: 2026-07-11 Jizixing Zhu -- derive explicit liquid, ice, total-water,
%           and gas-volume partials for the conservative storage assembler.
%   Change: 2026-07-12 Jizixing Zhu -- anchor phase relaxation to accepted Y.Theta_I
%           so Picard refreshes do not apply one time-step phase change repeatedly.
%   Change: Reduced this routine to endpoint preparation and workspace writeback.
%   Change: 2026-07-13 Jizixing Zhu -- read Cw from the shared physical
%           constants used by both retention and thermal liquid flow.
%   Change: 2026-07-14 Jizixing Zhu -- align the documented auxiliary contract
%           with the strict-SFCC endpoint fields actually written here.

    if nargin < 2 || isempty(dt_s)
        dt_s = 0;
    end
    requireWorkspace(built);
    phys = equations.physicalConstants();
    if nargin < 3 || isempty(L_f)
        L_f = fusionLatentHeat(built.control.switches, phys);
    end

    nCell = readDomainField(built.domain, 'nCell', 1);
    for c = 1:nCell
        built = refreshColumn(built, c, L_f, dt_s, phys);
    end
end

function built = refreshColumn(built, cellIndex, L_f, dt_s, phys)
% Refresh one vertical soil column.
    NL = built.domain.NL;
    nD = readDomainField(built.domain.grid, 'nD', 2);

    hTrial = soil.water.getNodeColumn(built.workspace.step.trial.h, cellIndex);
    TTrial = soil.water.getNodeColumn(built.workspace.step.trial.T, cellIndex);
    hFrezTrial = soil.water.getNodeColumn(built.workspace.step.trial.h_frez, cellIndex);
    TTcritTrial = soil.water.getNodeColumn(built.workspace.derived.TT_CRIT, cellIndex);
    dhFrezDh = soil.water.getNodeColumn(built.workspace.derived.DhFrezDh, cellIndex);
    dhFrezDT = soil.water.getNodeColumn(built.workspace.derived.DhFrezDT, cellIndex);

    [COR, CORh, hCorrected] = surfaceTensionCorrectedHead(hTrial, TTrial, built.control.switches, phys);
    state = struct();
    state.h = endpointFromNodes(hCorrected, NL, nD);
    state.T_C = endpointFromNodes(TTrial, NL, nD);
    state.h_frez = endpointFromNodes(hFrezTrial, NL, nD);
    state.TT_CRIT = endpointFromNodes(TTcritTrial, NL, nD);
    state.theta_i_prev = acceptedIceColumn(built, cellIndex);
    state.dt_s = dt_s;
    state.EHCAP = getEndpointColumn(built.workspace.derived.EHCAP, cellIndex);
    state.c_unsat = getEndpointColumn(built.workspace.derived.c_unsat, cellIndex);
    state.Lambda_eff = getEndpointColumn(built.workspace.derived.Lambda_eff, cellIndex);
    state.dz_cm = built.domain.dz(:) * ones(1, nD);

    soilColumn = selectSoilParameterColumn(built.parameters.soil, cellIndex);
    endpoint = soil.water.computeEndpointHydraulics( ...
        soilColumn, built.control.switches, state, L_f, phys);
    rawState = struct( ...
        'h_raw', endpointFromNodes(hTrial, NL, nD), ...
        'COR', endpointFromNodes(COR, NL, nD), ...
        'CORh', endpointFromNodes(CORh, NL, nD));
    freezePartials = struct( ...
        'dh_frez_dh', endpointFromNodes(dhFrezDh, NL, nD), ...
        'dh_frez_dT', endpointFromNodes(dhFrezDT, NL, nD));
    storagePartials = soil.water.computePhaseStoragePartials( ...
        endpoint, rawState, built.control.switches, freezePartials, phys);

    d = built.workspace.derived;
    d.COR = soil.water.setNodeColumn(d.COR, cellIndex, COR);
    d.CORh = soil.water.setNodeColumn(d.CORh, cellIndex, CORh);
    d.Gamma_hh = soil.water.setNodeColumn(d.Gamma_hh, cellIndex, ...
        endpointGammaToNodes(endpoint.gamma_hh, NL, nD));
    d.Theta_LL = soil.water.setEndpointColumn(d.Theta_LL, cellIndex, endpoint.theta_l);
    d.Theta_UU = soil.water.setEndpointColumn(d.Theta_UU, cellIndex, endpoint.theta_u);
    d.Theta_II = soil.water.setEndpointColumn(d.Theta_II, cellIndex, endpoint.theta_i);
    d.Theta_U = soil.water.setEndpointColumn(d.Theta_U, cellIndex, endpoint.theta_u);
    d.Theta_I = soil.water.setEndpointColumn(d.Theta_I, cellIndex, endpoint.theta_i);
    d.Theta_V = soil.water.setEndpointColumn(d.Theta_V, cellIndex, endpoint.theta_v);
    d.Theta_g = soil.water.setEndpointColumn(d.Theta_g, cellIndex, endpoint.theta_v);
    d.Se = soil.water.setEndpointColumn(d.Se, cellIndex, endpoint.Se_air);
    d.KL_h = soil.water.setEndpointColumn(d.KL_h, cellIndex, endpoint.KL_h);
    d.KfL_h = soil.water.setEndpointColumn(d.KfL_h, cellIndex, endpoint.KfL_h);
    d.DTheta_LLh = soil.water.setEndpointColumn(d.DTheta_LLh, cellIndex, endpoint.dtheta_l_dh);
    d.DTheta_UUh = soil.water.setEndpointColumn(d.DTheta_UUh, cellIndex, endpoint.dtheta_u_dh);
    d.Ratio_ice = soil.water.setEndpointColumn(d.Ratio_ice, cellIndex, endpoint.ratio_ice);
    d.DallAmicoThermalTime = soil.water.setEndpointColumn( ...
        d.DallAmicoThermalTime, cellIndex, endpoint.dallamico_thermal_time);
    d.DallAmicoThetaLStar = soil.water.setEndpointColumn( ...
        d.DallAmicoThetaLStar, cellIndex, endpoint.dallamico_theta_l_star);
    d.DallAmicoPhaseSource = soil.water.setEndpointColumn( ...
        d.DallAmicoPhaseSource, cellIndex, endpoint.dallamico_phase_source);
    d.DallAmicoLimiterCode = soil.water.setEndpointColumn( ...
        d.DallAmicoLimiterCode, cellIndex, endpoint.dallamico_limiter_code);
    d.DThetaLiquidDh = soil.water.setEndpointColumn( ...
        d.DThetaLiquidDh, cellIndex, storagePartials.dtheta_liquid_dh);
    d.DThetaLiquidDT = soil.water.setEndpointColumn( ...
        d.DThetaLiquidDT, cellIndex, storagePartials.dtheta_liquid_dT);
    d.DThetaWaterEquivDh = soil.water.setEndpointColumn( ...
        d.DThetaWaterEquivDh, cellIndex, storagePartials.dtheta_water_equiv_dh);
    d.DThetaWaterEquivDT = soil.water.setEndpointColumn( ...
        d.DThetaWaterEquivDT, cellIndex, storagePartials.dtheta_water_equiv_dT);
    d.DThetaIceDh = soil.water.setEndpointColumn( ...
        d.DThetaIceDh, cellIndex, storagePartials.dtheta_ice_dh);
    d.DThetaIceDT = soil.water.setEndpointColumn( ...
        d.DThetaIceDT, cellIndex, storagePartials.dtheta_ice_dT);
    d.DThetaGasDh = soil.water.setEndpointColumn( ...
        d.DThetaGasDh, cellIndex, storagePartials.dtheta_gas_dh);
    d.DThetaGasDT = soil.water.setEndpointColumn( ...
        d.DThetaGasDT, cellIndex, storagePartials.dtheta_gas_dT);
    built.workspace.derived = d;
end

function thetaI = acceptedIceColumn(built, cellIndex)
% Read the immutable accepted ice state for this nonlinear trial.
    if isfield(built.Y, 'Theta_I') && ~isempty(built.Y.Theta_I)
        thetaI = getEndpointColumn(built.Y.Theta_I, cellIndex);
    else
        thetaI = getEndpointColumn(built.workspace.derived.Theta_II, cellIndex);
    end
end

function A = endpointFromNodes(x, NL, nD)
% Build [NL x nD] endpoint arrays from an [NN x 1] node vector.
    A = zeros(NL, nD);
    for j = 1:nD
        A(:, j) = x(j:(j + NL - 1));
    end
end

function gammaNode = endpointGammaToNodes(gamma, NL, nD)
% Map endpoint Gamma_hh back to the NN node vector.
    gammaNode = [gamma(1:NL, 1); gamma(NL, 2:nD).'];
end

function [COR, CORh, hCorrected] = surfaceTensionCorrectedHead(h, T_C, switches, phys)
% Apply the de Vries-style surface-tension correction to a trial head vector.
    COR = ones(size(h));
    CORh = zeros(size(h));
    if isfield(switches, 'lhThmrl') && switches.lhThmrl == 1
        CORh(:) = phys.Cw;
        COR = exp(-CORh .* (T_C - phys.Tr));
        COR(COR == 0) = 1;
    end
    hCorrected = COR .* h;
end

function requireWorkspace(built)
% Guard that the caller allocated workspace before refreshing auxiliary fields.
    if ~isfield(built, 'workspace') || ~isfield(built.workspace, 'derived') ...
            || ~isfield(built.workspace, 'step')
        error('updateWaterContentAuxiliary:missingWorkspace', ...
            'Call init.allocateWorkspace before soil.water.updateWaterContentAuxiliary.');
    end
end

function L_f = fusionLatentHeat(switches, phys)
% Return fusion latent heat only when soil-freezing physics is enabled.
    if isfield(switches, 'lFreezing') && switches.lFreezing == 1
        L_f = phys.lambda_fus;
    else
        L_f = 0;
    end
end

function value = readDomainField(s, name, defaultValue)
% Read a struct field with fallback.
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end

function E = getEndpointColumn(A, cellIndex)
% Read one layer-endpoint matrix from single-column or regional arrays.
    if ismatrix(A)
        E = A;
    else
        E = A(:, :, cellIndex);
    end
end

function soilColumn = selectSoilParameterColumn(soilParams, cellIndex)
% Select one regional soil-parameter column without changing scalar metadata.
    soilColumn = soilParams;
    for name = {'theta_s','theta_r','alpha','n','m','porosity','field_capacity','ks'}
        if isfield(soilColumn, name{1})
            soilColumn.(name{1}) = selectColumnIfRegional(soilColumn.(name{1}), cellIndex);
        end
    end
    if isfield(soilColumn, 'modes') && isfield(soilColumn.modes, 'freeze') ...
            && isfield(soilColumn.modes.freeze, 'imped')
        soilColumn.modes.freeze.imped = selectColumnIfRegional( ...
            soilColumn.modes.freeze.imped, cellIndex);
    end
end

function value = selectColumnIfRegional(value, cellIndex)
% Return A(:,cellIndex) for [nVertical x nCell] arrays; leave vectors/scalars.
    if isnumeric(value) && ~isvector(value) && size(value, 2) >= cellIndex
        value = value(:, cellIndex);
    end
end
