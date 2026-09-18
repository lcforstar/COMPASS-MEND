function built = refreshThermalGasTransport(built, dt_s, refreshMode)
%REFRESHTHERMALGASTRANSPORT  B2 pre-solve thermal/vapor/dry-air transport-coefficient refresh.
%
%   built = soil.refreshThermalGasTransport(built, dt_s) marshals built into a legacy-compatible
%   working frame, calls the 9 EXISTING kernels in COMPASS.m:1685-1694 order (thermal/vapor/dry-air
%   transport coefficients), and scatters the produced fields into built.workspace.derived and
%   built.workspace.flux.V_A. This is a REUSE marshaling adapter -- it does NOT reimplement the
%   thermal/gas physics; the kernels are read-only. It is a cross-domain soil pre-solve adapter,
%   so it lives at top-level +soil (not +thermal/+gas).
%
%   Order (pinned to COMPASS.m:1685-1694): calculateThermalLiquidConductivity; Density_V;
%   calculateThermalLiquidDispersion; Latent; Density_DA; calculateThermalConductivityCapacity;
%   calculateGasConductivity; calculateVaporVariables; calculateGasDispersivity.
%
%   Notes: EHCAP/c_unsat/Lambda_eff (#6) are the B1<->B2 coupling anchor (updateWaterContentAuxiliary
%   reads them for Dall'Amico at dt_s>0). DRHODAz is scattered as NL x 1 (layer/edge). V_A -> flux.
%   dt_s>0 is required (Density_DA/calculateThermalLiquidDispersion divide by Delt_t).
%   Optional refreshMode is 'full' (default) or 'internal_energy'.  The
%   latter evaluates only the phase-dependent density, wetting-heat, latent-
%   heat, and heat-capacity kernels required by U(h,T,P_g) inversion; it does
%   not update conductivity or dispersivity caches that the state function
%   never reads.
%
%   Reference: docs/superpowers/specs/2026-07-07-stage4-B2-thermal-gas-transport-design.md.
%   Oracle: tests/test_B2_thermal_gas_oracle.m.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-10 Jizixing Zhu -- map legacy h/T/P_g from accepted substep-start Y; workspace.step.previous is hOLD/TOLD/P_gOLD for A3 extrapolation only.
%   Change: 2026-07-13 Jizixing Zhu -- scatter gas pressure mobility through
%           the Kg field without changing its computed values.
%   Change: 2026-07-13 Jizixing Zhu -- require a positive trial time step and
%           honor the configured thermal-vapor diffusion switch.
%   Change: 2026-07-13 Jizixing Zhu -- retain Alpha_Lg in the derived cache
%           so gas dispersivity is observable and calibratable.
%   Change: 2026-07-14 Jizixing Zhu -- retain dW_dh and dW_dTheta so the
%           heat-storage Jacobian receives the complete wetting-heat chain rule.
%   Change: 2026-07-17 Jizixing Zhu -- add an internal-energy-only refresh
%           mode that reuses the production kernels and skips unused transport.
    if nargin < 3 || isempty(refreshMode)
        refreshMode = 'full';
    end
    refreshMode = validatestring(refreshMode, {'full', 'internal_energy'}, ...
        mfilename, 'refreshMode');
    assert(nargin >= 2 && isscalar(dt_s) && isfinite(dt_s) && dt_s > 0, ...
        'refreshThermalGasTransport:invalidTimeStep', ...
        'dt_s must be a finite positive scalar.');
    frame = assembleThermalGasFrame(built, dt_s);
    % --- 9 kernels, COMPASS.m:1685-1694 order (read-only oracles, unchanged) ---
    SV = frame.SoilVariables; TC = frame.TransportCoefficient; cfg = frame.cfg;
    SP = frame.SoilProperties; ThC = frame.ThermalConductivity; VG = frame.VanGenuchten;
    if strcmp(refreshMode, 'full')
        TC = soil.thermal.conductivity.calculateThermalLiquidConductivity(SV, cfg, TC);
    end
    SV = soil.gas.Density_V(SV, cfg);
    TC = soil.thermal.conductivity.calculateThermalLiquidDispersion(SV, SP, VG, cfg, dt_s, TC);
    L  = soil.thermal.Latent(SV, cfg);
    SV = soil.gas.Density_DA(SV, cfg, dt_s);
    TC = soil.thermal.conductivity.calculateThermalConductivityCapacity(SV, SP, ThC, VG, cfg, L, TC);
    if strcmp(refreshMode, 'full')
        TC = soil.gas.conductivity.calculateGasConductivity(SV, SP, VG, TC, cfg);
        TC = soil.gas.conductivity.calculateVaporVariables(SV, SP, VG, TC, cfg);
        TC = soil.gas.conductivity.calculateGasDispersivity(SV, TC, cfg);
        built = scatterThermalGasFrame(built, SV, TC, L);
    else
        built = scatterInternalEnergyFrame(built, SV, TC, L);
    end
end

function frame = assembleThermalGasFrame(built, dt_s)
% Marshal built -> legacy working frame. Only PRIMARY inputs + params + switches; intra-chain
% products (RHOV/RHODA/MU_W/ZETA/k_g/Sa/L) are produced by the kernels. Paths per the B2 plan
% ASSEMBLE SOURCE MAP (file:line-grounded). INDEPENDENT of tests/+wofixture/thermalGasFixture (do
% NOT share an assemble helper -- the oracle relies on the two assemblies being separate).
    NL = built.domain.NL; NN = built.domain.NN;
    tr = built.workspace.step.trial; Y = built.Y;
    d  = built.workspace.derived; S = built.parameters.soil; sw = built.control.switches;

    SV = struct();
    SV.hh = tr.h(:); SV.TT = tr.T(:); SV.h = Y.h(:); SV.T = Y.T(:);   % accepted substep-start h/T (legacy h/T)
    SV.h_frez = tr.h_frez(:); SV.hh_frez = tr.h_frez(:);
    SV.P_gg = tr.p_gas(:); SV.P_g = Y.p_gas(:);
    SV.Theta_LL = d.Theta_LL; SV.Theta_g = d.Theta_g; SV.Theta_UU = d.Theta_UU;
    SV.Theta_V = d.Theta_V; SV.Theta_II = d.Theta_II; SV.Se = d.Se;
    SV.KL_h = d.KL_h; SV.KfL_h = d.KfL_h; SV.Ratio_ice = d.Ratio_ice; SV.Gamma_hh = d.Gamma_hh(:);
    SV.Theta_L = built.Y.Theta_L;                       % prev-step liquid theta (buildInitialState:144)
    SV.POR = S.porosity(:); SV.Imped = S.modes.freeze.imped(:);
    SV.XWILT = S.wilting_point(:); SV.XK = S.xk(:); SV.XSOC = S.msoc_volume(:);
    SV.DVT_Switch = double(isfield(sw, 'lDVT_Switch') && sw.lDVT_Switch == 1);
    if isfield(sw,'lSoilairefc') && sw.lSoilairefc == 1
        SV.KaT_Switch = 1;                              % StartInit:62 (air-only)
    else
        SV.KaT_Switch = [];                             % StartInit:58
    end

    % SoilProperties. SSUR/fc are legacy scalar constants with NO new-container path yet:
    %   source: legacy SoilProperties.m:78-79 and run_stemmus_initialize.m:232-233. built.parameters.soil
    %   currently has no formal field. These are documented constants injected here so the B2 adapter
    %   can fill the legacy working frame the reused kernels expect.
    %   TODO after B2/B3: promote SSUR/fc into built.parameters.soil(.thermal) and read from there.
    SP = struct('Ks', S.ks(:), 'Bulk_density', S.bulk_density(:), 'SSUR', 1e5, 'fc', 0.022);

    ThC = struct('HCAP', S.hcap_phase, 'SF', S.sf_phase, 'TCA', S.tca, ...
        'GA1', S.ga1, 'GA2', S.ga2, 'GB1', S.gb1, 'GB2', S.gb2, 'HCD', S.heat_capacity_dry, ...
        'ZETA0', S.zeta0, 'CON0', S.con0, 'PS1', S.ps1, 'PS2', S.ps2, ...
        'TPS1', S.tps1, 'TPS2', S.tps2, 'TCON0', S.tcon0, 'FEHCAP', S.fehcap, 'TCON_dry', S.tcon_dry);
    ThC.TCON_s = localTconS(S, sw);                     % Johansen vs Farouki per lThermCond (setInitialThermal:52-56)

    VG = struct('Theta_s', S.theta_s(:), 'Theta_r', S.theta_r(:), ...
        'Alpha', S.alpha(:), 'n', S.n(:), 'm', S.m(:));

    cfg = struct();
    cfg.Grid = struct('NL', NL, 'NN', NN, 'nD', built.domain.grid.nD, 'DeltZ', built.domain.dz);
    for f = {'lKLT_Switch','lDVT_Switch','lW_Chg','lThmrlCondCap','lThermCond','lSWCC','lSoilairefc','lThmrlefc','lDEBUG_NAN_TRACE'}
        if isfield(sw, f{1}), cfg.(f{1}) = sw.(f{1}); else, cfg.(f{1}) = 0; end
    end

    frame = struct('SoilVariables', SV, 'TransportCoefficient', struct(), ...
        'SoilProperties', SP, 'ThermalConductivity', ThC, 'VanGenuchten', VG, 'cfg', cfg, 'dt_s', dt_s);
end

function v = localTconS(S, sw)
% TCON_s selection mirrors setInitialThermal.m:52-56 (Farouki method vs Johansen/default).
    method = 0; if isfield(sw,'lThermCond'), method = sw.lThermCond; end
    if method == 4 && isfield(S,'tcon_s_farouki')
        v = S.tcon_s_farouki;                           % Farouki
    else
        v = S.tcon_s;                                   % Johansen / default
    end
end

function built = scatterThermalGasFrame(built, SV, TC, L)
% Scatter produced fields to new-frame targets. derived.* for all except V_A -> flux.V_A.
    d = built.workspace.derived;
    d.KL_T = TC.KL_T;
    d.RHOV = SV.RHOV(:); d.DRHOVh = SV.DRHOVh(:); d.DRHOVT = SV.DRHOVT(:);
    d.W = TC.W; d.WW = TC.WW;
    d.dW_dh = TC.dW_dh; d.dW_dTheta = TC.dW_dTheta;
    d.MU_W = TC.MU_W; d.D_Ta = TC.D_Ta;
    d.L = L(:);
    d.RHODA = SV.RHODA(:); d.Xaa = SV.Xaa(:); d.XaT = SV.XaT(:); d.Xah = SV.Xah(:); d.DRHODAt = SV.DRHODAt(:);
    d.DRHODAz = SV.DRHODAz(:);                          % legacy 1xNL -> NLx1
    d.EHCAP = TC.EHCAP; d.c_unsat = TC.c_unsat; d.Lambda_eff = TC.Lambda_eff; d.ZETA = TC.ZETA; d.ETCON = TC.ETCON;
    d.Sa = TC.Sa; d.k_g = TC.k_g;
    d.D_V = TC.D_V; d.Eta = TC.Eta; d.D_A = TC.D_A(:);
    d.Alpha_Lg = TC.Alpha_Lg; d.D_Vg = TC.D_Vg; d.Kg = TC.Kg;
    built.workspace.derived = d;
    built.workspace.flux.V_A = TC.V_A;                  % V_A -> flux, NOT derived
end

function built = scatterInternalEnergyFrame(built, SV, TC, L)
% Scatter only auxiliary fields read by the U(h,T,P_g) state function.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    d = built.workspace.derived;
    d.RHOV = SV.RHOV(:);
    d.DRHOVh = SV.DRHOVh(:);
    d.DRHOVT = SV.DRHOVT(:);
    d.RHODA = SV.RHODA(:);
    d.Xaa = SV.Xaa(:);
    d.XaT = SV.XaT(:);
    d.Xah = SV.Xah(:);
    d.L = L(:);
    d.WW = TC.WW;
    d.dW_dh = TC.dW_dh;
    d.dW_dTheta = TC.dW_dTheta;
    d.c_unsat = TC.c_unsat;
    built.workspace.derived = d;
end
