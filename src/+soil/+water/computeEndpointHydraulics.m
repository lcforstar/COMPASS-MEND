function endpoint = computeEndpointHydraulics(soilParams, switches, state, L_f, phys)
%COMPUTEENDPOINTHYDRAULICS  Vectorized soil-water endpoint hydraulic kernel.
%
%   endpoint = soil.water.computeEndpointHydraulics(soilParams, switches,
%   state, L_f, phys) computes the layer-endpoint auxiliary fields used by
%   UpdateSoilWaterContent:
%       Gamma_hh -> Theta_m -> Theta_UU -> Theta_II -> Theta_LL ->
%       Dall'Amico phase relaxation -> DTheta -> Se -> KL_h/KfL_h.
%
%   Inputs:
%     soilParams [in] struct built.parameters.soil, bottom-up [NL x 1]
%     switches   [in] struct built.control.switches
%     state      [in] struct endpoint arrays, each [NL x nD]:
%                    .h [cm], .T_C [degC], .h_frez [cm],
%                    .TT_CRIT [degC], .theta_i_prev [-],
%                    optional .dt_s [s], .EHCAP, .c_unsat, .Lambda_eff,
%                    .dz_cm [cm]
%     L_f        [in] [J/kg] latent heat coefficient. L_f==0 disables
%                    ice/latent-heat effects inside this same formula path.
%     phys       [in] struct physical constants from equations.physicalConstants()
%
%   Output:
%     endpoint [out] struct of [NL x nD] auxiliary arrays:
%                    theta_l/theta_u/theta_i, theta_v, Se_air, KL_h,
%                    KfL_h, DTheta_LLh, DTheta_UUh, ratio_ice,
%                    gamma_hh, and Dall'Amico diagnostics.
%
%   Role:
%     This is the canonical vectorized replacement for the water-content and
%     hydraulic-conductivity auxiliary part of UpdateSoilWaterContent. Frozen
%     and unfrozen cases use one formula path; the only numerical difference is
%     L_f. When L_f==0, ice content and ice impedance are zero by the same
%     equations rather than by a separate method.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: Replaced scalar endpoint helper with a vectorized kernel.
%   Change: 2026-06-28 Jizixing Zhu — KfL_T uses heaviside1 half-value at the
%           freezing threshold; lSFCC=0 path restores the legacy case-1 DTheta_UUh.
%   Change: 2026-07-11 Jizixing Zhu -- add positive, semantically named local
%           retention derivatives while preserving legacy DTheta outputs.
%   Change: 2026-07-12 Jizixing Zhu -- add the ClimaLand-style variably saturated
%           branch: augmented water storage is theta_s + Ss*h for h >= 0,
%           its pressure-head capacity is Ss, and conductivity reaches Ks.
%   Change: 2026-07-12 Jizixing Zhu -- keep Ss*h as total elastic storage while
%           capping physical liquid-ice pore water at theta_s.
%   Change: 2026-07-14 Jizixing Zhu -- remove stale KfL_T claims from the
%           endpoint contract; the current strict-SFCC kernel does not return it.
%   Change: 2026-07-14 Jizixing Zhu -- remove the unused TT_CRIT read and an
%           identical conditional left behind by the retired coefficient path.

    if nargin < 5 || isempty(phys)
        phys = equations.physicalConstants();
    end
    if nargin < 4 || isempty(L_f)
        L_f = 0;
    end

    h = state.h;
    T_C = state.T_C;
    hFrez = state.h_frez;
    thetaIPrev = readStateArray(state, 'theta_i_prev', zeros(size(h)));
    p = endpointSoilParameters(soilParams, size(h, 1), size(h, 2));
    LfActive = (L_f ~= 0);

    gamma = calculateGammaArray(h, p);
    thetaM = calculateThetaMArray(gamma, p);
    thetaU = calculateThetaUUArray(thetaM, gamma, h, p, switches, LfActive);
    [thetaI, dThetaIExplicitT] = calculateThetaIIArray(thetaIPrev, T_C, h, p, switches, LfActive);
    thetaL = calculateThetaLLArray(thetaU, thetaI, thetaM, gamma, h, hFrez, T_C, p, switches, LfActive, phys);
    [thetaU, thetaI, thetaL] = closeIceMassBalanceArray( ...
        thetaU, thetaL, p, LfActive, phys);

    dall = emptyDallAmico(size(h));
    if shouldApplyDallAmico(switches, state, L_f)
        rhoC = readStateArray(state, 'EHCAP', NaN(size(h)));
        cUnsat = readStateArray(state, 'c_unsat', NaN(size(h)));
        bad = ~isfinite(rhoC) | rhoC <= 0;
        rhoC(bad) = cUnsat(bad);
        phase = soil.water.conductivity.calculateDallAmicoPhaseChangeArray( ...
            thetaL, thetaI, T_C, p.theta_s, p.theta_r, p.alpha, p.n, p.m, ...
            rhoC, state.Lambda_eff, state.dz_cm, state.dt_s);
        compressedWater = max(thetaU - p.theta_s, 0);
        thetaL = phase.theta_l_next;
        thetaI = phase.theta_i_next;
        thetaU = thetaL + thetaI * phys.RHOI / phys.RHOL + compressedWater;
        dall.thermal_time_s = phase.thermal_time_s;
        dall.theta_l_star = phase.theta_l_star;
        dall.source_liquid_s = phase.source_liquid_s;
        dall.limiter_code = phase.limiter_code;
    end

    [dThetaUdh, dThetaLdh] = calculateDerivativeArrays(thetaL, thetaM, gamma, h, hFrez, p, switches, LfActive);
    physicalPartials = calculatePhysicalDerivativeArrays( ...
        thetaL, thetaM, gamma, h, hFrez, p, switches, LfActive);
    Se = calculateSeArray(thetaL, h, hFrez, p, switches, LfActive);
    KL = calculateKLArray(thetaL, T_C, p, phys);
    ratioIce = calculateRatioIceArray(thetaU, thetaI, phys);
    KfL = KL .* 10.^(-p.imped .* ratioIce);
    thetaV = max(p.porosity - thetaL - thetaI, 1e-14);

    endpoint = struct();
    endpoint.theta_l = thetaL;
    endpoint.theta_u = thetaU;
    endpoint.theta_i = thetaI;
    endpoint.theta_v = thetaV;
    endpoint.Se_air = Se;
    endpoint.KL_h = KL;
    endpoint.KfL_h = KfL;
    endpoint.dtheta_l_dh = dThetaLdh;
    endpoint.dtheta_u_dh = dThetaUdh;
    endpoint.dtheta_water_equiv_dbase_head = physicalPartials.dtheta_u_dbase;
    endpoint.dtheta_liquid_dbase_head = physicalPartials.dtheta_l_dbase;
    endpoint.dtheta_liquid_dfreeze_head = physicalPartials.dtheta_l_dfreeze;
    endpoint.dtheta_ice_explicit_dT = dThetaIExplicitT;
    endpoint.ratio_ice = ratioIce;
    endpoint.gamma_hh = gamma;
    endpoint.dallamico_thermal_time = dall.thermal_time_s;
    endpoint.dallamico_theta_l_star = dall.theta_l_star;
    endpoint.dallamico_phase_source = dall.source_liquid_s;
    endpoint.dallamico_limiter_code = dall.limiter_code;
end

function p = endpointSoilParameters(soilParams, NL, nD)
% Broadcast layer parameters to endpoint arrays.
    oneRow = ones(1, nD);
    p = struct();
    p.theta_s = soilParams.theta_s(:) * oneRow;
    p.theta_r = soilParams.theta_r(:) * oneRow;
    p.alpha = soilParams.alpha(:) * oneRow;
    p.n = soilParams.n(:) * oneRow;
    p.m = soilParams.m(:) * oneRow;
    p.porosity = soilParams.porosity(:) * oneRow;
    p.field_capacity = soilParams.field_capacity(:) * oneRow;
    p.ks = soilParams.ks(1:NL) * oneRow;
    p.imped = soilParams.modes.freeze.imped(1:NL) * oneRow;
    p.Ss = soilParams.modes.retention.Ss .* ones(NL, nD);
    p.h_trans = soilParams.modes.retention.h_trans;
    p.hd = soilParams.constants.hd;
    p.hm = soilParams.constants.hm;
end

function gamma = calculateGammaArray(h, p)
% Vectorized Gamma_hh.
    gamma = ones(size(h));
    absh = abs(h);
    gamma(absh >= abs(p.hd)) = 0;
    mask = (absh >= abs(p.hm)) & (absh < abs(p.hd));
    gamma(mask) = log(abs(p.hd) ./ absh(mask)) ./ log(abs(p.hd) / abs(p.hm));
end

function thetaM = calculateThetaMArray(gamma, p)
% Vectorized Theta_m capillary upper bound.
    thetaM = gamma .* p.theta_r + (p.theta_s - gamma .* p.theta_r) .* ...
        (1 + abs(p.alpha .* (-1)).^p.n).^(-p.m);
    hi = thetaM >= p.porosity;
    lo = ~hi & thetaM <= p.theta_s;
    thetaM(hi) = p.porosity(hi);
    thetaM(lo) = p.theta_s(lo);
end

function thetaU = calculateThetaUUArray(thetaM, gamma, h, p, switches, LfActive)
% Total water content before liquid-ice closure.
    if ~isSwitchOn(switches, 'lSWCC')
        thetaU = calculateClappThetaArray(h, p);
        return;
    end

    if LfActive && isSwitchOn(switches, 'lSFCC')
        sat = h >= 0;
        dry = ~sat & h <= p.hd;
        zero = ~sat & ~dry & gamma == 0;
        vg = ~sat & ~dry & gamma ~= 0;
        thetaU = zeros(size(h));
        % ClimaLand soil_hydrology_parameterizations.jl:121-125 inverts to
        % augmented theta = theta_s + Ss*h on the positive-pressure branch.
        thetaU(sat) = p.theta_s(sat) + p.Ss(sat) .* h(sat);
        thetaU(dry) = p.theta_r(dry);
        thetaU(zero) = 0;
        thetaVg = calculateThetaCase0Array(thetaM, h, gamma, p);
        thetaU(vg) = thetaVg(vg);
    else
        thetaU = calculateThetaCase2Array(h, p);
        sat = h >= 0;
        % ClimaLand soil_hydrology_parameterizations.jl:121-125.
        thetaU(sat) = p.theta_s(sat) + p.Ss(sat) .* h(sat);
        thetaU(h <= -1e7) = p.theta_r(h <= -1e7);
    end
end

function [thetaI, dThetaIExplicitT] = calculateThetaIIArray(thetaIPrev, T_C, h, p, switches, LfActive)
% Ice content pre-estimate.
    dThetaIExplicitT = zeros(size(h));
    if ~LfActive
        thetaI = zeros(size(h));
        return;
    end

    thetaI = thetaIPrev;
    thetaI(h <= -1e7) = 0;
    if isSwitchOn(switches, 'lSWCC') && ~isSwitchOn(switches, 'lSFCC')
        Tf1 = 273.15 + 1;
        Tf2 = 273.15 - 3;
        TK = T_C + 273.15;
        thetaI(TK > Tf1) = 0;
        mid = TK <= Tf1 & TK >= Tf2;
        phase = pi() * (TK(mid) - 0.5 * Tf1 - 0.5 * Tf2) / (Tf1 - Tf2);
        frac = 0.5 * (1 - sin(phase));
        thetaI(mid) = max(0, min(1, frac)) .* p.field_capacity(mid);
        dThetaIExplicitT(mid) = -0.5 * cos(phase) * pi() / (Tf1 - Tf2) .* ...
            p.field_capacity(mid);
        thetaI(TK < Tf2) = p.field_capacity(TK < Tf2);
    end
end

function thetaL = calculateThetaLLArray(thetaU, thetaI, thetaM, gamma, h, hFrez, T_C, p, switches, LfActive, phys)
% Liquid water content before liquid-ice mass closure.
    if ~isSwitchOn(switches, 'lSWCC')
        thetaL = calculateClappThetaArray(h, p);
        return;
    end

    if LfActive && isSwitchOn(switches, 'lSFCC')
        heatTerm = h + hFrez;
        allLiquid = heatTerm >= 0;
        dry = ~allLiquid & heatTerm <= p.hd;
        vg = ~allLiquid & ~dry;
        thetaL = zeros(size(h));
        thetaL(allLiquid) = min(thetaU(allLiquid), p.theta_s(allLiquid));
        thetaL(dry) = p.theta_r(dry);
        thetaVg = calculateThetaCase0Array(thetaM, heatTerm, gamma, p);
        thetaL(vg) = thetaVg(vg);
    elseif LfActive
        Tf1 = 273.15 + 1;
        heatTerm = h + hFrez;
        phaseWater = min(thetaU, p.theta_s);
        thetaL = phaseWater - thetaI * phys.RHOI / phys.RHOL;
        sat = h >= 0 & T_C + 273.15 > Tf1;
        dry = h <= -1e7;
        warm = ~sat & ~dry & T_C + 273.15 > Tf1;
        thetaL(sat) = p.theta_s(sat);
        thetaL(dry) = p.theta_r(dry);
        thetaWarm = calculateThetaCase2Array(heatTerm, p);
        thetaL(warm) = thetaWarm(warm);
        thetaL(thetaL <= 0.06) = 0.06;
    else
        thetaL = min(thetaU, p.theta_s);
    end
end

function [thetaU, thetaI, thetaL] = closeIceMassBalanceArray(thetaU, thetaL, p, LfActive, phys)
% Close pore-water phases without converting elastic storage into ice.
    phaseWater = max(0, min(thetaU, p.theta_s));
    if ~LfActive
        thetaL = phaseWater;
        thetaI = zeros(size(thetaL));
    else
        thetaL = max(0, min(thetaL, phaseWater));
        thetaI = max(0, (phaseWater - thetaL) * phys.RHOL / phys.RHOI);
    end
end

function [dThetaUdh, dThetaLdh] = calculateDerivativeArrays(thetaL, thetaM, gamma, h, hFrez, p, switches, LfActive)
% Hydraulic storage derivatives.
    if ~isSwitchOn(switches, 'lSWCC')
        dThetaUdh = calculateClappDerivativeArray(h, p);
        dThetaLdh = calculateClappDerivativeArray(h + hFrez, p);
        return;
    end

    if LfActive && isSwitchOn(switches, 'lSFCC') && isSwitchOn(switches, 'lThmrlefc')
        heatTerm = h + hFrez;
        dThetaUdh = calculateDThetaCase0Array(heatTerm, thetaM, gamma, p);
        dThetaUdh((heatTerm <= p.hd) | ((h >= -1) & (hFrez >= 0))) = 0;
        dThetaLdh = calculateDThetaCase0Array(h, thetaM, gamma, p);
        dThetaLdh((h >= -1) | (heatTerm <= p.hd) | (gamma == 0)) = 0;
    elseif LfActive && ~isSwitchOn(switches, 'lSFCC')
        % lSFCC=0 simplified-freeze path: keep the legacy case-1 analytic
        % DTheta_UUh (theta_m, heatTerm = h + h_frez), matching the scalar
        % calculateDTheta_UUh/calculateDTheta(1,...) under lThmrlefc=1.
        heatTerm = h + hFrez;
        dThetaUdh = calculateDThetaCase1Array(heatTerm, thetaM, p);
        dThetaUdh(h >= -1e-6 | h <= -1e7) = 0;
        dThetaLdh = dThetaUdh;
        dThetaLdh(h <= -1e7 | thetaL <= 0.06) = 0;
    else
        % L_f==0 (no freezing) stays on the standard van Genuchten derivative.
        dThetaUdh = calculateDThetaCase2Array(h, p);
        dThetaUdh(h >= -1e-6 | h <= -1e7) = 0;
        dThetaLdh = dThetaUdh;
    end
end

function partials = calculatePhysicalDerivativeArrays(thetaL, thetaM, gamma, h, hFrez, p, switches, LfActive)
% Return positive local derivatives with explicit mathematical meanings.
% These fields are independent of the inherited DTheta_LLh/DTheta_UUh names.
    if ~isSwitchOn(switches, 'lSWCC')
        c = calculateClappDerivativeArray(h, p);
        partials = packPhysicalPartials(c, c, zeros(size(c)));
        return;
    end

    if LfActive && isSwitchOn(switches, 'lSFCC')
        heatTerm = h + hFrez;
        [cUbase, ~] = strictRetentionPartials(h, h, thetaM, gamma, p);
        [cLbase, cLfreeze] = strictRetentionPartials(h, heatTerm, thetaM, gamma, p);
        saturated = h >= 0;
        allLiquid = heatTerm >= 0;
        cUbase(h <= p.hd) = 0;
        % ClimaLand soil_hydrology_parameterizations.jl:149-150.
        cUbase(saturated) = p.Ss(saturated);
        cLbase(heatTerm <= p.hd) = 0;
        cLfreeze(heatTerm <= p.hd) = 0;
        cLbase(allLiquid) = cUbase(allLiquid);
        cLbase(allLiquid & saturated) = 0;
        cLfreeze(allLiquid) = 0;
        partials = packPhysicalPartials(cUbase, cLbase, cLfreeze);
        return;
    end

    cUbase = positiveVgDerivative(h, p.theta_s, p.theta_r, p.alpha, p.n, p.m);
    saturated = h >= 0;
    cUbase(h <= -1e7) = 0;
    % ClimaLand soil_hydrology_parameterizations.jl:149-150 gives
    % dpsi/dtheta=1/Ss, hence dtheta/dpsi=Ss in a head-primary solve.
    cUbase(saturated) = p.Ss(saturated);
    % The simplified-freeze branch partitions theta_u with an explicit
    % sinusoidal ice fraction; its head derivative therefore equals the
    % total-water derivative. Temperature partitioning is added later.
    cLbase = cUbase;
    cLbase(saturated) = 0;
    cLbase(thetaL <= 0.06 & LfActive) = 0;
    partials = packPhysicalPartials(cUbase, cLbase, zeros(size(cLbase)));
end

function [dThetaDbase, dThetaDfreeze] = strictRetentionPartials(baseHead, heatTerm, thetaM, gamma, p)
% Differentiate theta(gamma(baseHead), heatTerm) by the chain rule.
    pow = 1 + abs(p.alpha .* heatTerm).^p.n;
    slope = (thetaM - gamma .* p.theta_r) .* p.m .* p.n .* p.alpha .* ...
        abs(p.alpha .* heatTerm).^(p.n - 1) .* pow.^(-p.m - 1);
    slope(~isfinite(slope)) = 0;

    dGammaDbase = zeros(size(baseHead));
    transition = abs(baseHead) >= abs(p.hm) & abs(baseHead) < abs(p.hd);
    dGammaDbase(transition) = -1 ./ (baseHead(transition) .* log(abs(p.hd / p.hm)));

    g1 = (1 + abs(p.alpha .* (-1)).^p.n).^(-p.m);
    thetaMraw = gamma .* p.theta_r + (p.theta_s - gamma .* p.theta_r) .* g1;
    dThetaMDgamma = p.theta_r .* (1 - g1);
    clamped = thetaMraw >= p.porosity | thetaMraw <= p.theta_s;
    dThetaMDgamma(clamped) = 0;
    G = pow.^(-p.m);
    gammaTerm = dGammaDbase .* (p.theta_r + (dThetaMDgamma - p.theta_r) .* G);

    dThetaDbase = gammaTerm + slope;
    dThetaDfreeze = slope;
    dThetaDbase = max(dThetaDbase, 0);
    dThetaDfreeze = max(dThetaDfreeze, 0);
end

function dtheta = positiveVgDerivative(h, thetaS, thetaR, alpha, n, m)
% van Genuchten 1980 SSSAJ 44:892 Eq. 21, differentiated for signed h < 0.
    dtheta = zeros(size(h));
    mask = h < 0;
    x = abs(alpha(mask) .* h(mask));
    dtheta(mask) = (thetaS(mask) - thetaR(mask)) .* m(mask) .* n(mask) .* ...
        alpha(mask) .* x.^(n(mask) - 1) .* (1 + x.^n(mask)).^(-m(mask) - 1);
    dtheta = max(dtheta, 0);
end

function partials = packPhysicalPartials(cUbase, cLbase, cLfreeze)
% Package local phase-retention derivatives before raw h/T chain rules.
    partials = struct('dtheta_u_dbase', cUbase, ...
        'dtheta_l_dbase', cLbase, 'dtheta_l_dfreeze', cLfreeze);
end

function Se = calculateSeArray(thetaL, h, hFrez, p, switches, LfActive)
% Effective saturation used by gas and vapor auxiliary fields.
    Se = thetaL ./ p.porosity;
    if isSwitchOn(switches, 'lSWCC')
        if LfActive && isSwitchOn(switches, 'lSFCC')
            heatTerm = h + hFrez;
            sat = h >= -1 & hFrez >= 0;
            dry = ~sat & heatTerm <= p.hd;
            Se(sat) = 1;
            Se(dry) = 0;
        else
            Se(thetaL <= 0.06 | h <= -1e7) = 0;
        end
    else
        Se(h + hFrez <= -1e7 | h <= -1e7) = 0;
    end
    Se = max(0, min(1, Se));
end

function KL = calculateKLArray(thetaL, T_C, p, phys)
% Mualem-VG liquid conductivity with temperature-viscosity correction.
% NOTE (Jizixing Zhu, 2026-07-03): phys.R (8.314468) here is the canonical gas constant; the legacy
% kernel calculateHydraulicConductivityArray.m hardcodes the STEMMUS literal 8.31441 in the same
% MU_WN/MU_W term. That ~7e-6-rel constant gap induces a benign ~ppm KfL_h difference (0 at the 20C
% reference; -3.4e-6 rel at 2C). Deliberately NOT changed to match the legacy literal (decision
% 2026-07-03): the new chain keeps the correct phys.R. Verified in tests/test_runStepEbal_parity.m.
    MU_WN = phys.MU_W0 * exp(phys.MU1 / (phys.R * (phys.Tr + 133.3)));
    MU_W = phys.MU_W0 .* exp(phys.MU1 ./ (phys.R * (T_C + 133.3)));
    MU_W(T_C < -20) = 3.71e-2;
    MU_W(T_C > 150) = 1.81e-3;
    SeK = (thetaL - p.theta_r) ./ max(p.theta_s - p.theta_r, 1e-12);
    SeK = max(0, min(1, SeK));
    KL = (MU_WN ./ MU_W) .* p.ks .* sqrt(SeK) .* (1 - (1 - SeK.^(1 ./ p.m)).^p.m).^2;
    KL(KL <= 1e-20) = 1e-20;
    KL(isnan(KL)) = 0;
    KL(SeK == 0) = 0;
end

function ratioIce = calculateRatioIceArray(thetaU, thetaI, phys)
% Ice-to-water ratio used by the impedance correction.
    ratioIce = zeros(size(thetaU));
    nz = thetaU ~= 0;
    ratioIce(nz) = phys.RHOI * thetaI(nz) ./ (phys.RHOL * thetaU(nz));
end

function theta = calculateThetaCase0Array(thetaM, heatTerm, gamma, p)
% Strict SFCC retention curve.
    theta = gamma .* p.theta_r + (thetaM - gamma .* p.theta_r) ./ ...
        (1 + abs(p.alpha .* heatTerm).^p.n).^p.m;
    theta(heatTerm <= p.hd) = 0;
end

function theta = calculateThetaCase2Array(h, p)
% Standard van Genuchten retention curve.
    theta = p.theta_r + (p.theta_s - p.theta_r) ./ (1 + abs(p.alpha .* h).^p.n).^p.m;
end

function dtheta = calculateDThetaCase0Array(heatTerm, thetaM, gamma, p)
% Strict SFCC derivative with the legacy solver sign convention.
    pow = 1 + abs(p.alpha .* heatTerm).^p.n;
    A = (-p.theta_r) ./ (abs(heatTerm) .* log(abs(p.hd / p.hm))) .* (1 - pow.^(-p.m));
    B = p.alpha .* p.n .* p.m .* (thetaM - gamma .* p.theta_r);
    C = pow.^(-p.m - 1);
    D = abs(p.alpha .* heatTerm).^(p.n - 1);
    dtheta = A - B .* C .* D;
    dtheta(heatTerm <= p.hd) = 0;
end

function dtheta = calculateDThetaCase2Array(h, p)
% Standard VG derivative with the legacy solver sign convention.
    dtheta = zeros(size(h));
    mask = h < -1e-6 & h > -1e7;
    x = abs(p.alpha(mask) .* h(mask));
    positive = (p.theta_s(mask) - p.theta_r(mask)) .* p.m(mask) .* p.n(mask) .* ...
        p.alpha(mask) .* x.^(p.n(mask) - 1) .* (1 + x.^p.n(mask)).^(-p.m(mask) - 1);
    dtheta(mask) = -max(positive, 0);
end

function dtheta = calculateDThetaCase1Array(heatTerm, thetaM, p)
% Legacy case-1 VG derivative (theta_m form) with the solver sign convention.
    pow = 1 + abs(p.alpha .* heatTerm).^p.n;
    A = (thetaM - p.theta_r) .* p.alpha .* p.n;
    B = abs(p.alpha .* heatTerm).^(p.n - 1) .* (-p.m);
    C = pow.^(-p.m - 1);
    dtheta = A .* B .* C;
end

function tf = shouldApplyDallAmico(switches, state, L_f)
% True when the vectorized kernel should apply Dall'Amico relaxation.
    tf = (L_f ~= 0) ...
        && isfield(state, 'dt_s') && state.dt_s > 0 ...
        && isSwitchOn(switches, 'lThmrlefc') ...
        && getLogicalSwitch(switches, 'lDallAmicoPhaseChange', false) ...
        && isfield(state, 'Lambda_eff') && any(isfinite(state.Lambda_eff(:)) & state.Lambda_eff(:) > 0) ...
        && isfield(state, 'dz_cm') && any(isfinite(state.dz_cm(:)) & state.dz_cm(:) > 0) ...
        && hasHeatCapacity(state);
end

function tf = hasHeatCapacity(state)
% True if either heat-capacity endpoint field can support thermal-time relaxation.
    tf = false;
    if isfield(state, 'EHCAP')
        tf = tf || any(isfinite(state.EHCAP(:)) & state.EHCAP(:) > 0);
    end
    if isfield(state, 'c_unsat')
        tf = tf || any(isfinite(state.c_unsat(:)) & state.c_unsat(:) > 0);
    end
end

function dall = emptyDallAmico(sz)
% Zero Dall'Amico diagnostic arrays for inactive phase relaxation.
    dall = struct();
    dall.thermal_time_s = zeros(sz);
    dall.theta_l_star = zeros(sz);
    dall.source_liquid_s = zeros(sz);
    dall.limiter_code = zeros(sz);
end

function value = readStateArray(state, name, defaultValue)
% Read optional endpoint state array.
    if isfield(state, name) && ~isempty(state.(name))
        value = state.(name);
    else
        value = defaultValue;
    end
end

function theta = calculateClappThetaArray(~, ~)
% Reserved Clapp-Hornberger interface.
    theta = NaN; %#ok<NASGU>
    error('computeEndpointHydraulics:unsupportedClapp', ...
        'lSWCC=0 requires canonical Clapp-Hornberger parameters before this branch can run.');
end

function dtheta = calculateClappDerivativeArray(~, ~)
% Reserved Clapp-Hornberger derivative interface.
    dtheta = NaN; %#ok<NASGU>
    error('computeEndpointHydraulics:unsupportedClapp', ...
        'lSWCC=0 requires canonical Clapp-Hornberger parameters before this branch can run.');
end

function tf = isSwitchOn(switches, name)
% Read a legacy logical switch from a struct.
    tf = isfield(switches, name) && switches.(name) == 1;
end

function tf = getLogicalSwitch(switches, name, defaultValue)
% Read an optional logical switch with a default.
    if isfield(switches, name)
        tf = logical(switches.(name));
    else
        tf = defaultValue;
    end
end
