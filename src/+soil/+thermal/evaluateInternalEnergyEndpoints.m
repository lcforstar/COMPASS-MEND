function [U, CTh, CTT, CTa] = evaluateInternalEnergyEndpoints(estate, ederiv, ctl)
%EVALUATEINTERNALENERGYENDPOINTS  Unified internal-energy state function.
%
%   [U, CTh, CTT, CTa] = soil.thermal.evaluateInternalEnergyEndpoints(...
%       estate, ederiv, ctl)
%   returns the endpoint volumetric internal energy U [J/cm³] and its exact
%   partial derivatives CTh = ∂U/∂h, CTT = ∂U/∂T, CTa = ∂U/∂P_g at fixed
%   phase composition, using the same multiphase energy definition as the
%   heat-equation time residual (§4.3.1).  All endpoint arrays are NL×2
%   (bottom-up element indexing: DeltZ(1) = deepest layer).
%
%   This is the SINGLE thermodynamic kernel.  Every consumer — Picard
%   coefficient refresh, time residual, column-transform parameters, and
%   temperature inversion — MUST call this function.  No copy of the
%   internal-energy formula may exist outside this file.
%
%   Input structs:
%     estate — per-endpoint thermodynamic state
%       .T            [°C]             NL×2 endpoint temperature
%       .Theta_LL     [cm³/cm³]        NL×2 liquid water content
%       .Theta_II     [cm³/cm³]        NL×2 ice content
%       .Theta_g      [cm³/cm³]        NL×2 gas-filled pore fraction
%       .RHOV         [g/cm³]          NN×1 nodal vapor density
%       .RHODA        [g/cm³]          NN×1 nodal dry-air density
%       .L_vap        [J/g]            NN×1 nodal vaporization energy
%       .WW           [J/g]            NL×2 wetting heat
%       .c_unsat      [J/cm³/K]        NL×2 bulk sensible heat capacity
%
%     ederiv — phase-amount and gas-density derivatives
%       .DThetaLiquidDh, .DThetaLiquidDT   [1/cm], [1/K]        NL×2
%       .DThetaIceDh,    .DThetaIceDT      [1/cm], [1/K]        NL×2
%       .DThetaGasDh,    .DThetaGasDT      [1/cm], [1/K]        NL×2
%       .DRHOVh, .DRHOVT       [g/cm⁴], [g/cm³/K]               NN×1
%       .Xah,    .XaT,  .Xaa   [g/cm⁴], [g/cm³/K], [s²/cm²]     NN×1
%       .dW_dh,  .dW_dTheta    [J/g/cm], [J/g]                  NL×2
%
%     ctl — runtime controls
%       .L_f          [J/kg]    per-step latent heat of fusion
%                               (lambda_fus when freezing, 0 otherwise)
%       .validate     logical   enable full input validation (default true)
%
%   Returns:
%     U     [J/cm³]       endpoint volumetric internal energy (NL×2)
%     CTh   [J·cm⁻³·cm⁻¹] ∂U/∂h  at fixed T,P_g (NL×2)
%     CTT   [J·cm⁻³·K⁻¹]  ∂U/∂T  at fixed h,P_g (NL×2)
%     CTa   [J·cm⁻³·(g⁻¹·cm·s²)]  ∂U/∂P_g at fixed h,T (NL×2)
%           P_g internal coordinate is g·cm⁻¹·s⁻² (≈ 0.1 Pa).
%
%   Reference: docs/compass_soil_water_heat_physics.md §2.6 and §4.3.1.
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%
%   See also: soil.thermal.assembleNodalInternalEnergy,
%             soil.thermal.invertTemperatureFromInternalEnergy

    arguments
        estate  struct
        ederiv  struct
        ctl     struct
    end

    % --- physical constants (fixed) ---
    phys = equations.physicalConstants();
    c_a  = phys.c_a;          % [J/g/K] specific heat of dry air
    c_L  = phys.c_L;          % [J/g/K] specific heat of liquid water
    c_i  = 2.0455;            % [J/g/K] specific heat of ice
    RHOL = phys.RHOL;         % [g/cm³] density of liquid water
    RHOI = phys.RHOI;         % [g/cm³] density of ice
    L0   = phys.L0;           % [J/g]  latent heat of vaporization at 0°C

    % --- runtime control (caller-supplied) ---
    L_f = ctl.L_f;            % [J/kg] per-step latent heat of fusion

    T           = estate.T;
    Theta_LL    = estate.Theta_LL;
    Theta_II    = estate.Theta_II;
    Theta_g     = estate.Theta_g;
    RHOV        = estate.RHOV;
    RHODA       = estate.RHODA;
    L_vap       = estate.L_vap;
    WW          = estate.WW;
    c_unsat     = estate.c_unsat;

    DThetaLiquidDh = ederiv.DThetaLiquidDh;
    DThetaLiquidDT = ederiv.DThetaLiquidDT;
    DThetaIceDh    = ederiv.DThetaIceDh;
    DThetaIceDT    = ederiv.DThetaIceDT;
    DThetaGasDh    = ederiv.DThetaGasDh;
    DThetaGasDT    = ederiv.DThetaGasDT;
    DRHOVh   = ederiv.DRHOVh;
    DRHOVT   = ederiv.DRHOVT;
    Xah      = ederiv.Xah;
    XaT      = ederiv.XaT;
    Xaa      = ederiv.Xaa;
    dW_dh     = ederiv.dW_dh;
    dW_dTheta = ederiv.dW_dTheta;

    doValidate = true;
    if isfield(ctl, 'validate'), doValidate = logical(ctl.validate); end

    [NL, nD] = size(T);
    if doValidate
        validateInternalEnergyInputs();
    end

    U   = zeros(NL, nD);
    CTh = zeros(NL, nD);
    CTT = zeros(NL, nD);
    CTa = zeros(NL, nD);

    for ML = 1:NL
        for ND = 1:nD
            MN = ML + ND - 1;

            tt    = T(ML, ND);
            lv    = L_vap(MN);
            rhov  = RHOV(MN);
            rhoda = RHODA(MN);

            % Per-phase specific energies [J/g] and volumetric energies [J/cm³].
            % Reference state: liquid water at 0°C.
            vaporEnergyPerMass = c_L * tt + lv;
            airEnergyPerMass   = c_a * tt;
            liquidEnergyPerVol = RHOL * (c_L * tt - WW(ML,ND));
            iceEnergyPerVol    = RHOI * (c_i * tt - L_f * 1e-3);
            gasEnergyPerVol    = rhov * vaporEnergyPerMass ...
                               + rhoda * airEnergyPerMass;

            % ---- internal energy U (exact state function) ----
            U(ML,ND) = c_unsat(ML,ND) * tt ...
                     + L0 * Theta_g(ML,ND) * rhov ...
                     - RHOI * (L_f * 1e-3) * Theta_II(ML,ND) ...
                     - RHOL * WW(ML,ND) * Theta_LL(ML,ND);

            % ---- CTh = ∂U/∂h|_{T,P_g} ----
            phaseStorageDh = liquidEnergyPerVol * DThetaLiquidDh(ML,ND) ...
                           + iceEnergyPerVol    * DThetaIceDh(ML,ND) ...
                           + gasEnergyPerVol    * DThetaGasDh(ML,ND);
            vaporDensityStorageDh = Theta_g(ML,ND) * vaporEnergyPerMass * DRHOVh(MN);
            airDensityStorageDh   = Theta_g(ML,ND) * airEnergyPerMass   * Xah(MN);
            CTh(ML,ND) = phaseStorageDh + vaporDensityStorageDh + airDensityStorageDh;
            % Wetting-heat chain rule: -ρ_l θ_l ∂W/∂h  (§4.3.1d, Eq. 4.3.4a)
            CTh(ML,ND) = CTh(ML,ND) - RHOL * Theta_LL(ML,ND) ...
                * (dW_dh(ML,ND) + dW_dTheta(ML,ND) * DThetaLiquidDh(ML,ND));

            % ---- CTT = ∂U/∂T|_{h,P_g} ----
            phaseStorageDT = liquidEnergyPerVol * DThetaLiquidDT(ML,ND) ...
                           + iceEnergyPerVol    * DThetaIceDT(ML,ND) ...
                           + gasEnergyPerVol    * DThetaGasDT(ML,ND);
            vaporDensityStorageDT = Theta_g(ML,ND) * vaporEnergyPerMass * DRHOVT(MN);
            airDensityStorageDT   = Theta_g(ML,ND) * airEnergyPerMass   * XaT(MN);
            CTT(ML,ND) = c_unsat(ML,ND) + phaseStorageDT ...
                       + vaporDensityStorageDT + airDensityStorageDT;

            % ---- CTa = ∂U/∂P_g|_{h,T} ----
            CTa(ML,ND) = Theta_g(ML,ND) * airEnergyPerMass * Xaa(MN);
        end
    end

    return;

    % ---- nested input validation -----------------------------------------
    function validateInternalEnergyInputs()
        assert(nD == 2, ...
            'evaluateInternalEnergyEndpoints:badShape', ...
            'Endpoint arrays must be NL×2 (nD=%d).', nD);
        assert(all(isfinite(T(:))), ...
            'evaluateInternalEnergyEndpoints:nonFinite', ...
            'estate.T contains non-finite values.');
        sz = [NL, nD];
        nn  = NL + 1;
        szN = [nn, 1];
        flds2 = {'Theta_LL','Theta_II','Theta_g','WW','c_unsat', ...
                 'DThetaLiquidDh','DThetaLiquidDT','DThetaIceDh','DThetaIceDT', ...
                 'DThetaGasDh','DThetaGasDT','dW_dh','dW_dTheta'};
        for i = 1:numel(flds2)
            fn = flds2{i};
            assert(isfield(estate, fn) || isfield(ederiv, fn), ...
                'evaluateInternalEnergyEndpoints:missingField', ...
                'Required field "%s" not found in estate or ederiv.', fn);
            val = getField(estate, ederiv, fn);
            assert(isequal(size(val), sz), ...
                'evaluateInternalEnergyEndpoints:badShape', ...
                'Field "%s" must be NL×2; got %s.', fn, mat2str(size(val)));
            assert(all(isfinite(val(:))), ...
                'evaluateInternalEnergyEndpoints:nonFinite', ...
                'Field "%s" contains non-finite values.', fn);
        end
        fldsN = {'RHOV','RHODA','L_vap','DRHOVh','DRHOVT','Xah','XaT','Xaa'};
        for i = 1:numel(fldsN)
            fn = fldsN{i};
            val = getField(estate, ederiv, fn);
            assert(isequal(size(val), szN), ...
                'evaluateInternalEnergyEndpoints:badShape', ...
                'Field "%s" must be NN×1; got %s.', fn, mat2str(size(val)));
            assert(all(isfinite(val(:))), ...
                'evaluateInternalEnergyEndpoints:nonFinite', ...
                'Field "%s" contains non-finite values.', fn);
        end
        assert(isfinite(L_f) && isscalar(L_f) && L_f >= 0, ...
            'evaluateInternalEnergyEndpoints:badLf', ...
            'ctl.L_f must be a finite nonnegative scalar; got %g.', L_f);
    end

    function v = getField(s1, s2, name)
        if isfield(s1, name), v = s1.(name);
        elseif isfield(s2, name), v = s2.(name);
        else, v = [];
        end
    end
end
