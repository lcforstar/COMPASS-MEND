function [RHS, EnergyMatrices] = calculateBoundaryConditions( ...
        boundary, EnergyMatrices, waterFlux, meteo, state, ...
        precipitation, evaporation, dt_s, RHS, latentHeat, kt)
%CALCULATEBOUNDARYCONDITIONS  Apply the retained heat-equation endpoint rows.
%
%   The active typed boundary adapter maps every production request to one of
%   the following mathematical conditions:
%
%     bottom temperature: T_1 = T_b                          (Dirichlet)
%     bottom energy flux:  R_1 <- R_1 + q_U,bottom^in        (Neumann)
%     top temperature:     T_N = T_s                          (Dirichlet)
%     top energy flux:     R_N <- R_N + q_U,top^in            (Neumann)
%     top response:        q(T_N)=q_ref+dq/dT(T_N-T_ref)      (Robin)
%
%   Public energy flux is positive into the soil. The retained top kernel
%   stores BCT with the opposite sign, hence the code applies `R_N-BCT`.
%   Dirichlet elimination uses the actual strict-lower/superdiagonal entries:
%
%     R_2     <- R_2     - A(2,1) T_b,
%     R_(N-1) <- R_(N-1) - A(N-1,N) T_s.
%
%   Bottom code 4 is the named free-drainage internal-energy condition:
%
%     q_U,bottom^in = q_w,bottom^in rho_l c_l (T_b - 0 degC).
%
%   The water flux is frozen only within the current coupled Picard map and is
%   refreshed with h/T/P_g on the next map. This is the 1-D COMPASS form of
%   ClimaLand EnergyWaterFreeDrainage (local source boundary_conditions.jl,
%   lines 585-615). The top `otherwise` branch alone remains a compatibility
%   path; named atmospheric requests use the matched fixed heat flux.
%
%   Inputs and units:
%     boundary.NBCT/NBCTB  retained top/bottom heat codes
%     boundary.BCT/BCTB    degC for temperature, J cm^-2 s^-1 for flux
%     waterFlux.QMB        cm s^-1, legacy bottom Robin path only
%     precipitation       cm s^-1, legacy surface reconstruction only
%     evaporation         cm s^-1, legacy surface reconstruction only
%     RHS                 J cm^-2 s^-1 equation row, or degC after Dirichlet
%
%   Editor: Jizixing Zhu; Date: 2026-04-25.
%   Change: 2026-07-16 Jizixing Zhu -- rewrite the contract in equation order,
%           document signs/units, and separate typed production paths from
%           inactive legacy water-energy reconstruction branches.
%   Change: 2026-07-17 Jizixing Zhu -- promote bottom liquid sensible-energy export
%           to the named free_drainage_internal_energy Robin condition.
%   Change: 2026-07-19 Jizixing Zhu -- assemble the affine atmospheric
%           surface-energy response in the top diagonal and RHS.
%   Change: 2026-07-20 Jizixing Zhu -- align the natural heat-policy name
%           with the atmospheric water boundary.

    NN = numel(RHS);
    phys = equations.physicalConstants();
    c_a = phys.c_a;
    c_L = phys.c_L;
    RHOL = phys.RHOL;

    C5 = EnergyMatrices.C5;
    C5_a = EnergyMatrices.C5_a;

    % Bottom endpoint: node 1 is the bottom because the grid is bottom-up.
    if boundary.NBCTB == 1
        % Seasonal Tbh, frozen before the current soil trial.
        Tbottom = state.Tbh(kt);
        [RHS, C5, C5_a] = imposeBottomTemperature( ...
            RHS, C5, C5_a, Tbottom);
    elseif boundary.NBCTB == 2
        % Configured fixed bottom temperature.
        Tbottom = boundary.BCTB;
        [RHS, C5, C5_a] = imposeBottomTemperature( ...
            RHS, C5, C5_a, Tbottom);
    elseif boundary.NBCTB == 3
        % q_U,bottom^in is already positive into the column.
        RHS(1) = RHS(1) + boundary.BCTB;
    elseif boundary.NBCTB == 4
        % QMB is positive into the column; free drainage therefore has QMB<0.
        % Moving rho_l*c_l*QMB*T to the left gives the Robin diagonal below.
        assert(isfield(waterFlux, 'QMB') && isscalar(waterFlux.QMB) ...
                && isfinite(waterFlux.QMB), ...
            'ThermalBoundary:MissingBottomWaterFlux', ...
            ['free_drainage_internal_energy requires the finite bottom ', ...
             'water flux from the same coupled Picard map.']);
        C5(1,1) = C5(1,1) - c_L * RHOL * waterFlux.QMB;
    else
        error('ThermalBoundary:UnsupportedBottomCode', ...
            'Unsupported bottom heat boundary code %g.', boundary.NBCTB);
    end

    % Top endpoint: node NN is the soil surface.
    if boundary.NBCT == 1
        Ttop = boundary.BCT;
        if ~isreal(Ttop) || ~isfinite(Ttop)
            Ttop = meteo.Ta(kt);
        end
        [RHS, C5, C5_a] = imposeTopTemperature( ...
            RHS, C5, C5_a, Ttop, NN);
    elseif boundary.NBCT == 2
        % BCT uses the retained outward-positive top sign, so -BCT is the
        % public inward-positive total-energy flux.
        RHS(NN) = RHS(NN) - boundary.BCT;
    elseif boundary.NBCT == 3
        qRef = boundary.BCT;
        Tref = boundary.BCT_ref_temperature;
        dq_dT = boundary.BCT_dflux_dtemperature;
        assert(all(isfinite([qRef, Tref, dq_dT])) && dq_dT <= 1e-14, ...
            'ThermalBoundary:BadSurfaceEnergyResponse', ...
            ['The affine top response requires finite q_ref, T_ref, and ', ...
             'a nonpositive inward-flux temperature derivative.']);
        % A*T = R + q_ref + dq_dT*(T-T_ref).  Moving the trial-T term to
        % the left increases the diagonal because dq_dT is nonpositive.
        C5(NN,1) = C5(NN,1) - dq_dT;
        RHS(NN) = RHS(NN) + qRef - dq_dT * Tref;
    else
        % Legacy fallback only. Named atmospheric requests never
        % enter here; they use the frozen total-energy flux branch above.
        assert(isfield(meteo, 'r_a_SOIL') && isfield(meteo, 'Rn_SOIL'), ...
            'ThermalBoundary:MissingFallbackForcing', ...
            'Legacy heat fallback requires r_a_SOIL and Rn_SOIL.');
        aerodynamicResistance = meteo.r_a_SOIL(kt);
        sensibleHeat = 0.1200 * c_a * ...
            (state.T(NN) - meteo.Ta(kt)) / aerodynamicResistance;
        RHS(NN) = RHS(NN) ...
            + meteo.Rn_SOIL(kt) * 1e-4 ...
            - RHOL * latentHeat(NN) * evaporation(kt) ...
            - sensibleHeat ...
            + RHOL * c_L * (meteo.Ta(kt) * precipitation(kt) ...
                + boundary.DSTOR0 * state.T(NN) / dt_s);
    end

    EnergyMatrices.C5 = C5;
    EnergyMatrices.C5_a = C5_a;
end

function [RHS, C5, C5_a] = imposeBottomTemperature(RHS, C5, C5_a, Tbottom)
% Impose T(1)=Tbottom and eliminate the first-node contribution from row 2.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    RHS(1) = Tbottom;
    C5(1,1) = 1;
    RHS(2) = RHS(2) - C5_a(1) * Tbottom;
    C5(1,2) = 0;
    C5_a(1) = 0;
end

function [RHS, C5, C5_a] = imposeTopTemperature(RHS, C5, C5_a, Ttop, NN)
% Impose T(N)=Ttop and eliminate the surface-node term from row N-1.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    RHS(NN) = Ttop;
    C5(NN,1) = 1;
    RHS(NN-1) = RHS(NN-1) - C5(NN-1,2) * Ttop;
    C5(NN-1,2) = 0;
    C5_a(NN-1) = 0;
end
