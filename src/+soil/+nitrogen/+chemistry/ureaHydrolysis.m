function [Shys, Vpore, WFPS, DeltaUrea] = ureaHydrolysis(K_urea_hydrolysis, NL, Nurea, Vpore, WFPS, Theta_g, Theta_LL, Shys, Delt_t, TT, DeltZ)
%UREAHYDROLYSIS  Urea → NH4 hydrolysis with WFPS modulation + Q10 temperature correction.
%
% First-order reaction analytically integrated over Δt (avoids numerical issues
% when k_eff × Δt is large).
%
% Formulation:
%   k_eff(ML,ND) = 5 × WFPS(ML,ND) × K_urea_T(ML,ND)            [1/s]
%   K_urea_T     = K_urea_hydrolysis × Q10^((T − T_ref)/10)       [1/s]
%   ΔN (concentration drop over Δt):
%     ΔN = N_urea × (1 − exp(−k_eff × Δt))                        [μg/cm³]
%   Shys (rate, for ADE G1 source term):
%     Shys = ΔN / Δt                                              [μg/(cm³·s)]
%
% Inputs:
%   K_urea_hydrolysis [in]    [1/s, NL×1]       base hydrolysis rate (Excel /86400)
%   NL                [in]    [-]               number of layers
%   Nurea             [in]    [μg/cm³, NN×1]   prior-step urea concentration
%   Vpore             [inout] [-, NL×2]        pore volume = θ_g + θ_LL (overwritten)
%   WFPS              [inout] [-, NL×2]        water-filled pore space (overwritten)
%   Theta_g           [in]    [m³/m³, NL×2]    gas-phase content
%   Theta_LL          [in]    [m³/m³, NL×2]    liquid content
%   Shys              [inout] [μg/(cm³·s), NL×2]  hydrolysis rate (overwritten)
%   Delt_t            [in]    [s]              outer step length
%   TT                [in]    [°C, NN×1]       soil temperature at nodes
%
% Returns:
%   Shys      [μg/(cm³·s), NL×2]  hydrolysis rate (for ADE source term G1)
%   Vpore     [-, NL×2]
%   WFPS      [-, NL×2]
%   DeltaUrea [μg/cm³, NN×1]      urea concentration drop this Δt (= ΔN)
%
% Reference:
%   Lloyd & Taylor 1994 — Q10 temperature dependence
%   Liang Hao 2020 thesis — WFPS exponential-saturation modulation (factor 5)
%
% History:
%   2026-05-14: fixed unit bug — Shys was returning concentration-drop ΔN instead
%     of rate; now Shys = ΔN/Δt matching G1 source-term units. DeltaUrea (concentration
%     drop, for caller's Urea pool update) decoupled from Shys.

% DeltZ [cm, NLx1]: bottom-up layer thicknesses (required).
% Editor: Jizixing Zhu; Date: 2026-09-07.
% Change: 2026-09-07 Jizixing Zhu -- project endpoint reaction mass onto
% shared nodes using half-layer storage weights, rather than overwriting.
% sum(node_capacity .* DeltaUrea) = dt * sum(dz .* mean(Shys,2)).
    assert(numel(DeltZ)==NL && all(isfinite(DeltZ(:)) & DeltZ(:)>0) ...
        && isscalar(Delt_t) && isfinite(Delt_t) && Delt_t>0, ...
        'ureaHydrolysis:badGridOrTime','Positive layer thickness and time are required.');
    DeltZ=DeltZ(:);
    Q10 = 2.5;
    T_ref = 20;

    DeltaUrea = zeros(size(Nurea));

    for ML = 1:NL
        for ND = 1:2
            MN = ML + ND - 1;

            Vpore(ML, ND) = Theta_g(ML, ND) + Theta_LL(ML, ND);
            WFPS(ML, ND) = Theta_LL(ML, ND) ./ max(Vpore(ML, ND), eps);

            % Temperature-adjusted urea hydrolysis rate constant (Lloyd & Taylor 1994)
            K_urea_temp = K_urea_hydrolysis(ML) * Q10^((TT(MN) - T_ref) / 10);

            % First-order effective rate with WFPS modulation (Liang Hao 2020)
            k_eff = 5.0 * WFPS(ML, ND) * K_urea_temp;       % [1/s]

            % Analytical concentration drop over Δt
            dN = Nurea(MN) * (1 - exp(-k_eff * Delt_t));    % [μg/cm³]

            % Average rate over Δt (for ADE G1 NH4 source, μg/(cm³·s))
            Shys(ML, ND) = dN / max(Delt_t, eps);

            % Concentration drop for caller's urea-pool update
            DeltaUrea(MN) = DeltaUrea(MN) + dN*DeltZ(ML)/2;
        end
    end
    capacity=([DeltZ;0]+[0;DeltZ])/2;
    DeltaUrea=DeltaUrea./capacity;
end
