function TransportCoefficient = calculateGasDispersivity(SoilVariables, TransportCoefficient, cfg)
%CALCULATEGASDISPERSIVITY Compute gas pressure mobility and dispersion.
%   The gas mechanical-dispersion coefficient is
%       D_Vg = alpha_Lg*abs(V_A),
%   where V_A = -(k_g/MU_a)*dP_g/dz is the Darcy gas flux and alpha_Lg is
%   the saturation-dependent longitudinal dispersivity.
%
%   Inputs:
%     SoilVariables.P_gg   - Trial gas pressure [g cm^-1 s^-2], NN x 1.
%     TransportCoefficient.k_g - Intrinsic gas permeability [cm^2], NL x nD.
%     TransportCoefficient.Sa  - Actual gas saturation Theta_g/POR [-], NL x nD.
%     cfg                  - Grid dimensions and layer thickness [cm].
%
%   Updated outputs in TransportCoefficient:
%     Alpha_Lg - longitudinal gas dispersivity [cm], NL x nD.
%     D_Vg     - longitudinal gas dispersion [cm^2 s^-1], NL x 1.
%     V_A      - Darcy gas flux [cm s^-1], NL x 1.
%     Kg       - gas pressure mobility [cm^3 s g^-1], NL x nD.
%
%   The Zeng/STEMMUS correlation is retained because no universal
%   soil-gas dispersivity law exists. Sa is bounded to its physical interval
%   before evaluating the fifth-order correlation; this leaves all physical
%   states unchanged and prevents numerical water/ice overshoot from
%   producing an unphysical dispersivity. Alpha_Lg is exposed for diagnosis
%   and future site-scale calibration instead of remaining a hidden local.
%
%   Reference: Yu et al. (2021), Geosci. Model Dev. 14:7345-7376,
%   Supplement Sect. S2.2.7, Eqs. (S2.29)-(S2.30), DOI:
%   10.5194/gmd-14-7345-2021. Local copy:
%   data/references/soil_gas_transport/Yu_et_al_2021_STEMMUS_UEB_supplement.pdf.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- rename the gas mobility field to Kg without
%   changing the coefficient values or Darcy flux calculation.
%   Change: 2026-07-13 Jizixing Zhu -- bound actual gas saturation, expose
%   Alpha_Lg, and vectorize the documented STEMMUS gas-dispersion closure.

NL         = cfg.Grid.NL;
DeltZ      = cfg.Grid.DeltZ;
KaT_Switch = SoilVariables.KaT_Switch;

P_gg = SoilVariables.P_gg(:);
k_g  = TransportCoefficient.k_g;
Sa   = TransportCoefficient.Sa;

MU_a = equations.getPhysicalConstant('MU_a');  % g cm^-1 s^-1

if cfg.lDEBUG_NAN_TRACE
    fprintf('[calculateGasDispersivity] Input: P_gg NaN=%d, k_g NaN=%d\n', ...
        sum(isnan(P_gg(:))), sum(isnan(k_g(:))));
end

% Yu et al. (2021), Supplement Eq. (S2.30). The reference value is
% alpha_L,sat=0.078 m=7.8 cm. Sa is the actual gas-filled pore fraction
% Theta_g/POR, not the effective mobile-gas saturation used by k_g.
SaBounded = min(1, max(0, Sa));
alphaLSaturatedCm = 0.078 * 100;
Alpha_Lg = alphaLSaturatedCm .* ...
    (13.6 - 16 .* SaBounded + 3.4 .* SaBounded.^5);

Kg = k_g ./ MU_a;
KgBAR = 0.5 .* (Kg(:, 1) + Kg(:, 2));
Alpha_LgBAR = 0.5 .* (Alpha_Lg(:, 1) + Alpha_Lg(:, 2));
DPgDZ = diff(P_gg) ./ DeltZ(:);
V_A = -KgBAR .* DPgDZ;

D_Vg = zeros(NL, 1);
if KaT_Switch == 1
    D_Vg = Alpha_LgBAR .* abs(V_A);
end

if cfg.lDEBUG_NAN_TRACE
    outOfRange = isfinite(Sa) & (Sa < 0 | Sa > 1);
    if any(outOfRange(:))
        fprintf('[calculateGasDispersivity] bounded Sa at %d endpoint(s).\n', ...
            nnz(outOfRange));
    end
    badVelocity = find(isnan(V_A));
    for k = 1:numel(badVelocity)
        ML = badVelocity(k);
        fprintf('  [CRITICAL@ML=%d] V_A=NaN. KgBAR=%.2e, DPgDZ=%.2e\n', ...
            ML, KgBAR(ML), DPgDZ(ML));
    end
    fprintf('[calculateGasDispersivity] Output: V_A NaN=%d, D_Vg NaN=%d\n', ...
        sum(isnan(V_A(:))), sum(isnan(D_Vg(:))));
end

TransportCoefficient.Alpha_Lg = Alpha_Lg;
TransportCoefficient.D_Vg     = D_Vg;
TransportCoefficient.V_A      = V_A;
TransportCoefficient.Kg       = Kg;
end
