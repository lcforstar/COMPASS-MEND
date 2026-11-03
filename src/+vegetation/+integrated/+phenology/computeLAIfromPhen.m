function LAI = computeLAIfromPhen(cfg, pftType, phen, lai_epic, C_leaf_g_m2)
%COMPUTELAIFROMPHEN  Compute canopy LAI from phenology and the selected LAI engine.
%
%   Dispatches to one of two LAI engines based on cfg.Veg.LAI_engine:
%
%     'epic_phen_mix'   : LAI = (1 - phen) * ALAI_MIN + phen * EPIC_curve.
%                         This is the calibrated heat-unit/EPIC path. The
%                         phenology factor blends the empirical EPIC curve
%                         toward ALAI_MIN during dormancy and transition periods.
%
%     'leafc_sla_phen'  : LAI = LAI_full(C_leaf), with no final phen multiplier.
%                         In this path, phenology acts upstream by controlling
%                         C_leaf through allocation, onset transfer, and
%                         senescence. The LAI function only diagnoses canopy
%                         structure from the accepted prognostic leaf carbon
%                         pool.
%
%                         CLM5-style leafC-to-LAI conversion:
%                           if SLAm > 0:
%                             LAI_full = (SLA0 / SLAm) *
%                                        (exp(C_leaf * SLAm) - 1)
%                           else:
%                             LAI_full = SLA0 * C_leaf
%                         where SLA0 maps to CLM5 slatop and SLAm maps to
%                         dsladlai.
%
%   Returns:
%     LAI      [m²/m²]  Final canopy LAI for the current time step.
%
%   Inputs:
%     cfg          [in]  config struct (uses cfg.Veg.LAI_engine, cfg.Veg.ALAI_MIN,
%                        cfg.Veg.LAI_max)
%     pftType      [in]  PFT name string (e.g. 'Apple')
%     phen         [in]  [-] phenology factor in [0, 1]; used by 'epic_phen_mix'
%                        and kept for signature compatibility in 'leafc_sla_phen'
%     lai_epic     [in]  [m²/m²] EPIC frPHU S-curve LAI from cfg.Veg.lai_model_current
%                        (used by 'epic_phen_mix' engine; ignored by 'leafc_sla_phen')
%     C_leaf_g_m2  [in]  [g C/m²] current leaf C pool (used by 'leafc_sla_phen';
%                        ignored by 'epic_phen_mix')
%
%   References:
%     CLM5 CNVegStructUpdateMod.F90:152-156 (Thornton & Zimmerman 2007 Eq 3)
%     Williams JR (1989) Trans ASAE 32:497 (EPIC LAI S-curve)
%
%   Created 2026-05-26 (Phase D LAI engine factorization).

ALAI_MIN = cfg.Veg.ALAI_MIN;
LAI_max  = cfg.Veg.LAI_max;
engine   = lower(vegetation.getVegetationProperty(cfg.Veg, 'LAI_engine', {}, 'epic_phen_mix'));

switch engine
    case 'epic_phen_mix'
        % --- Engine 1: EPIC × phen linear blend ---
        % Preserves Williams 1989 EPIC GA-calibrated S-curve at full canopy;
        % uses phen to handle dormant/onset/offset boundaries cleanly.
        if isempty(lai_epic) || ~isfinite(lai_epic) || lai_epic < 0
            lai_epic = ALAI_MIN;
        end
        LAI = (1 - phen) * ALAI_MIN + phen * lai_epic;

    case 'leafc_sla_phen'
        % --- Engine 2: prognostic leaf C -> CLM5-style canopy LAI ---
        SLA0 = vegetation.getPFTValue(pftType, 'SLA0');     % [m²/g C] CLM5 slatop
        SLAm = vegetation.getPFTValue(pftType, 'SLAm');     % [m²/g C / LAI] CLM5 dsladlai

        % Use the prognostic C_leaf pool directly (no N_leaf x CN_L reverse-derive)
        leafC_g_m2 = max(0, C_leaf_g_m2);              % [g C / m²]

        % CLM5 CNVegStructUpdateMod.F90:152-156 maps leafc, slatop, and dsladlai to tlai.
        if SLAm > 1e-9
            % CLM5 Thornton 2007 Eq 3: integrated SLA(LAI) = SLA0 + SLAm × LAI
            LAI_full = (SLA0 / SLAm) * (exp(leafC_g_m2 * SLAm) - 1);
        else
            % Single-SLA fallback when no canopy SLA gradient is specified.
            LAI_full = SLA0 * leafC_g_m2;
        end
        LAI_full = max(0, min(LAI_full, LAI_max));

        % Phenology already acts upstream on C_leaf; do not multiply by phen again here.
        LAI = LAI_full;

    otherwise
        warning('computeLAIfromPhen:unknownEngine', ...
            'LAI_engine="%s" unknown, falling back to epic_phen_mix', engine);
        if isempty(lai_epic) || ~isfinite(lai_epic) || lai_epic < 0
            lai_epic = ALAI_MIN;
        end
        LAI = (1 - phen) * ALAI_MIN + phen * lai_epic;
end

end
