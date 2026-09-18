function p = lunaParams(jmaxb1_override)
%LUNAPARAMS  CLM5 LUNA v1.0 physical constants (1:1 from LunaMod.F90 line 55-78).
%
%   Inputs:
%     jmaxb1_override  [-]  optional, override Jmaxb1 with PFT-specific value
%                            (e.g. from vegetation.getPFTValue(pftType, 'Jmaxb1')).
%                            If omitted: default 0.0936 (CLM5 global calibrated).
%
%   Returns p struct with all LUNA constants. Single source of truth — never
%   inline these constants elsewhere. Values are FIXED (CLM5 hardcoded, not
%   tunable) except Jmaxb1 which is namelist-adjustable in CLM5 (per-PFT in
%   COMPASS via PFT table, see getPFTParameterSet group [G] line 65 etc.).
%
%   Reference:
%     Ali et al. (2016) GMD 9:587-606 — LUNA v1.0 model description
%     Xu et al. (2012) New Phytol. 193:30 — N allocation algorithm (origin of Fj25, NUEr25)
%     Lawrence et al. (2019) CLM5 Tech Note Ch.9 — formula numbering
%     Rogers (2014) Photosynthesis Research — Fc25 = 6.22 × 47.3 = 294.2
%     Evans (1989) Oecologia 78:9-19 — Cb chl N for light capture
%     Bernacchi et al. (2001) PCE 24:253 — Kc25/Ko25/Cp25
%     Suzuki et al. (2001) — Rubisco half-life ~7 days → Enzyme_turnover_daily=0.1
%
%   Source file: D:\OneDrive\CTSM-release-clm5.0.37\src\biogeophys\LunaMod.F90:55-78
%
%   Author: Zhu Jizixing, 2026-05-20 (M2-B.3 LUNA implementation)

    % Carbon conversion
    p.Cv = 1.2e-5 * 3600;             % [g C / (μmol CO2 · s) × s/h] → g C/(μmol CO2 · h)

    % Farquhar / Bernacchi (2001) constants at 25°C
    p.Kc25 = 40.49;                   % [Pa] Mechanis constant of CO2 for Rubisco
    p.Ko25 = 27840;                   % [Pa] Mechanis constant of O2 for Rubisco
    p.Cp25 = 4.275;                   % [Pa] CO2 compensation point at 25°C

    % N-to-capacity slopes (LUNA core constants, NOT tunable)
    p.Fc25 = 294.2;                   % [μmol CO2/(g N · s)] Vcmax25 per g N_cb (Rogers 2014)
    p.Fj25 = 1257.0;                  % [μmol e-/(g N · s)] Jmax25 per g N_et (Xu 2012)
    p.NUEr25 = 33.69;                 % [μmol C/(g N · s)] Respiration NUE (Xu 2012)
    p.Cb = 1.78;                      % [-] NUE chl for light capture (Evans 1989)

    % Reference conditions (used by NUEref)
    p.O2ref = 209460.0;               % [ppm] O2 in air
    p.CO2ref = 380.0;                 % [ppm] reference CO2
    p.forc_pbot_ref = 101325.0;       % [Pa] reference air pressure

    % Enzyme dynamics
    p.Q10Enz = 2.0;                   % [-] Q10 for enzyme decay rate
    p.Enzyme_turnover_daily = 0.1;    % [1/day] daily turnover at 25°C (Suzuki 2001)

    % Allocation baselines (Jmaxb1 is namelist-tunable in CLM5, per-PFT in COMPASS)
    p.Jmaxb0 = 0.0311;                % [-] baseline PNet
    p.Jmaxb1 = 0.0936;                % [-] CLM5 global calibrated default
                                       %     (per-PFT value comes via jmaxb1_override from
                                       %      PFT table getPFTParameterSet; apple is recalibrated there)
    if nargin >= 1 && ~isempty(jmaxb1_override) && isfinite(jmaxb1_override)
        p.Jmaxb1 = jmaxb1_override;   % override with PFT-specific value
    end
    p.Wc2Wjb0 = 0.8054;               % [-] baseline Wc/Wj ratio
    p.relhExp = 6.0999;               % [-] RH influence on J
    p.NMCp25 = 0.715;                 % [g C/(g N · day)] maintenance C cost coefficient
                                       % (assumes 80% of maintenance respiration for photo enzyme)

    % Structural / boundary
    p.SNC = 0.004;                    % [g N / g C] structural N concentration
    p.mp = 9.0;                       % [-] stomatal conductance slope (may be overridden by PFT)
    p.Trange1 = 5.0;                  % [°C] lower T limit for N optimization
    p.Trange2 = 42.0;                 % [°C] upper T limit for N optimization
    p.PARLowLim = 200.0;              % [μmol/m²/s] minimum PAR for N optimization
    p.minrelh = 0.25;                 % [-] minimum relative humidity for N optimization

    % Fallback ceilings/floors (LunaMod.F90 line 392-419)
    p.Vcmax25_fallback = 50.0;        % [μmol/m²/s] seed value & NaN/extreme fallback
    p.Jmax25_fallback  = 85.0;        % [μmol/m²/s] seed value & NaN/extreme fallback
    p.Vcmax25_max = 1000.0;           % [μmol/m²/s] hard upper bound (triggers fallback)
    p.Jmax25_max  = 2000.0;           % [μmol/m²/s] hard upper bound

end
