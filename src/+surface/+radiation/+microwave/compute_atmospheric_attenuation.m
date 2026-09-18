function atm = compute_atmospheric_attenuation(params, state)
%COMPUTE_ATMOSPHERE  L-band atmospheric opacity and up/downwelling brightness.
%
%   Returns:
%     atm.tau_atm   [-]   Zenith atmospheric opacity
%     atm.gossat    [-]   Path transmission = exp(-tau_atm/cosθ)
%     atm.tb_up     [K]   Upwelling atmospheric brightness temperature at TOA
%     atm.tb_dn     [K]   Downwelling brightness temperature reaching TOC
%     atm.T_aeq     [K]   Effective atmospheric equivalent temperature
%     atm.T_cosmic  [K]   Cosmic microwave background (constant 2.725 K)
%
%   Inputs:
%     params.calib.{iAtm, Inc}
%     params.site.{Z}                Site elevation [km]
%     state.Tair_K                    Surface air temperature [K]
%
%   Mode switch (params.calib.iAtm):
%     0 = none (tau=1, T_up=T_dn=0)         — backward compatible (default)
%     1 = Pellarin L-band (CMEM)            — recommended for 1.4 GHz
%     2 = Ulaby (CMEM)                       — wider 1–100 GHz support
%
%   Method (iAtm=1, Pellarin):
%       tau_atm = exp(-3.926 − 0.2211·Z − 0.00369·Tair)
%       T_aeq   = exp( 4.927 + 0.002195·Tair)
%       gossat  = exp(-tau_atm/cosθ)
%       Tb_dn   = T_aeq·(1 − gossat) + T_cosmic·gossat
%       Tb_up   = T_aeq·(1 − gossat)
%
%   References:
%     Pellarin, Wigneron, Calvet, Waldteufel 2003, IEEE-TGRS 41(9), 2098–2103
%       — empirical L-band sky temperature model (CATDS/SMOS heritage)
%     Source code: ECMWF CMEM v5 atm_sub.F90 (Apache 2.0)
%       data/references/CMEM/CMEM/atm_sub.F90 SUBROUTINE ATMPELLARIN
%
%   Notes:
%     - Pellarin coefficients are calibrated for L-band (1.4 GHz) only.
%     - tau_atm typical L-band clear sky ≈ 0.010–0.015 (so gossat ≈ 0.985–0.99
%       at 40° incidence). T_aeq ≈ 4–6 K.
%     - T_cosmic = 2.725 K (Penzias & Wilson 1965) added to downwelling sky
%       since CMEM omits it (CMEM is for satellites looking down — we are
%       computing TOA brightness so it matters here).

if ~isfield(params.calib, 'iAtm'), params.calib.iAtm = 0; end

atm.T_cosmic = 2.725;

switch params.calib.iAtm
    case 0
        atm.tau_atm = 0;
        atm.gossat  = 1;
        atm.tb_up   = 0;
        atm.tb_dn   = 0;
        atm.T_aeq   = NaN;
        return;

    case 1
        % Pellarin L-band (CMEM ATMPELLARIN, atm_sub.F90:46–59)
        Z_km    = params.site.Z;
        Tair_K  = state.Tair_K;
        atm.tau_atm = exp(-3.926 - 0.2211*Z_km - 0.00369*Tair_K);
        atm.T_aeq   = exp( 4.927 + 0.002195*Tair_K);

    case 2
        % Ulaby fixed table (CMEM ULABY_TAU_ATM)
        fghz = params.calib.Fre;
        ROSWAT = 7.5;                                  % g/m³ surface vapor density
        if fghz <= 11
            atm.tau_atm = 0.04;
        elseif abs((fghz - 15)/15) <= 0.1
            atm.tau_atm = 0.055 + 0.004*ROSWAT;
        elseif abs((fghz - 22.235)/22.235) <= 0.1
            atm.tau_atm = 0.11  + 0.048*ROSWAT;
        elseif abs((fghz - 35)/35) <= 0.1
            atm.tau_atm = 0.17  + 0.013*ROSWAT;
        elseif abs((fghz - 90)/90) <= 0.1
            atm.tau_atm = 0.17  + 0.06*ROSWAT;
        else
            atm.tau_atm = 0.04;                        % fallback to low-freq
        end
        atm.T_aeq = state.Tair_K;                       % Ulaby uses surface air T

    otherwise
        error('compute_atmospheric_attenuation: unsupported iAtm = %d', params.calib.iAtm);
end

atm.gossat = exp(-atm.tau_atm / cosd(params.calib.Inc));
atm.tb_dn  = atm.T_aeq * (1 - atm.gossat) + atm.T_cosmic * atm.gossat;
atm.tb_up  = atm.T_aeq * (1 - atm.gossat);
end
