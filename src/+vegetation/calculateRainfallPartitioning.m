function [CanopyInterception_KT, Throughfall_KT] = calculateRainfallPartitioning(P_KT, Delt_t, LAI, kV, S_L, F0)
%CALCULATERAINFALLPARTITIONING  Gash sparse-canopy rainfall interception rates.
%
%   [CanopyInterception_KT, Throughfall_KT] = vegetation.calculateRainfallPartitioning( ...
%       P_KT, Delt_t, LAI, kV, S_L, F0)
%   splits one outer-step precipitation rate into canopy interception and
%   throughfall using a Gash-type sparse canopy model. Pure: no reservoir state
%   is updated. Downstream vegetation.computeCanopyInterception adds CI to h2ocan,
%   applies drip, and computes fwet.
%
%   Returns:
%     CanopyInterception_KT  [cm/s]  Potential canopy interception rate
%     Throughfall_KT         [cm/s]  Throughfall rate; equals P_KT - CanopyInterception_KT
%
%   Inputs:
%     P_KT    [in]  [cm/s]     Precipitation rate for this step
%     Delt_t  [in]  [s]        Outer step length
%     LAI     [in]  [m2/m2]    Leaf area index
%     kV      [in]  [-]        Canopy extinction coefficient (canopy.kV)
%     S_L     [in]  [mm/LAI]   Canopy storage capacity per unit LAI
%     F0      [in]  [-]        Free-throughfall / wet-canopy ratio parameter
%                              used as f_ER = f_v * F0 in this implementation
%
%   Method:
%     Convert P_KT to step rainfall P_mm [mm]. If LAI <= 0 or P_mm <= 0, return
%     zero interception and full throughfall. Otherwise form canopy cover
%     f_v = 1 - exp(-kV*LAI), storage Sv = S_L*LAI [mm], and f_ER = f_v*F0.
%     When 0 < f_ER < f_v, use the Gash saturation rainfall
%     Pg1 = -Sv/f_ER * log(1 - f_ER/f_v); intercept f_v*P_mm below Pg1 and
%     f_v*Pg1 + f_ER*(P_mm - Pg1) above. Clip interception to [0, P_mm] and
%     convert rates back to cm/s. Does not evaporate or update h2ocan.
%
%   Reference:
%     Gash, J. H. C. (1979). An analytical model of rainfall interception by
%     forests. Quarterly Journal of the Royal Meteorological Society, 105,
%     43-55. https://doi.org/10.1002/qj.49710544304
%     [SOURCE UNVERIFIED against local PDF 2026-07-16] Equation symbols follow
%     the common sparse-canopy Gash implementation used in COMPASS; Pg1 and
%     f_ER wiring match the prior COMPASS/STEMMUS-style code path, not a
%     re-derivation from a checked local copy of Gash (1979).
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-16 Jizixing Zhu — English compass-coding-standard doc string
%           and inline comments; body unchanged.

    % Step rainfall [mm] for comparison with storage Sv [mm].
    P_mm = P_KT * Delt_t * 10;   % [mm] from [cm/s] * [s]

    if LAI <= 0 || P_mm <= 0
        CanopyInterception_KT = 0;
        Throughfall_KT        = P_KT;
        return;
    end

    % Gash geometric / storage parameters
    f_v  = 1 - exp(-kV * LAI);   % [-] canopy cover fraction
    Sv   = S_L * LAI;            % [mm] canopy storage capacity
    f_ER = f_v * F0;             % [-] wet-canopy evaporation-to-rainfall ratio proxy

    if f_ER <= 0 || f_ER >= f_v
        Ei_mm = f_v * P_mm;      % [mm] interception when f_ER is out of range
    else
        Pg1 = -Sv / f_ER * log(1 - f_ER / f_v);   % [mm] rainfall to saturate canopy
        if P_mm < Pg1
            Ei_mm = f_v * P_mm;
        else
            Ei_mm = f_v * Pg1 + f_ER * (P_mm - Pg1);
        end
    end

    Ei_mm = max(0, min(Ei_mm, P_mm));   % [mm] interception cannot exceed rainfall

    % Convert interception back to [cm/s]
    Ei_cms                = Ei_mm / 10 / Delt_t;
    CanopyInterception_KT = Ei_cms;
    Throughfall_KT        = P_KT - Ei_cms;
end
