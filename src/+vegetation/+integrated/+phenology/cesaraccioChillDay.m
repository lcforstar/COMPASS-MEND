function [Cd, Ca] = cesaraccioChillDay(Tmin, Tmax, Tmean, baset)
%CESARACCIOCHILLDAY  One-day chill-day (Cd) and anti-chill-day (Ca) increments.
%
%   Cesaraccio et al. (2004) triangular chill/anti-chill-day model as used by
%   CLM5-FruitTree for perennial fruit-tree bud-burst. Chilling effectiveness
%   depends on the diurnal temperature cycle (Tmin/Tmax), not just the daily mean:
%   Cd (<=0) accumulates "chilling" toward the chill requirement crequ; once the
%   requirement is met Ca (>=0) accumulates the post-chilling warmth that drives
%   bud-burst. Both are six-branch piecewise functions of Tmin/Tmax/Tmean vs the
%   base temperature.
%
%   Returns:
%     Cd  [°C·d]  chill-day increment (<= 0; accumulated into chill_day)
%     Ca  [°C·d]  anti-chill-day increment (>= 0; accumulated into anti_chill_day)
%
%   Inputs:
%     Tmin, Tmax  [°C]  daily minimum / maximum air temperature
%     Tmean       [°C]  daily mean air temperature
%     baset       [°C]  base temperature (apple 4 °C)
%
%   Method:
%     Verbatim port of CLM5-FruitTree CNPhenologyMod.F90:1824-1862 (O.Dombrowski;
%     Cesaraccio et al. 2004 Int J Biometeorol 48:67; Fan et al. 2015 CLM-Palm).
%     Branches use temperatures in degrees C. Guards Tmax==Tmin (zero diurnal range).
%
%   Author: Zhu Jizixing 2026-06-04 (stage-13b apple Cesaraccio chilling)

    Cd = 0; Ca = 0;
    dT = Tmax - Tmin;
    if dT <= 0
        % Zero diurnal range: fall back to mean-temperature chilling/forcing.
        if Tmean <= baset && Tmean >= 0
            Cd = -(baset - Tmean);          % below base, above freezing -> chilling
        elseif Tmean > baset
            Ca = Tmean - baset;             % above base -> forcing
        end
        return;
    end

    % --- Chill-day Cd (CNPhenologyMod.F90:1826-1837) ---
    if (0 <= baset) && (baset <= Tmin) && (Tmin <= Tmax)
        Cd = 0;
    elseif (0 <= Tmin) && (Tmin <= baset) && (baset < Tmax)
        Cd = -((Tmean - Tmin) - (Tmax - baset)^2 / (2*dT));
    elseif (0 <= Tmin) && (Tmin <= Tmax) && (Tmax <= baset)
        Cd = -(Tmean - Tmin);
    elseif (Tmin < 0) && (0 <= Tmax) && (Tmax <= baset)
        Cd = -(Tmax^2 / (2*dT));
    elseif (Tmin < 0) && (0 < baset) && (baset < Tmax)
        Cd = -(Tmax^2 / (2*dT)) - ((Tmax - baset)^2 / (2*dT));
    else
        Cd = 0;
    end

    % --- Anti-chill-day Ca (CNPhenologyMod.F90:1847-1858) ---
    if (0 <= baset) && (baset <= Tmin) && (Tmin <= Tmax)
        Ca = Tmean - baset;
    elseif (0 <= Tmin) && (Tmin <= baset) && (baset < Tmax)
        Ca = (Tmax - baset)^2 / (2*dT);
    elseif (0 <= Tmin) && (Tmin <= Tmax) && (Tmax <= baset)
        Ca = 0;
    elseif (Tmin < 0) && (0 <= Tmax) && (Tmax <= baset)
        Ca = 0;
    elseif (Tmin < 0) && (0 < baset) && (baset < Tmax)
        Ca = (Tmax - baset)^2 / (2*dT);
    else
        Ca = 0;
    end
end
