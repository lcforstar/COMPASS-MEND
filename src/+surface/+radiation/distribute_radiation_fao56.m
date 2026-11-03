function R_hourly = distribute_radiation_fao56(R_day, DOY, days_in_year, latitude_deg, is_watt, solar_hour_centers)
%DISTRIBUTE_RADIATION_FAO56  Distribute daily shortwave radiation to 24 hours.
%
%   R_hourly = surface.radiation.distribute_radiation_fao56(R_day, DOY,
%   days_in_year, latitude_deg, is_watt, solar_hour_centers) applies the
%   FAO-56/Liu-Jordan daylight weighting curve to split a daily shortwave total
%   or daily-mean shortwave flux into forcing-clock hourly values.
%
%   Inputs:
%     R_day         [MJ/m2/d] or [W/m2] daily shortwave amount or daily mean flux
%     DOY           [-] day of year, 1..366
%     days_in_year  [-] 365 or 366 for the Gregorian forcing-clock year
%     latitude_deg  [deg] latitude, positive north
%     is_watt       [-] true when R_day is daily-mean W/m2; false for MJ/m2/d
%     solar_hour_centers [24 x 1 h] apparent-solar-time centers corresponding
%                   to forcing-clock intervals 00:00-01:00 through 23:00-24:00
%
%   Returns:
%     R_hourly      [24 x 1] hourly shortwave values in the same convention as
%                   R_day. Nighttime values are zero; the discrete daylight
%                   weights are normalized to conserve the supplied daily total.
%
%   Reference:
%     Liu and Jordan 1960 Solar Energy 4:1-19; FAO-56 Eq. 1.6.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-08-11 Jizixing Zhu -- use the configured Gregorian year
%     length, map apparent-solar centers onto forcing-clock intervals, and
%     normalize the 24 discrete weights for exact daily conservation.

    if days_in_year ~= 365 && days_in_year ~= 366
        error('distribute_radiation_fao56:badYearLength', ...
            'days_in_year must be 365 or 366, got %g.', days_in_year);
    end
    solar_hour_centers = solar_hour_centers(:);
    if numel(solar_hour_centers) ~= 24 || any(~isfinite(solar_hour_centers))
        error('distribute_radiation_fao56:badSolarHourCenters', ...
            'solar_hour_centers must contain 24 finite apparent-solar hours.');
    end

    phi   = deg2rad(latitude_deg);
    delta = 0.409 * sin((2*pi/days_in_year) * DOY - 1.39);   % [rad] solar declination

    % Sunset hour angle with polar-day/night clipping.
    cos_arg = -tan(phi) * tan(delta);
    cos_arg = max(-1, min(1, cos_arg));
    omega_s = acos(cos_arg);                         % [rad]

    omega = (pi/12) .* (solar_hour_centers - 12);    % [rad] noon omega = 0

    % Liu-Jordan hourly fraction: positive in daylight, zero at night.
    rt = zeros(24, 1);
    day_mask = abs(omega) < omega_s;
    if any(day_mask)
        num   = cos(omega(day_mask)) - cos(omega_s);
        denom = sin(omega_s) - omega_s * cos(omega_s);
        rt(day_mask) = (pi/24) .* (num ./ denom);
    end
    weight_sum = sum(rt);
    if weight_sum > 0
        rt = rt ./ weight_sum;
    end

    if is_watt
        % Daily-mean W/m2 -> hourly W/m2 with an exact 24-hour mean.
        R_hourly = R_day .* rt .* 24;
    else
        % Daily MJ/m2/d -> hourly MJ/m2: split the daily total directly.
        R_hourly = R_day .* rt;
    end

    R_hourly = max(0, R_hourly);
end
