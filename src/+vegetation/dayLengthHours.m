function hours = dayLengthHours(latitude_deg, doy)
%DAYLENGTHHOURS Astronomical daylight [h] at latitude [deg] and day of year.
% Scalar site latitude and day-of-year inputs; polar night/day are bounded
% by the sunset-angle geometry. The established COMPASS phenology solar
% declination approximation is shared with LUNA; this introduces no new fit.
% Source: integrated.phenology.advanceDay prior local dayLength.
% Editor: Wenjie Liu; Date: 2026-09-30.
    assert(isscalar(latitude_deg)&&isfinite(latitude_deg)&&abs(latitude_deg)<=90 ...
        &&isscalar(doy)&&isfinite(doy)&&doy>=1&&doy<=366, ...
        'vegetation:dayLengthHours:invalidInput','Invalid latitude or day of year.');
    decl = deg2rad(23.45)*sin(2*pi/365*(doy-81));
    cosOmega = min(1,max(-1,-tan(deg2rad(latitude_deg))*tan(decl)));
    hours = 24/pi*acos(cosOmega);
end
