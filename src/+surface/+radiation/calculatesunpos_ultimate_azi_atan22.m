function [solarz, altitude, azi_south, azi_north, R, hour_angle, sunrise_local, sunset_local, day_length, delta] = calculatesunpos_ultimate_azi_atan22(...
    inyear, inmon, inday, gmtime, xlat, xlon, method_choice)
%CALCULATESUNPOS_ULTIMATE_AZI_ATAN2 High-precision solar position and sunrise/sunset calculator
%   Based on subsolar point and atan2 method for solar position calculation
%   Converted from Taiping Zhang's Fortran code
%
%   Input parameters:
%     inyear: Year (4-digit year), e.g., 1998, 2020 [unit: year]
%     inmon: Month (1-12) [unit: month]
%     inday: Day (1-28/29/30/31) [unit: day]
%     gmtime: GMT time (Greenwich Mean Time), decimal hours, e.g., 15.2167 [unit: hour]
%     xlat: Latitude (geographic latitude), positive for northern hemisphere, negative for southern hemisphere [unit: degree]
%     xlon: Longitude (geographic longitude), positive for east longitude, negative for west longitude [unit: degree]
%     method_choice: Hour angle calculation method selection (0=based on longitude difference, 1=based on true solar time) [dimensionless]
%
%   Output parameters:
%     solarz: Solar zenith angle (angle between sun direction and zenith direction) [unit: degree]
%     altitude: Solar altitude angle (angle between sun direction and horizon) [unit: degree]
%     azi_south: Solar azimuth angle (south as 0 deg reference, clockwise increase) [unit: degree]
%     azi_north: Solar azimuth angle (north as 0 deg reference, clockwise increase) [unit: degree]
%     R: Earth-Sun distance (instantaneous distance between Earth and Sun) [unit: astronomical unit AU]
%     hour_angle: Solar hour angle (angular distance of sun relative to local meridian) [unit: degree]
%     sunrise_time: Sunrise time (moment when sun's upper edge touches horizon) [unit: hour, GMT time]
%     sunset_time: Sunset time (moment when sun's upper edge touches horizon) [unit: hour, GMT time]
%     day_length: Day length (time interval from sunrise to sunset) [unit: hour]
%
%   Reference:
%     A solar azimuth formula that renders circumstantial treatment
%     unnecessary without compromising mathematical rigor:
%     Mathematical setup, application and extension of a formula based on
%     the subsolar point and atan2 function
%     
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

% ======================== INPUT VALIDATION ========================
if nargin < 7
    method_choice = 1; % Default method: use true solar time based method
end

% Validate inputs (preserving original validation logic)
% validateattributes(inyear, {'numeric'}, {'scalar', 'integer', '>=', 1900});
validateattributes(inmon, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 12});
validateattributes(inday, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 31});
validateattributes(gmtime, {'numeric'}, {'scalar', '>=', 0, '<=', 24});
validateattributes(xlat, {'numeric'}, {'scalar', '>=', -90, '<=', 90});
validateattributes(xlon, {'numeric'}, {'scalar', '>=', -180, '<=', 180});
validateattributes(method_choice, {'numeric'}, {'scalar', '>=', 0, '<=', 1});

% ======================== CONSTANTS ========================
rpd = pi / 180;                    % Degrees to radians [rad/deg]
UTC_offset = round(xlon / 15);     % Local timezone offset [hours], auto from longitude
central_meridian = UTC_offset * 15; % Timezone central meridian [deg]

% ======================== DATE PROCESSING ========================
% Original calendar calculation logic
[~, julday] = calculateCalendarData(inyear);
dyear = inyear - 2000;              % Years since 2000 [years]
dayofyr = julday(inmon) + inday;    % Day of year [day]
xleap = calculateLeapYears(dyear);  % Leap year count [years]

% ======================== ASTRONOMICAL CALCULATIONS ========================
% Original astronomical calculations
[n, L, g] = calculateAstronomicalParameters(dyear, dayofyr, gmtime, xleap);
[g_rad, lambda_rad, ~, epsilon_rad] = calculateRadianParameters(g, L, n);
[~, alpha, delta] = calculateSolarCoordinates(lambda_rad, epsilon_rad, rpd);
[R, EoT] = calcRandEoT(g_rad, L, alpha);

% ======================== SOLAR GEOMETRY ========================
sunlat = delta;  % Solar latitude (declination) [deg]
sunlon = -15.0 * (gmtime - 12.0 + EoT * 4 / 60);  % Solar longitude [deg]

% ======================== HOUR ANGLE CALCULATION ========================
hour_angle = calculateHourAngle(method_choice, xlon, sunlon, gmtime, EoT);

% ======================== SOLAR POSITION CALCULATION ========================
[solarz, azi_south, azi_north] = calculateSolarPosition(xlat, xlon, sunlat, sunlon, R, rpd);
altitude = 90 - solarz;  % Solar altitude angle [deg]

% ======================== SUNRISE/SUNSET CALCULATION ========================
% Original iterative calculation logic
[sunrise_gmt, sunset_gmt, day_length, ~] = calculateSunriseSunset_Iterative(...
    delta, xlat, xlon, EoT, rpd, UTC_offset, central_meridian, inyear, inmon, inday);

% ======================== TIME CONVERSION ========================
% Original time conversion logic
[sunrise_local, sunset_local] = correctTimeConversion(...
    sunrise_gmt, sunset_gmt, UTC_offset, inyear, inmon, inday);

% ======================== RESULTS DISPLAY ========================
% Original display logic with formatting
%displayResults(method_choice, hour_angle, solarz, altitude, azi_south, azi_north, R, ...
    %sunrise_gmt, sunset_gmt, sunrise_local, sunset_local, day_length, gmtime, iter_info);

end

%% ======================== HELPER FUNCTIONS ========================
% All helper functions below maintain EXACTLY the same computational logic as your original code
% Only formatting and comments have been improved for readability

function [nday, julday] = calculateCalendarData(inyear)
%CALCULATECALENDARDATA Calculate calendar data
%   Calculate monthly days and cumulative days for specified year
%
%   Input parameters:
%     inyear: Year [unit: year]
%
%   Output parameters:
%     nday: Monthly days array [unit: day]
%     julday: Cumulative days array [unit: day]
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    % Monthly days (non-leap year)
    nday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    
    % Leap year judgment: divisible by 4 but not by 100, or divisible by 400
    if (mod(inyear, 100) ~= 0 && mod(inyear, 4) == 0) || ...
       (mod(inyear, 100) == 0 && mod(inyear, 400) == 0)
        nday(2) = 29;  % Leap year February has 29 days
    else
        nday(2) = 28;  % Common year February has 28 days
    end
    
    % Calculate cumulative days
    julday = zeros(1, 13);
    julday(1) = 0;
    for i = 2:13
        julday(i) = julday(i-1) + nday(i-1);
    end
end

function xleap = calculateLeapYears(dyear)
%CALCULATELEAPYEARS Calculate leap year count
%   Calculate number of leap years from 2000 to specified year
%
%   Input parameters:
%     dyear: Year difference relative to 2000 [unit: year]
%
%   Output parameters:
%     xleap: Leap year count [unit: year]
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    if dyear <= 0
        xleap = floor(dyear / 4); % Leap year calculation before 2000
    else
        if mod(dyear, 4) == 0
            xleap = floor(dyear / 4);  % Leap year
        else
            xleap = floor(dyear / 4) + 1; % Common year, "+1" is correction for year 2000
        end
    end
end

function [n, L, g] = calculateAstronomicalParameters(dyear, dayofyr, gmtime, xleap)
%CALCULATEASTRONOMICALPARAMETERS Calculate basic astronomical parameters
%   Calculate solar orbit parameters based on J2000.0 epoch
%
%   Input parameters:
%     dyear: Year difference relative to 2000 [unit: year]
%     dayofyr: Day of year [unit: day]
%     gmtime: GMT time [unit: hour]
%     xleap: Leap year count correction [unit: year]
%
%   Output parameters:
%     n: Julian day number since J2000.0 [unit: day]
%     L: Solar mean longitude [unit: degree]
%     g: Solar mean anomaly [unit: degree]
     % n: Julian day number since J2000.0 epoch (January 1, 2000 12:00 UT)
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    n = -1.5 + dyear * 365.0 + xleap * 1.0 + dayofyr + gmtime / 24;
    
    % L: Solar mean longitude (mean longitude position of sun on ecliptic)
    L = mod(280.460 + 0.9856474 * n, 360.0);
    
    % g: Solar mean anomaly (mean angular position of sun on elliptical orbit)
    g = mod(357.528 + 0.9856003 * n, 360.0);
end

function [g_rad, lambda_rad, epsilon, epsilon_rad] = calculateRadianParameters(g, L, n)
%CALCULATERADIANPARAMETERS Calculate radian-based astronomical parameters
%   Convert angle parameters to radians and calculate related parameters
%
%   Input parameters:
%     g: Solar mean anomaly [unit: degree]
%     L: Solar mean longitude [unit: degree]
%     n: Julian day number [unit: day]
%
%   Output parameters:
%     g_rad: Solar mean anomaly [unit: radian]
%     lambda_rad: Solar true longitude [unit: radian]
%     epsilon: Obliquity of the ecliptic [unit: degree]
%     epsilon_rad: Obliquity of the ecliptic [unit: radian]
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

   % Convert to radians
    rpd = pi / 180;
    g_rad = g * rpd;
    % lambda_rad: Solar true longitude (true longitude considering orbital eccentricity correction)
    lambda_rad = mod(L + 1.915 * sin(g_rad) + 0.020 * sin(2 * g_rad), 360.0) * rpd;
    % epsilon: Obliquity of the ecliptic (angle between ecliptic plane and equatorial plane)
    epsilon = 23.439 - 0.0000004 * n;
    epsilon_rad = epsilon * rpd;
end

function [alpha_rad, alpha, delta] = calculateSolarCoordinates(lambda_rad, epsilon_rad, rpd)
%CALCULATESOLARCOORDINATES Calculate solar coordinates
%   Calculate right ascension and declination of sun
%
%   Input parameters:
%     lambda_rad: Solar true longitude [unit: radian]
%     epsilon_rad: Obliquity of the ecliptic [unit: radian]
%     rpd: Degree to radian conversion factor [dimensionless]
%
%   Output parameters:
%     alpha_rad: Solar right ascension [unit: radian]
%     alpha: Solar right ascension [unit: degree]
%     delta: Solar declination [unit: degree]
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    % alpha_rad: Solar right ascension (longitude in equatorial coordinate system)
    % Use atan2 function to calculate right ascension, which can correctly handle angle calculation in all quadrants
    alpha_rad = atan2(cos(epsilon_rad) * sin(lambda_rad), cos(lambda_rad));
    % alpha: Solar right ascension (convert to degrees)
    alpha = mod(alpha_rad / rpd, 360.0);
    % delta: Solar declination (latitude in equatorial coordinate system)
    delta = asin(sin(epsilon_rad) * sin(lambda_rad)) / rpd;
end

function [R, EoT] = calcRandEoT(g_rad, L, alpha)
%CALCULATEEARTHSUNDISTANCE Calculate Earth-Sun distance and equation of time
%   Calculate distance between Earth and Sun and time difference between true and mean solar time
%
%   Input parameters:
%     g_rad: Solar mean anomaly [unit: radian]
%     L: Solar mean longitude [unit: degree]
%     alpha: Solar right ascension [unit: degree]
%
%   Output parameters:
%     R: Earth-Sun distance [unit: astronomical unit AU]
%     EoT: Equation of time [unit: degree]
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    % R: Earth-Sun distance (calculated based on elliptical orbit model)
    R = 1.00014 - 0.01671 * cos(g_rad) - 0.00014 * cos(2 * g_rad);
    % EoT: Equation of time (difference between true solar time and mean solar time)
    EoT = mod((L - alpha) + 180.0, 360.0) - 180.0;
end

function hour_angle = calculateHourAngle(method_choice, xlon, sunlon, gmtime, EoT)
%CALCULATEHOURANGLE Calculate solar hour angle
%   Calculate solar hour angle using two different methods
%
%   Input parameters:
%     method_choice: Calculation method selection (0=longitude difference, 1=true solar time)
%     xlon: Observer longitude [unit: degree]
%     sunlon: Solar longitude [unit: degree]
%     gmtime: GMT time [unit: hour]
%     EoT: Equation of time [unit: degree]
%
%   Output parameters:
%     hour_angle: Solar hour angle [unit: degree]
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    if method_choice == 0
        % Method 0: Hour angle calculation based on longitude difference
        hour_angle = xlon - sunlon;
    else
        % Method 1: Hour angle calculation based on true solar time
        % Convert GMT time to Beijing time (UTC+8)
        beijing_time_hours = gmtime + 8;
        
        % Calculate local mean solar time (considering longitude time difference)
        local_mean_time = beijing_time_hours + (xlon - 120.0) / 15.0;
        
        % Convert equation of time to hours (15 degrees/hour)
        EoT_hours = EoT / 15.0;
        
        % Calculate true solar time
        true_solar_time = local_mean_time + EoT_hours;
        
        % Calculate hour angle (hour angle = (true solar time - 12) * 15)
        hour_angle = (true_solar_time - 12.0) * 15.0;
    end
    
    % Adjust hour angle to reasonable range [-180, 180] degrees
    hour_angle = mod(hour_angle + 180.0, 360.0) - 180.0;
end

function [solarz, azi_south, azi_north] = calculateSolarPosition(xlat, xlon, sunlat, sunlon, R, rpd)
%CALCULATESOLARPOSITION Calculate solar position angles
%   Calculate solar zenith angle and azimuth angle
%
%   Input parameters:
%     xlat: Observer latitude [unit: degree]
%     xlon: Observer longitude [unit: degree]
%     sunlat: Solar latitude [unit: degree]
%     sunlon: Solar longitude [unit: degree]
%     R: Earth-Sun distance [unit: astronomical unit AU]
%     rpd: Degree to radian conversion factor [dimensionless]
%
%   Output parameters:
%     solarz: Solar zenith angle [unit: degree]
%     azi_south: Solar azimuth angle (south reference) [unit: degree]
%     azi_north: Solar azimuth angle (north reference) [unit: degree]
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    % Convert to radians
    PHIo = xlat * rpd;    % Observer latitude [unit: radian]
    PHIs = sunlat * rpd;  % Solar latitude [unit: radian]
    LAMo = xlon * rpd;    % Observer longitude [unit: radian]
    LAMs = sunlon * rpd;  % Solar longitude [unit: radian]
    
    % Calculate direction cosines
    % Sx, Sy, Sz: Solar direction cosines (components of unit vector pointing to sun in east-north-zenith coordinate system) [dimensionless]
    Sx = R * cos(PHIs) * sin(LAMs - LAMo);% East component
    Sy = R * cos(PHIo) * sin(PHIs) - sin(PHIo) * cos(PHIs) * cos(LAMs - LAMo);% North component
    Sz = R * sin(PHIo) * sin(PHIs) + cos(PHIo) * cos(PHIs) * cos(LAMs - LAMo);% Zenith direction component
   
    % Calculate solar zenith angle and azimuth angle
    solarz = acos(Sz) / rpd;% Solar zenith angle (angle between sun direction and zenith direction) [unit: degree]
    % azi_south: Solar azimuth angle (south as 0 deg reference, clockwise increase)
    % Use atan2 function to calculate azimuth angle, which can correctly handle angle calculation in all quadrants
    azi_south = atan2(-Sx, -Sy) / rpd;% [unit: degree]
    % azi_north: Solar azimuth angle (north as 0 deg reference, clockwise increase)
    azi_north = mod(azi_south + 180, 360); % [unit: degree]
end

function [sunrise_gmt, sunset_gmt, day_length, iter_info] = calculateSunriseSunset_Iterative(...
    delta, xlat, xlon, EoT, rpd, UTC_offset, central_meridian, inyear, inmon, inday)
%CALCULATESUNRISESUNSET_ITERATIVE  Refine sunrise and sunset times iteratively.
%
%   Inputs:
%     delta             [deg] solar declination
%     xlat/xlon         [deg] site latitude/longitude
%     EoT               [deg] equation of time
%     rpd               [rad/deg] degree-to-radian conversion factor
%     UTC_offset        [h] local-standard offset used by the legacy routine
%     central_meridian  [deg] local-standard central meridian
%     inyear/inmon/inday  [-] calendar date
%
%   Returns:
%     sunrise_gmt       [h] sunrise time in GMT decimal hours
%     sunset_gmt        [h] sunset time in GMT decimal hours
%     day_length        [h] daylight duration
%     iter_info         [-] iteration counts for sunrise and sunset
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    
    % Polar day/night check (original logic)
    lat_rad = xlat * rpd;
    delta_rad = delta * rpd;
    h0_rad = -0.833 * rpd;
    
    cos_omega = (sin(h0_rad) - sin(lat_rad) * sin(delta_rad)) / (cos(lat_rad) * cos(delta_rad));
    
    if cos_omega > 1
        sunrise_gmt = NaN; sunset_gmt = NaN; day_length = 0;
        iter_info = [0, 0]; return;
    elseif cos_omega < -1
        sunrise_gmt = NaN; sunset_gmt = NaN; day_length = 24;
        iter_info = [0, 0]; return;
    end
    
    % Initial estimate (original logic)
    [sunrise_approx, sunset_approx] = initialSunEvents(delta, xlat, xlon, EoT, rpd, UTC_offset, central_meridian);
    
    % Iteration parameters (original values)
    max_iter = 5;
    tol = 0.0002778; % ~1 second
    
    % Sunrise iteration (original logic)
    sunrise_old = sunrise_approx;
    for iter = 1:max_iter
        [delta_new, EoT_new] = solarParamsAtTime(inyear, inmon, inday, sunrise_old, xlat, xlon);
        sunrise_new = calculateSunriseTime(delta_new, xlat, xlon, EoT_new, rpd, UTC_offset, central_meridian);
        
        if abs(sunrise_new - sunrise_old) < tol
            sunrise_gmt = sunrise_new;
            break;
        end
        sunrise_old = sunrise_new;
        if iter == max_iter
            sunrise_gmt = sunrise_new;
        end
    end
    
    % Sunset iteration (original logic)
    sunset_old = sunset_approx;
    for iter = 1:max_iter
        [delta_new, EoT_new] = solarParamsAtTime(inyear, inmon, inday, sunset_old, xlat, xlon);
        sunset_new = calculateSunsetTime(delta_new, xlat, xlon, EoT_new, rpd, UTC_offset, central_meridian);
        
        if abs(sunset_new - sunset_old) < tol
            sunset_gmt = sunset_new;
            break;
        end
        sunset_old = sunset_new;
        if iter == max_iter
            sunset_gmt = sunset_new;
        end
    end
    
    % Day length calculation (original logic)
    if sunset_gmt < sunrise_gmt
        day_length = (sunset_gmt + 24) - sunrise_gmt;
    else
        day_length = sunset_gmt - sunrise_gmt;
    end
    
    iter_info = [iter, iter];
end

function [sunrise, sunset] = initialSunEvents(delta, lat, lon, EoT, rpd, UTC_offset, central_meridian)
    % Original initial sun events calculation
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    lat_rad = lat * rpd;
    delta_rad = delta * rpd;
    h0_rad = -0.833 * rpd;
    
    cos_omega = (sin(h0_rad) - sin(lat_rad) * sin(delta_rad)) / (cos(lat_rad) * cos(delta_rad));
    
    if cos_omega > 1 || cos_omega < -1
        sunrise = NaN; sunset = NaN; return;
    end
    
    half_day_length_rad = acos(cos_omega);
    half_day_length_h = half_day_length_rad / 0.2617993878;
    
    EoT_hours = EoT / 15.0;
    longitude_correction = (lon - central_meridian) / 15.0;
    
    local_noon_gmt = 12.0 - UTC_offset - longitude_correction - EoT_hours;
    
    sunrise = local_noon_gmt - half_day_length_h;
    sunset = local_noon_gmt + half_day_length_h;
    
    sunrise = mod(sunrise, 24);
    sunset = mod(sunset, 24);
end

function [delta, EoT] = solarParamsAtTime(inyear, inmon, inday, gmtime, ~, ~)
    % Original solar parameters at specific time
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    [~, julday] = calculateCalendarData(inyear);
    dyear = inyear - 2000;
    dayofyr = julday(inmon) + inday;
    xleap = calculateLeapYears(dyear);
    
    [n, L, g] = calculateAstronomicalParameters(dyear, dayofyr, gmtime, xleap);
    [g_rad, lambda_rad, ~, epsilon_rad] = calculateRadianParameters(g, L, n);
    [~, alpha, delta] = calculateSolarCoordinates(lambda_rad, epsilon_rad, pi/180);
    [~, EoT] = calcRandEoT(g_rad, L, alpha);
end

function sunrise_time = calculateSunriseTime(delta, lat, lon, EoT, rpd, UTC_offset, central_meridian)
    % Original sunrise time calculation
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    [sunrise_time, ~] = initialSunEvents(delta, lat, lon, EoT, rpd, UTC_offset, central_meridian);
end

function sunset_time = calculateSunsetTime(delta, lat, lon, EoT, rpd, UTC_offset, central_meridian)
    % Original sunset time calculation
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    [~, sunset_time] = initialSunEvents(delta, lat, lon, EoT, rpd, UTC_offset, central_meridian);
end

function [sunrise_beijing, sunset_beijing] = correctTimeConversion(sunrise_gmt, sunset_gmt, UTC_offset, ~, ~, ~)
    % Original time conversion logic
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    if isnan(sunrise_gmt)
        sunrise_beijing = NaN; sunset_beijing = NaN; return;
    end
    
    sunrise_beijing_raw = sunrise_gmt + UTC_offset;
    sunset_beijing_raw = sunset_gmt + UTC_offset;
    
    if sunrise_beijing_raw >= 24
        sunrise_beijing = sunrise_beijing_raw - 24;
    elseif sunrise_beijing_raw < 0
        sunrise_beijing = sunrise_beijing_raw + 24;
    else
        sunrise_beijing = sunrise_beijing_raw;
    end
    
    if sunset_beijing_raw >= 24
        sunset_beijing = sunset_beijing_raw - 24;
    elseif sunset_beijing_raw < 0
        sunset_beijing = sunset_beijing_raw + 24;
    else
        sunset_beijing = sunset_beijing_raw;
    end
    
    if sunset_beijing < sunrise_beijing
        sunset_beijing = sunset_beijing + 24;
    end
end

function displayResults(~, hour_angle, solarz, altitude, azi_south, azi_north, R, ...
    sunrise_gmt, sunset_gmt, sunrise_beijing, sunset_beijing, day_length, gmtime, iter_info) %#ok<DEFNU>
%DISPLAYRESULTS  Print a human-readable solar-position diagnostic summary.
%
%   This debug helper is not called by the initialization forcing path; it is
%   kept for manual solar-position checks.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    fprintf('\n=== High-precision solar position with iterative sunrise/sunset ===\n');
    fprintf('Input GMT time: %.2f h\n', gmtime);
    fprintf('Computed hour angle: %.4f deg\n', hour_angle);
    fprintf('Iterations: sunrise %d, sunset %d\n', iter_info(1), iter_info(2));

    fprintf('\nSolar position parameters:\n');
    fprintf('   Solar zenith angle: %.4f deg\n', solarz);
    fprintf('   Solar altitude angle: %.4f deg\n', altitude);
    fprintf('   Solar azimuth angle: %.4f deg (south-reference convention)\n', azi_south);
    fprintf('   Solar azimuth angle: %.4f deg (north-reference convention)\n', azi_north);
    fprintf('   Earth-Sun distance: %.6f AU\n', R);

    if isnan(sunrise_gmt)
        if day_length == 0
            fprintf('\nPolar night: no sunrise or sunset.\n');
        else
            fprintf('\nPolar day: the sun remains above the horizon.\n');
        end
        fprintf('Day length: %.1f h\n', day_length);
        return;
    end

    fprintf('\n1. GMT decimal-hour format:\n');
    fprintf('   Sunrise: %.4f h\n', sunrise_gmt);
    fprintf('   Sunset: %.4f h\n', sunset_gmt);
    fprintf('   Day length: %.4f h\n', day_length);

    fprintf('\n2. Local-time decimal-hour format:\n');
    fprintf('   Sunrise: %.4f h\n', sunrise_beijing);
    fprintf('   Sunset: %.4f h\n', sunset_beijing);
    fprintf('   Day length: %.4f h\n', day_length);

    fprintf('\n3. GMT clock format:\n');
    fprintf('   Sunrise: %s\n', formatTime(sunrise_gmt));
    fprintf('   Sunset: %s\n', formatTime(sunset_gmt));

    fprintf('\n4. Local clock format:\n');
    fprintf('   Sunrise: %s\n', formatTime(sunrise_beijing));
    fprintf('   Sunset: %s\n', formatTime(sunset_beijing));

    day_hours = floor(day_length);
    day_minutes = round((day_length - day_hours) * 60);
    fprintf('   Day length: %d h %d min\n', day_hours, day_minutes);

    fprintf('\n=== Iteration precision check ===\n');
    if iter_info(1) > 1 || iter_info(2) > 1
        fprintf('Iterative refinement completed; expected precision is within about 1 s.\n');
    else
        fprintf('Fast convergence: one iteration reached the precision criterion.\n');
    end
end
function time_str = formatTime(hours)
%FORMATTIME  Convert decimal hours to HH:MM:SS text.
%
%   Inputs:
%     hours     [h] decimal hour value
%
%   Returns:
%     time_str  char formatted as HH:MM:SS
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    total_seconds = round(hours * 3600);
    h = floor(total_seconds / 3600);
    m = floor(mod(total_seconds, 3600) / 60);
    s = mod(total_seconds, 60);
    time_str = sprintf('%02d:%02d:%02d', h, m, s);
end
