function [sun_altitude, sun_azimuth] = calculateSunPosition(lat, lon, year, month, day, hour, minute)
%CALCULATESUNPOSITION 计算指定地点、日期和时间的太阳高度角和方位角
%   使用Bourges算法计算太阳位置，基于经纬度、日期和时间输入
%
%   输入参数：
%     lat    : 纬度（度），北纬为正
%     lon    : 经度（度），东经为正
%     year   : 年份
%     month  : 月份
%     day    : 日期
%     hour   : 小时（24小时制）
%     minute : 分钟
%
%   输出参数：
%     sun_altitude : 太阳高度角（度），地平线为0°，天顶为90°
%     sun_azimuth  : 太阳方位角（度），正北为0°，顺时针增加
%
%   参考算法：Bourges太阳位置算法
%   编写规范遵循MATLAB编程最佳实践[7,8](@ref)

    % 输入参数验证
    validateattributes(lat, {'numeric'}, {'scalar', '>=', -90, '<=', 90});
    validateattributes(lon, {'numeric'}, {'scalar', '>=', -180, '<=', 180});
    validateattributes(year, {'numeric'}, {'integer', 'scalar', 'positive'});
    validateattributes(month, {'numeric'}, {'integer', 'scalar', '>=', 1, '<=', 12});
    validateattributes(day, {'numeric'}, {'integer', 'scalar', '>=', 1, '<=', 31});
    validateattributes(hour, {'numeric'}, {'scalar', '>=', 0, '<=', 23});
    validateattributes(minute, {'numeric'}, {'scalar', '>=', 0, '<=', 59});

    %% 计算年积日（一年中的第几天）
    start_of_year = datenum(year, 1, 1);
    current_date = datenum(year, month, day);
    day_of_year = current_date - start_of_year;
    
    %% 计算太阳赤纬角（Bourges算法）
    % 计算n0参数（年份校正）
    n0 = 78.801 + (0.2422 * (year - 1969)) - fix(0.25 * (year - 1969));

    % 计算日角b（弧度）
    b = 2 * pi * (day_of_year - n0 - 1) / 365.2422;

    % 计算太阳赤纬角（度）
    solar_declination = 0.3723 + ...
                        23.2567 * sin(b) + ...
                        0.1149 * sin(2 * b) - ...
                        0.1712 * sin(3 * b) - ...
                        0.758 * cos(b) + ...
                        0.3656 * cos(2 * b) + ...
                        0.02010 * cos(3 * b);

    % 转换为弧度
    solar_declination_rad = deg2rad(solar_declination);
    lat_rad = deg2rad(lat);
    
    %% 计算太阳时角
    % 计算时差Ts（分）
    Ts = 0.0028 - 1.9587 * sin(b) + 9.9059 * sin(2 * b) - ...
         7.0924 * cos(b) - 0.6882 * cos(2 * b);
    
    % 计算平太阳时（小时）
    mean_solar_time = hour + (minute - (120.0 - lon) * 4.0) / 60.0;
    
    % 计算真太阳时（小时）
    true_solar_time = mean_solar_time + Ts / 60;
    
    % 计算时角（弧度）
    hour_angle = (true_solar_time - 12.0) * pi / 12.0;
    
    
    %% 计算太阳高度角
    sin_altitude = sin(solar_declination_rad) * sin(lat_rad) + ...
                   cos(solar_declination_rad) * cos(lat_rad) * cos(hour_angle);
    sun_altitude_rad = asin(sin_altitude);
    sun_altitude = rad2deg(sun_altitude_rad);
    
    %% 计算太阳方位角
    % 计算方位角余弦值
    cos_azimuth = (sin(sun_altitude_rad) * sin(lat_rad) - sin(solar_declination_rad)) / ...
                  (cos(sun_altitude_rad) * cos(lat_rad));
    
    % 确保余弦值在有效范围内[-1, 1]
    cos_azimuth = max(min(cos_azimuth, 1), -1);
    
    % 计算方位角（弧度）
    sun_azimuth_rad = acos(cos_azimuth);
    sun_azimuth = rad2deg(sun_azimuth_rad);
    
    % 方位角校正：以南为0°转换为以北为0°
    sun_azimuth = 180 - sun_azimuth;
    
    % 根据时角调整方位角范围
    if hour_angle > 0
        sun_azimuth = 360 - sun_azimuth;
    end
    
    % 确保方位角在0-360度范围内
    sun_azimuth = mod(sun_azimuth, 360); 
end

% 示例使用代码
% [altitude, azimuth] = vegetation.hydraulics.calculateSunPosition(30, 120, 2024, 6, 1, 12, 0)