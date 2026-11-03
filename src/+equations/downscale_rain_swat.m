function R_hourly_24 = downscale_rain_swat(R_day, T_dur, t_peak_ratio, start_hour)
%DOWNSCALE_RAIN_SWAT  Downscale daily rainfall to 24 hourly amounts.
%
%   R_hourly_24 = equations.downscale_rain_swat(R_day, T_dur, t_peak_ratio,
%   start_hour) uses the SWAT double-exponential storm pattern to distribute a
%   daily rainfall amount over a finite storm duration.
%
%   Inputs:
%     R_day         [cm] daily rainfall total
%     T_dur         [h] storm duration, optional; random 2..12 h when omitted
%     t_peak_ratio  [-] peak-time fraction in [0,1], optional; random triangular
%                   SWAT-style draw when omitted
%     start_hour    [1..24] start hour index, optional; random when omitted
%
%   Returns:
%     R_hourly_24   [24 x 1, cm] hourly rainfall totals. The output sum is
%                   normalized to R_day.
%
%   Reference:
%     SWAT theoretical documentation, Eqs. 1-67 through 1-82.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    % Treat trace rain below 0.01 mm as no effective rainfall.
    if R_day <= 0.001
        R_hourly_24 = zeros(24, 1);
        return;
    end

    if nargin < 2 || isempty(T_dur)
        T_dur = 2 + rand() * 10;
    end
    T_dur = min(24, T_dur);

    if nargin < 3 || isempty(t_peak_ratio)
        % SWAT triangular distribution for storm peak timing.
        tp_L = 0.05;
        tp_U = 0.95;
        tp_M = 0.25;

        rnd1 = rand();
        threshold = (tp_M - tp_L) / (tp_U - tp_L);

        if rnd1 <= threshold
            t_peak_ratio = tp_L + sqrt(rnd1 * (tp_U - tp_L) * (tp_M - tp_L));
        else
            t_peak_ratio = tp_U - sqrt((1 - rnd1) * (tp_U - tp_L) * (tp_U - tp_M));
        end
    end

    % Decay constants such that endpoint intensity is one percent of peak.
    d1 = t_peak_ratio / 4.605;
    d2 = (1.0 - t_peak_ratio) / 4.605;

    delta1 = d1 * T_dur;
    delta2 = d2 * T_dur;
    T_peak_time = t_peak_ratio * T_dur;

    % Integrated double-exponential area = 0.99 * i_mx * (delta1 + delta2).
    i_mx = R_day / (0.99 * (delta1 + delta2));     % [cm/h]

    time_steps = 0:1:ceil(T_dur);
    cum_R = zeros(length(time_steps), 1);

    for k = 1:length(time_steps)
        t = time_steps(k);
        if t <= T_peak_time
            cum_R(k) = i_mx * delta1 * (exp((t - T_peak_time)/delta1) - exp(-T_peak_time/delta1));
        else
            cum_peak = i_mx * delta1 * (1 - exp(-T_peak_time/delta1));
            cum_R(k) = cum_peak + i_mx * delta2 * (1 - exp((T_peak_time - t)/delta2));
        end
    end

    storm_hourly = max(0, diff(cum_R));

    % Normalize exactly to the input daily rainfall amount.
    sum_storm = sum(storm_hourly);
    if sum_storm > 0
        storm_hourly = storm_hourly * (R_day / sum_storm);
    end

    R_hourly_24 = zeros(24, 1);
    storm_len = length(storm_hourly);

    if nargin < 4 || isempty(start_hour)
        max_start = 24 - storm_len + 1;
        if max_start < 1
            max_start = 1;
        end
        start_hour = randi([1, max_start]);
    end

    end_hour = start_hour + storm_len - 1;
    if end_hour > 24
        valid_len = 24 - start_hour + 1;
        R_hourly_24(start_hour:24) = storm_hourly(1:valid_len);
        lost_rain = sum(storm_hourly(valid_len+1:end));
        R_hourly_24(24) = R_hourly_24(24) + lost_rain;
    else
        R_hourly_24(start_hour:end_hour) = storm_hourly;
    end
end