function SoilVariables = calculateBottomTempBoundary(SoilVariables, TransportCoefficient, cfg, KT, NoTime_KT)
%CALCULATEBOTTOMTEMPBOUNDARY  Deep-soil harmonic temperature boundary condition.
%
%   Evaluates the semi-infinite solid harmonic solution at the configured total
%   depth and writes hourly (Tbh) and daily (Tb) bottom-temperature arrays.
%   The surface boundary is owned by the canopy energy balance; this function
%   supplies only the lower thermal Dirichlet condition.
%
%   Returns:
%     SoilVariables.Tbh(KT)  [degC] hourly bottom temperature for step KT
%     SoilVariables.Tb       [degC, nDay x 1] daily bottom-temperature state
%
%   Inputs:
%     SoilVariables     [inout] carries Tb (daily cache) and receives Tbh/Tb
%     TransportCoefficient [in] Lambda_eff [internal thermal conductivity]
%                               and c_unsat [J cm^-3 K^-1] volumetric heat capacity
%     cfg.Grid.Tot_Depth  [cm] total soil column depth for the bottom boundary
%     cfg.Env.Tav         [degC] mean annual surface temperature
%     cfg.Env.At          [degC] annual surface temperature amplitude
%     cfg.Time            [struct] start_time, end_time, Delt_t [s] step length
%     KT                  [-] outer-step index (1-based)
%     NoTime_KT           [-] hourly elapsed step index (NoTime counter)
%
%   Method:
%     T(z,t) = Tav + At * exp(-z / De) * cos(2*pi*(JD-JD0)/365 - z/De)
%     De     = sqrt(P * lambda_eff / (pi * C_vol))
%     where P = 365 * 86400 [s] is the annual period, De [cm] is the damping
%     depth, JD [-] is the Julian day of the current step, and JD0 [-] anchors
%     1 January to the 200th Julian day. The harmonic formula is evaluated once
%     per day; within a day the hourly boundary value is flat-spread (deep-soil
%     diurnal variation is negligible).
%
%   Reference:
%     Carslaw, H. S. and Jaeger, J. C. (1959), Conduction of Heat in Solids,
%       2nd ed., Oxford University Press, Ch. 2, Sect. 2.6, p. 81,
%       Eq. (4): temperature in a semi-infinite solid with periodic surface
%       temperature.
%     The same analytical form is adopted by CLM5 (SoilTemperatureMod.F90),
%     Noah-MP, JULES, ORCHIDEE, and CABLE for the deep-soil lower thermal
%     boundary condition.
%   Editor: Jizixing Zhu; Date: 2026-07-13.

    Tot_Depth  = cfg.Grid.Tot_Depth;     % [cm] column bottom depth
    Tav        = cfg.Env.Tav;            % [degC] mean annual surface temperature
    At         = cfg.Env.At;             % [degC] annual amplitude
    Lambda_eff = TransportCoefficient.Lambda_eff;  % internal thermal conductivity units
    c_unsat    = TransportCoefficient.c_unsat;     % [J cm^-3 K^-1] volumetric heat capacity
    Tb         = SoilVariables.Tb;                 % [degC, nDay x 1] daily state cache

    % Column-mean thermal properties across the nD end-point dimension.
    Lambda_effBAR = mean(Lambda_eff, 2);           % [NL x 1]
    c_unsatBAR    = mean(c_unsat, 2);              % [NL x 1]
    Lambda_effAvg = mean(Lambda_effBAR);           % [-] column average
    c_unsatAvg    = mean(c_unsatBAR);              % [-] column average

    % Julian-day series and annual damping depth.
    % jd_series and JD0 depend only on cfg.Time (invariant across steps);
    % they are cached in SoilVariables so the series is built only once.
    % De depends on TransportCoefficient (via Lambda_effAvg / c_unsatAvg)
    % and is recomputed every call.
    if ~isfield(SoilVariables, 'Tb_jd_series') || isempty(SoilVariables.Tb_jd_series)
        SoilVariables.Tb_jd_series = equations.generateJDSeries( ...
            cfg.Time.start_time, cfg.Time.end_time);
        SoilVariables.Tb_JD0 = equations.calculateJulianDay200th( ...
            year(cfg.Time.start_time));
    end
    jd_series = SoilVariables.Tb_jd_series;        % [-] Julian day of each simulation day
    JD0       = SoilVariables.Tb_JD0;               % [-] Julian day of day 200 (anchor)
    De = (365 * 86400 * Lambda_effAvg / (pi() * c_unsatAvg))^0.5;  % [cm] damping depth

    % Convert the hourly NoTime_KT counter to a daily index.
    if NoTime_KT <= 0
        idx = 1;
    else
        idx = floor((NoTime_KT - 1) * cfg.Time.Delt_t / 86400) + 1;
    end
    idx = min(max(idx, 1), numel(jd_series));

    % Harmonic bottom temperature at the current day.
    Tb(idx) = Tav ...
        + At * exp(-Tot_Depth / De) ...
        * cos(2*pi*(jd_series(idx) - JD0) / 365 - Tot_Depth / De);

    % Flat-spread to hourly: the harmonic formula is evaluated daily;
    % deep-soil diurnal temperature change is negligible.
    SoilVariables.Tbh(KT) = Tb(idx);   % [degC] hourly bottom temperature at step KT
    SoilVariables.Tb      = Tb;        % [degC, nDay x 1] daily bottom-temperature state
end
