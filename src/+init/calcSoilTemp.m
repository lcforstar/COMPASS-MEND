function soilTemperature = calcSoilTemp(varargin)
    %CALCSOILTEMP  Initialize the soil temperature profile (two modes).
    %
    % Mode 1: simple linear gradient (original method)
    %   soilTemperature = calcSoilTemp(T01, T02, delta, domainZ)
    %
    % Mode 2: interpolate from measured data
    %   soilTemperature = calcSoilTemp(depths_measured, T_measured, depths_target)
    %
    % Inputs (mode 1):
    %   T01      - surface temperature [degC]
    %   T02      - bottom temperature [degC]
    %   delta    - total depth [cm]
    %   domainZ  - target depth [cm]
    %
    % Inputs (mode 2):
    %   depths_measured - measurement depths [cm]
    %   T_measured      - measured temperatures [degC]
    %   depths_target   - target depths [cm]
    %
    % Returns:
    %   soilTemperature - interpolated temperature [degC]
    %
    % Mode 1 is a simple linear-gradient interpolation; mode 2 uses linear
    % interpolation + extrapolation and supports multiple measurement points.
    %
    % Reference:
    %   Zeng & Su (2013): STEMMUS (Simultaneous Transfer of Energy, Mass and
    %   Momentum in Unsaturated Soil), Section 5.4.3.
    %
    % Editor: Jizixing Zhu; Date: 2026-06-27.
    % Change: 2026-06-27 Jizixing Zhu -- comments/messages translated to English.

    if nargin == 4
        % ===== mode 1: simple linear gradient =====
        T01 = varargin{1};
        T02 = varargin{2};
        delta = varargin{3};
        domainZ = varargin{4};

        soilTemperature = T02 + domainZ * (T01 - T02) / delta;

    elseif nargin == 3
        % ===== mode 2: interpolate from measured data =====
        depths_measured = varargin{1};
        T_measured = varargin{2};
        depths_target = varargin{3};

        % force column vectors
        depths_measured = depths_measured(:);
        T_measured = T_measured(:);
        depths_target = depths_target(:);

        % linear interpolation + extrapolation
        soilTemperature = init.interpLinearWithExtrapolation(depths_measured, T_measured, depths_target);

    else
        error('calcSoilTemp:badArgCount', ...
            'Wrong number of inputs: expected 4 (mode 1) or 3 (mode 2).');
    end

end
