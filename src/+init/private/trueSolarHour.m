function solarHour = trueSolarHour(hourAngleDeg)
%TRUESOLARHOUR  Apparent solar hour [0..24) from solar hour angle [deg].
%   hour angle = 0 at solar noon, so solar time = 12 + hourAngle/15.
%   hourAngleDeg comes from surface.radiation.calculatesunpos (it already folds in the
%   equation of time), so no manual EoT/longitude handling is needed here.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    solarHour = mod(hourAngleDeg / 15 + 12, 24);
end
