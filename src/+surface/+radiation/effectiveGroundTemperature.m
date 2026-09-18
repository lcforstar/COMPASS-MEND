function temperature_C = effectiveGroundTemperature( ...
        bare_temperature_C, snow_temperature_C, snow_fraction, ...
        bare_emissivity, snow_emissivity)
%EFFECTIVEGROUNDTEMPERATURE  Preserve tile thermal emission in one RT solve.
%
%   Returns the radiative temperature whose emission with the area-weighted
%   emissivity equals the area-weighted emission of bare soil and snow:
%
%     eps_eff T_eff^4 = (1-f) eps_b T_b^4 + f eps_s T_s^4.
%
%   Temperatures are in degC on input and output. The function does not mix
%   prognostic states; it only supplies the equivalent lower-boundary
%   emission required by the single-canopy RTMt solve.
%
% Editor: Jizixing Zhu; Date: 2026-07-18.

    f = min(max(double(snow_fraction), 0), 1);
    epsBare = min(max(double(bare_emissivity), 0), 1);
    epsSnow = min(max(double(snow_emissivity), 0), 1);
    epsEffective = (1 - f) .* epsBare + f .* epsSnow;
    kelvinBare = bare_temperature_C + 273.15;
    kelvinSnow = snow_temperature_C + 273.15;
    emission = (1 - f) .* epsBare .* kelvinBare.^4 ...
        + f .* epsSnow .* kelvinSnow.^4;
    temperature_C = (emission ./ max(epsEffective, eps)).^0.25 - 273.15;
end
