function Soil = applySurfaceProperties(Soil, state, parameters, spectral)
%APPLYSURFACEPROPERTIES  Blend snow cover into the SCOPE ground surface.
%
%   Snow fraction replaces the effective lower-boundary optical reflectance
%   and thermal emissivity used by the single canopy RT solve. Bare-soil and
%   snow tile properties are retained explicitly for the subsequent tile
%   energy partition. The prognostic snow temperature and the bare-soil skin
%   temperature/history remain separate; neither Ts nor GAM is area-blended.
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.
%   Change: 2026-07-17 Jizixing Zhu -- leave vapor partitioning to the
%     energy-balance tile adapter and always expose the snow diagnostics.
%   Change: 2026-07-18 Jizixing Zhu -- retain explicit bare/snow optical and
%     thermal properties and stop mixing prognostic snow with bare Ts/GAM.

    d = surface.snow.diagnoseState(state, parameters);
    Soil.snow = d;
    Soil.Ts_bare = Soil.Ts;
    Soil.GAM_bare = Soil.GAM;
    Soil.refl_bare = Soil.refl;
    Soil.rs_thermal_bare = Soil.rs_thermal;
    Soil.refl_snow = parameters.albedo .* ones(size(Soil.refl));
    Soil.rs_thermal_snow = 1 - parameters.emissivity;
    if isfield(spectral, 'IwlT') && ~isempty(spectral.IwlT)
        Soil.refl_snow(spectral.IwlT) = Soil.rs_thermal_snow;
    end
    f = d.cover_fraction;
    if f <= 0, return; end

    Soil.refl = (1 - f) .* Soil.refl_bare + f .* Soil.refl_snow;
    Soil.rs_thermal = (1 - f) .* Soil.rs_thermal_bare ...
        + f .* Soil.rs_thermal_snow;
end
