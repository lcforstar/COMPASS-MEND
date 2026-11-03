function registerSnow(reg, cfg) %#ok<INUSD>
%REGISTERSNOW  ClimaLand-compatible bulk-snow state and flux outputs.
% Output: snow.nc
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-17 Jizixing Zhu -- expose conservative drainage when the last
%   solid snow phase melts out.
% Change: 2026-07-18 Jizixing Zhu -- expose the enthalpy carried by terminal liquid
%   drainage so the snow and receiving-surface energy ledgers remain closed.
% Change: 2026-07-18 Jizixing Zhu -- document terminal drainage of solid
%   tails below the coupled latent-energy resolution.
    g = 'snow';
    scalar(reg, g, 'snow_water_equivalent', ...
        'Snow water equivalent', 'cm');
    scalar(reg, g, 'snow_liquid_water', ...
        'Retained liquid water in snow', 'cm');
    scalar(reg, g, 'snow_ice_water_equivalent', ...
        'Ice water equivalent in snow', 'cm');
    scalar(reg, g, 'snow_internal_energy', ...
        'Bulk snow internal energy', 'J m-2');
    scalar(reg, g, 'snow_temperature', ...
        'Bulk snow temperature', 'degC');
    scalar(reg, g, 'snow_depth', 'Snow depth', 'm');
    scalar(reg, g, 'snow_density', 'Bulk snow density', 'kg m-3');
    scalar(reg, g, 'snow_cover_fraction', 'Snow-covered ground fraction', '1');
    rateAmount(reg, g, 'snowfall', 'Post-canopy snowfall');
    rateAmount(reg, g, 'rainfall', 'Post-canopy liquid rainfall');
    rateAmount(reg, g, 'snow_runoff', 'Liquid drainage from snow');
    rateAmount(reg, g, 'snow_complete_meltout', ...
        ['Terminal liquid drainage when solid snow is exhausted or below ', ...
         'model energy resolution']);
    scalar(reg, g, 'snow_complete_meltout_energy_flux', ...
        'Terminal meltwater enthalpy flux; positive into soil', 'W m-2');
    scalar(reg, g, 'snow_complete_meltout_energy', ...
        'Terminal meltwater enthalpy transferred this step', 'J m-2');
    rateAmount(reg, g, 'snow_vapor_loss', ...
        'Snow evaporation or sublimation; negative denotes deposition');
    scalar(reg, g, 'frozen_surface_water_to_snow', ...
        'Liquid depression storage transferred into snow this step', 'mm');
    scalar(reg, g, 'snow_soil_heat_flux', ...
        'Conductive snow-to-soil heat flux; positive downward', 'W m-2');
    scalar(reg, g, 'snow_surface_energy_flux', ...
        'Unclipped ClimaLand snow energy flux; positive toward atmosphere', 'W m-2');
    scalar(reg, g, 'snow_applied_energy_flux', ...
        'Applied clipped ClimaLand snow energy flux; positive toward atmosphere', 'W m-2');
    scalar(reg, g, 'snow_phase_change_flux', ...
        'Liquid phase-change tendency; ClimaLand outward-positive convention', 'cm s-1');
    scalar(reg, g, 'snow_water_balance_residual', ...
        'Bulk-snow total-water step residual', 'mm');
    scalar(reg, g, 'snow_energy_balance_residual', ...
        'Bulk-snow energy step residual', 'J m-2');
end

function scalar(reg, group, name, longName, units)
% Register one scalar time series.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    reg.register(struct('name',name, 'long_name',longName, 'units',units, ...
        'dims',{{'time'}}, 'group',group, 'transform','scalar'));
end

function rateAmount(reg, group, prefix, longName)
% Register a canonical [cm/s] rate and accepted-step [mm] amount.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    scalar(reg, group, [prefix '_flux'], [longName ' flux'], 'cm s-1');
    scalar(reg, group, [prefix '_amount'], [longName ' amount'], 'mm');
end
