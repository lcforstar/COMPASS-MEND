function built = buildWaterOnlyState(raw, plan)
%BUILDWATERONLYSTATE  Assemble a standalone laboratory Richards column.
%
%   built = init.buildWaterOnlyState(raw, plan) constructs only the common
%   vertical domain, soil-water controls, typed water boundary, prescribed
%   boundary driver, soil parameters, and accepted soil initial state. It does
%   not build meteorological forcing, vegetation, roots, plant C/N, canopy
%   radiation, SCOPE libraries, management sources, or their runtime services.
%   The returned water state uses the same parameter and solver containers as
%   the field pathway, so both contexts call the same Richards kernel.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-12 Jizixing Zhu -- pass the site-selected soil.csv vertical mapping into stage B.

    assert(isfield(plan, 'is_laboratory') && plan.is_laboratory, ...
        'buildWaterOnlyState:badPlan', 'A laboratory initialization plan is required.');

    built = struct();
    gridCfg = init.buildGrid(raw.config);
    built.domain = init.buildDomain(gridCfg);
    built.domain.grid = gridCfg;
    site = raw.config.site;
    elevation = NaN;
    if isfield(site, 'elevation'), elevation = site.elevation; end
    built.domain.site = struct('latitude', site.latitude, ...
        'longitude', site.longitude, 'elevation', elevation);

    built.control = struct();
    built.control.simulation = plan;
    built.control.switches = init.buildSwitches(raw.config);
    built.control.numerics = init.buildNumerics(raw.config);
    built.control.time = init.buildTime(raw.config);
    built.boundary = init.buildBoundary(raw.config, plan);

    rawBoundary = [];
    if isfield(raw, 'water_boundary'), rawBoundary = raw.water_boundary; end
    built.forcing = struct('water_boundary', ...
        init.buildWaterBoundaryForcing(rawBoundary));

    built.switches = built.control.switches;
    built.numerics = built.control.numerics;
    built.time = built.control.time;
    built.grid = built.domain.grid;

    if height(raw.soil) ~= gridCfg.NS
        error('buildWaterOnlyState:soilRowMismatch', ...
            'soil.csv row count (%d) must match grid.NS (%d).', height(raw.soil), gridCfg.NS);
    end
    soilParams = init.buildSoilParameters(raw.soil, built.domain, raw.config);
    bd = init.buildBulkDensity(raw.initial, soilParams, built.domain, raw.config);
    soilParams.bulk_density = bd.layer;
    soilParams.porosity = bd.porosity;
    soilParams = init.buildMoistureChar(soilParams);
    soilParams = init.buildSoilThermal(soilParams);
    soilParams.constants = soil.water.constants();

    built.parameters = struct('soil', soilParams);
    [built.parameters.soil, soilModeState] = ...
        init.buildSoilModes(built, raw.config);
    built.p = struct();
    built.Y = init.buildInitialState(raw.initial, built, raw.config, bd, plan);
    [built.Y, soilModeState] = init.applySoilModesToInitialState( ...
        built.Y, soilModeState, built, raw.config);
    built.Y.soil_modes = soilModeState;
    built.Y.canopy_air = struct();
    built.Y.scope_soil = struct();
end
