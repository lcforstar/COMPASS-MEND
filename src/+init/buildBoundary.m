function bc = buildBoundary(rawConfig, plan)
%BUILDBOUNDARY  Normalize boundary config into canonical and legacy views.
%
%   bc = init.buildBoundary(rawConfig) builds bc.water, the canonical typed
%   water-boundary policy, and typed heat/dry-air policies, while retaining
%   flat legacy BoundarySettings aliases used by transitional kernels. This is
%   pure config translation: it does not evaluate forcing or update state.
%   Fertilizer and irrigation are management sources, so their schedules belong
%   to built.sources rather than built.boundary.
%
%   Inputs:
%     rawConfig  struct  raw.config (uses .boundary)
%
%   Returns bc with:
%     water:   one site mode per top/bottom endpoint normalized into
%              canonical bc.water plus transitional NBCh_top/NBCh_bottom and
%              h_top/h_bottom aliases. Public fixed flux is [cm/s], positive
%              into the soil domain; sign conversion occurs in the adapter.
%     heat:    named bc.heat policy plus generated NBCT/T compatibility aliases
%     air:     named bc.dry_air policy plus generated NBCP/P compatibility aliases
%     solute:  NBC_{NH4,NO3,UREA}_top/bottom (type 1-3), BC_*_top/bottom [mol/m3]
%     air-rt:  BtmPg [Pa*10] internal bottom air pressure, UseHydrostaticAirBottom [logical]
%
%   Notes:
%     - urea has no boundary transport in site.toml (it enters via fertilizer);
%       the urea BC fields are set to the legacy BoundarySettings defaults.
%     - inactive configured values remain neutral placeholders; runtime forcing
%       is frozen by materializeSoilBoundaryTrial.
%   Reference: legacy BoundarySettings (field names, defaults).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-06-28 Jizixing Zhu — add bc.BtmPg + bc.UseHydrostaticAirBottom (air-pressure runtime BC).
%   Change: 2026-07-01 Jizixing Zhu — build bc.Fertilizer schedule from raw.config.fertilization (replaces legacy setBoundaryCondition fert block).
%   Change: 2026-07-02 Jizixing Zhu — bc.Fertilizer is plain schedule data (no closure); once-only applied-state moved to runtime.management (checkpoint/resume safe).
%   Change: 2026-07-02 Jizixing Zhu — moved the fertilizer schedule OUT of boundary to built.sources.fertilization_schedule (management source, not a BC); boundary no longer has Fertilizer.

%   Change: 2026-07-11 Jizixing Zhu -- canonical typed water BC with strict legacy adapter.
%   Change: 2026-07-11 Jizixing Zhu -- use one flat typed water schema inside the shared [boundary] table.
%   Change: 2026-07-11 Jizixing Zhu -- accept one constant/variable head/flux mode per endpoint; series interpolation is derived automatically.
%   Change: 2026-07-11 Jizixing Zhu -- allow laboratory water-only assembly without inactive heat, air, or solute keys.
%   Change: 2026-07-15 Jizixing Zhu -- add typed heat and dry-air policies;
%           retain flat fields only as compatibility aliases.
%   Change: 2026-07-16 Jizixing Zhu -- replace public heat/dry-air integer
%           codes with named kinds and generate legacy codes only here.
%   Change: 2026-07-17 Jizixing Zhu -- normalize the atmospheric surface-water
%           storage capacity used by the complementarity boundary policy.
%   Change: 2026-07-17 Jizixing Zhu -- require the ClimaLand-style bottom energy
%           outflow to be paired with the same free-drainage water policy.

    b  = rawConfig.boundary;
    bc = struct();
    laboratoryOnly = nargin >= 2 && isfield(plan, 'is_laboratory') ...
        && logical(plan.is_laboratory);

    % Water has one canonical spec. Flat fields remain temporary solver aliases.
    [bc.water, legacyWater] = soil.water.boundary.normalizeSpec(b);
    bc.NBCh_top    = legacyWater.NBCh_top;
    bc.NBCh_bottom = legacyWater.NBCh_bottom;
    bc.h_top       = legacyWater.h_top;
    bc.h_bottom    = legacyWater.h_bottom;
    bc.surface_water = struct( ...
        'capacity_cm', readNonnegativeScalar(b, ...
            'surface_water_capacity_cm', 0));

    if laboratoryOnly
        return;
    end

    % Heat and dry air use named public policies. Integer codes are confined
    % to the compatibility view consumed by the retained matrix adapters.
    [bc.heat, legacyHeat] = soil.thermal.boundary.normalizeSpec(b);
    if strcmp(bc.heat.bottom.kind, 'free_drainage_internal_energy')
        assert(strcmp(bc.water.bottom.kind, 'free_drainage'), ...
            'buildBoundary:incompatibleFreeDrainageEnergy', ...
            ['heat_bottom_kind=free_drainage_internal_energy requires ', ...
             'water_bottom_kind=free_drainage so heat and water use the ', ...
             'same bottom mass flux.']);
    end
    bc.NBCT_top = legacyHeat.NBCT_top;
    bc.NBCT_bottom = legacyHeat.NBCT_bottom;
    bc.T_top = legacyHeat.T_top;
    bc.T_bottom = legacyHeat.T_bottom;
    [bc.dry_air, legacyAir] = soil.gas.boundary.normalizeSpec(b, true);
    bc.NBCP_top = legacyAir.NBCP_top;
    bc.NBCP_bottom = legacyAir.NBCP_bottom;
    bc.P_top = legacyAir.P_top;
    bc.P_bottom = legacyAir.P_bottom;

    % Solute BC codes remain transitional legacy inputs.
    map = { ...
        'nh4_top','NBC_NH4_top'; 'nh4_bottom','NBC_NH4_bottom'; ...
        'nh4_top_value','BC_NH4_top'; 'nh4_bottom_value','BC_NH4_bottom'; ...
        'no3_top','NBC_NO3_top'; 'no3_bottom','NBC_NO3_bottom'; ...
        'no3_top_value','BC_NO3_top'; 'no3_bottom_value','BC_NO3_bottom'};
    for i = 1:size(map, 1)
        bc.(map{i, 2}) = b.(map{i, 1});
    end

    % urea BC: absent in site.toml (urea enters via fertilizer) -> legacy defaults
    bc.NBC_UREA_top = 2; bc.NBC_UREA_bottom = 3;
    bc.BC_UREA_top  = 0; bc.BC_UREA_bottom  = 0;

    % runtime air-pressure boundary helpers (legacy setBoundaryCondition.m:106,127):
    % BtmPg is the internal bottom air pressure [Pa*10]; the bottom Dirichlet air BC
    % is built from surface pressure + a hydrostatic increase when enabled.
    bc.BtmPg = 89821.33497 * 10;          % [Pa*10] STEMMUS internal bottom air pressure
    bc.UseHydrostaticAirBottom = true;
    % Fertilizer is a management SOURCE (not a boundary condition) -> built by
    % init.buildSources as built.sources.fertilization_schedule. Do NOT add it here.
end

function value = readNonnegativeScalar(source, name, fallback)
% Read one optional nonnegative scalar without inventing a site value.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isfield(source, name) && ~isempty(source.(name))
        value = double(source.(name));
    end
    assert(isnumeric(value) && isscalar(value) && isfinite(value) ...
            && value >= 0, ...
        'buildBoundary:badSurfaceWaterCapacity', ...
        'boundary.%s must be one finite nonnegative value [cm].', name);
end
