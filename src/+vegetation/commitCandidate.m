function built = commitCandidate(built, candidate, kt)
%COMMITCANDIDATE Install the accepted vegetation-domain state candidate.
%
%   Inputs:
%     built     [struct] Legacy state shell after soil/surface acceptance.
%     candidate [struct] Vegetation water, C/N, phenology, and acclimation
%               candidates for the accepted outer interval.
%     kt        [-] Accepted one-based outer-step index.
%
%   Returns:
%     built [struct] Shell with accepted vegetation states and LUNA memory.
%
%   Physical contract:
%     Plant hydraulic storage, canopy interception water, vegetation state,
%     plant C/N, and delayed acclimation memories are vegetation-owned.  LUNA
%     is intentionally advanced only now, preserving the established one-step
%     lag and preventing rejected soil trials from contaminating plant memory.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    required = {'vegetation','canopy_water','leaf_resp_Ta', ...
        'luna_climate_entry'};
    for i = 1:numel(required)
        assert(isstruct(candidate) && isfield(candidate, required{i}) ...
                && ~isempty(candidate.(required{i})), ...
            'vegetation:commitCandidate:missingCandidate', ...
            'Vegetation candidate.%s is required at accepted commit.', ...
            required{i});
    end
    if ~isempty(candidate.plant_cn)
        built = init.scatterNitrogenVariables(built, candidate.plant_cn, kt);
    end
    plantWater = candidate.plant_hydraulic_W;
    if ~isempty(plantWater)
        assert(numel(plantWater) == 3 && all(isfinite(plantWater(:))), ...
            'vegetation:commitCandidate:badHydraulicStorage', ...
            'plant_hydraulic_W must contain finite root, stem, and leaf water.');
        built.Y.plant_hydraulic.W_root = plantWater(1);
        built.Y.plant_hydraulic.W_stem = plantWater(2);
        built.Y.plant_hydraulic.W_leaf = plantWater(3);
    end
    built.Y.canopy_water = candidate.canopy_water;
    built.Y.vegetation = candidate.vegetation;
    built.Y.plant_acclimation.leaf_resp.Ta_history(kt) = ...
        candidate.leaf_resp_Ta;
    lunaBuffer = vegetation.integrated.LUNA.Climate10dBuffer.fromStruct( ...
        built.Y.plant_acclimation.luna_climate_buffer);
    lunaBuffer.push(kt, candidate.luna_climate_entry);
    [lunaStateNext, lunaBufferNext] = ...
        vegetation.integrated.LUNA.runDailyLuna(struct('built', built), kt, lunaBuffer);
    built.Y.plant_acclimation.luna_state = lunaStateNext;
    built.Y.plant_acclimation.luna_climate_buffer = lunaBufferNext;
end
