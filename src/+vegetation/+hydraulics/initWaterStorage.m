function W_state = initWaterStorage(plantHydro, LAI_now) %#ok<INUSD>
%INITWATERSTORAGE  Initialize V2 prognostic plant water storage [W_root; W_stem; W_leaf].
%
%   Strategy: 90% of W_capacity at t=0 (well-watered start).
%   Spin-up of 1-2 days lets W reach quasi-equilibrium with atmospheric demand.
%
%   Inputs:
%     plantHydro  struct  Must contain .W_capacity [3×1, m]
%     LAI_now     [m²/m²] (unused for now; reserved for future LAI scaling)

    if ~isfield(plantHydro, 'W_capacity') || numel(plantHydro.W_capacity) ~= 3
        error('initWaterStorage:NoWCapacity', ...
            'plantHydro.W_capacity (3×1) required. Call mergeV2Params first.');
    end

    init_saturation = 0.9;
    W_state = plantHydro.W_capacity(:) * init_saturation;
end
