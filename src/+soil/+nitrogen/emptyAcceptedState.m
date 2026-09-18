function state = emptyAcceptedState()
%EMPTYACCEPTEDSTATE  Allocate the schema for accepted soil-nitrogen memory.
%
%   The actual array dimensions are materialized from the accepted soil
%   endpoint by materializeAcceptedState.  This split is deliberate: static
%   kinetic parameters remain in built.parameters.soil_n, while only history
%   required by the discrete N balance belongs to Y.soil_nitrogen.
%
%   Returns:
%     state [out] struct Accepted N-history schema: prior liquid/ice water
%                 [m3 m^-3], solver variables, and two-level flux ledger.
%
%   Physical meaning:
%     This container contains no mineral-N pool itself.  NH4, NO3, urea,
%   granule, ice N, and O2 remain owned by their existing soil state fields;
%   keeping history separate avoids a duplicate concentration state.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    state = struct( ...
        'schema_version', 1, ...
        'theta_liquid_prev', [], ...
        'theta_ice_prev', [], ...
        'variables', struct(), ...
        'flux', struct(), ...
        'initialized_from_legacy_checkpoint', false);
end
