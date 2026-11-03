function [gas, audit] = materializeOxygenCoordinate(gas, phase)
%MATERIALIZEOXYGENCOORDINATE Restore the O2 owner without advancing physics.
% GAS is the accepted scalar soil_gas struct. PHASE contains same-shaped
% bottom-up node arrays (NN-by-cell): temperature_C [deg C], gas_pressure_Pa
% and vapor_pressure_Pa [Pa]. The accepted carrier phase must be rebuilt
% before calling this adapter; no initialization or transport is performed.
%
% Existing y_O2 [-] is preserved bitwise and derived C_O2 [mol/cm3 gas] is
% rebuilt with species.dryMoleFractionConcentration. If y_O2 is absent/empty,
% the legacy gas concentration supplies y_O2 = C_O2 / ((Pg-Pv)/(R*T)*1e-6).
% That branch retains the original C_O2 exactly; the new coordinate recovers
% it to floating-point roundoff. All other GAS fields remain unchanged.
% AUDIT marks coordinate migration only, never equivalence of trajectories
% obtained with different oxygen equations. Invalid coordinates or phase
% shapes reject; concentrations and fractions are never clamped.
%
% Source: species.dryMoleFractionConcentration and the coordinate conversion
% introduced in init.initializeSoilAux on 2026-10-23.
% Editor: Wenjie Liu; Date: 2026-10-23.
    assert(isstruct(gas) && isscalar(gas) && isstruct(phase) && isscalar(phase), ...
        'materializeOxygenCoordinate:container', 'Gas and phase must be scalar structs.');
    fields={'temperature_C','gas_pressure_Pa','vapor_pressure_Pa'};
    assert(all(isfield(phase,fields)), ...
        'materializeOxygenCoordinate:phase', 'Accepted T, Pg and Pv are required.');
    shape=size(phase.temperature_C);
    for i=1:numel(fields)
        value=phase.(fields{i});
        assert(isnumeric(value) && isreal(value) && ~isempty(value) ...
            && ismatrix(value) && isequal(size(value),shape), ...
            'materializeOxygenCoordinate:phaseShape', ...
            'Accepted phase fields must be same-shaped nonempty node arrays.');
    end
    dryDensity=soil.gas.species.dryMoleFractionConcentration( ...
        phase.gas_pressure_Pa,phase.temperature_C,ones(shape),phase.vapor_pressure_Pa);
    existingOwner=isfield(gas,'y_O2') && ~isempty(gas.y_O2);
    audit=struct('schema_version',1,'coordinate','dry_gas_mole_fraction', ...
        'coordinate_migration',~existingOwner,'owner_preserved_bitwise',existingOwner, ...
        'trajectory_equivalence_claimed',false, ...
        'legacy_roundtrip_error_mol_cm3',[]);
    if existingOwner
        validateCoordinate(gas.y_O2,shape,'y_O2');
        % Do not reshape, normalize or overwrite the accepted owner.
        gas.C_O2=soil.gas.species.dryMoleFractionConcentration( ...
            phase.gas_pressure_Pa,phase.temperature_C,gas.y_O2,phase.vapor_pressure_Pa);
        audit.policy='owner_preserved_derived_concentration_rebuilt';
    else
        assert(isfield(gas,'C_O2'), ...
            'materializeOxygenCoordinate:missingOxygen', 'Legacy oxygen concentration is required.');
        validateCoordinate(gas.C_O2,shape,'C_O2');
        candidate=gas.C_O2./dryDensity;
        reconstructed=soil.gas.species.dryMoleFractionConcentration( ...
            phase.gas_pressure_Pa,phase.temperature_C,candidate,phase.vapor_pressure_Pa);
        gas.y_O2=candidate;
        audit.policy='legacy_C_O2_to_y_O2_coordinate_migration';
        audit.legacy_roundtrip_error_mol_cm3=max(abs(reconstructed(:)-gas.C_O2(:)));
    end
end

function validateCoordinate(value,shape,name)
% Validate shape and sign; the named species conversion validates y <= 1.
    assert(isnumeric(value) && isreal(value) && isequal(size(value),shape) ...
        && all(isfinite(value(:))) && all(value(:)>=0), ...
        'materializeOxygenCoordinate:coordinate', ...
        'Accepted %s must be a finite nonnegative node array matching phase.',name);
end
