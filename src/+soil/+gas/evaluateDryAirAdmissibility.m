function status = evaluateDryAirAdmissibility( ...
        density, storage, referenceDensity, referenceStorage, remainingFraction)
%EVALUATEDRYAIRADMISSIBILITY Check the physical dry-air feasible domain.
%
%   Dry-air density rho_a [g cm^-3] and dry-air storage
%
%     B_a = rho_a (theta_g + H_c theta_l) [g cm^-3]
%
%   must remain nonnegative. When reference values are supplied, the trial
%   must also retain `remainingFraction` of every positive reference value.
%   This is the fraction-to-boundary condition used by the coupled h/U/P_g
%   line search; it scales the complete correction and never clips a state.
%
%   Inputs:
%     density, storage                 current constitutive values [g cm^-3]
%     referenceDensity/referenceStorage feasible iteration values [g cm^-3]
%     remainingFraction                fraction in [0,1); zero checks only
%                                      the closed physical domain
%
%   Returns:
%     status.admissible                true only for finite feasible values
%     status.minimum_density_g_cm3     minimum nodal rho_a [g cm^-3]
%     status.minimum_storage_g_cm3     minimum endpoint B_a [g cm^-3]
%     status.failure_reason            stable diagnostic identifier
%
%   The constitutive values must be produced by the active production
%   Density_DA/calculateDryAirParameters path. This function introduces no
%   alternative vapor-pressure formula or thermodynamic closure.
%
%   Editor: Jizixing Zhu; Date: 2026-08-12.

    if nargin < 3 || isempty(referenceDensity)
        referenceDensity = [];
    end
    if nargin < 4 || isempty(referenceStorage)
        referenceStorage = [];
    end
    if nargin < 5 || isempty(remainingFraction)
        remainingFraction = 0;
    end

    density = density(:);
    storage = storage(:);
    assert(~isempty(density) && ~isempty(storage) ...
            && isscalar(remainingFraction) && isfinite(remainingFraction) ...
            && remainingFraction >= 0 && remainingFraction < 1, ...
        'evaluateDryAirAdmissibility:badInput', ...
        ['Dry-air density/storage must be nonempty and remainingFraction ', ...
         'must be a finite scalar in [0,1).']);

    densityFloor = zeros(size(density));
    storageFloor = zeros(size(storage));
    if ~isempty(referenceDensity) || ~isempty(referenceStorage)
        referenceDensity = referenceDensity(:);
        referenceStorage = referenceStorage(:);
        assert(isequal(size(referenceDensity), size(density)) ...
                && isequal(size(referenceStorage), size(storage)) ...
                && all(isfinite(referenceDensity)) ...
                && all(isfinite(referenceStorage)) ...
                && all(referenceDensity >= 0) ...
                && all(referenceStorage >= 0), ...
            'evaluateDryAirAdmissibility:badReference', ...
            ['Reference dry-air density/storage must be finite, ', ...
             'nonnegative, and shape-compatible.']);
        densityFloor = remainingFraction * referenceDensity;
        storageFloor = remainingFraction * referenceStorage;
    end

    finiteValues = all(isfinite(density)) && all(isfinite(storage));
    densityFeasible = finiteValues && all(density >= densityFloor);
    storageFeasible = finiteValues && all(storage >= storageFloor);
    if ~finiteValues
        failureReason = 'nonfinite_dry_air_state';
    elseif ~densityFeasible
        failureReason = 'negative_or_boundary_crossing_dry_air_density';
    elseif ~storageFeasible
        failureReason = 'negative_or_boundary_crossing_dry_air_storage';
    else
        failureReason = '';
    end

    status = struct( ...
        'admissible', densityFeasible && storageFeasible, ...
        'minimum_density_g_cm3', min(density), ...
        'minimum_storage_g_cm3', min(storage), ...
        'failure_reason', failureReason, ...
        'remaining_fraction', remainingFraction);
end
