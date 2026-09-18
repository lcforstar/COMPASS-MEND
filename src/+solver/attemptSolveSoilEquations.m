function [model, outcome] = attemptSolveSoilEquations(model, kt, dt, solveCtl)
%ATTEMPTSOLVESOILEQUATIONS  Solve one switched scalar soil trial substep.
%
%   [model, outcome] = solver.attemptSolveSoilEquations(model, kt, dt, solveCtl) solves ONE trial
%   substep. Water is always solved by modified Picard; heat and dry air are added when their
%   master switches are enabled. The four supported combinations are water, water+heat,
%   water+dry-air, and water+heat+dry-air. The common sequence is: A3 extrapolate the
%   guess -> A5 dual-porosity sink -> A4 Picard/Anderson setup -> for KIT=1:NIT { apply the accelerated
%   guess -> refresh all enabled equation coefficients once -> assemble their common tridiagonal
%   matrices -> solve water -> optionally solve dry air -> optionally solve heat -> commit the
%   Anderson map -> test the coupled accept criterion }. The water nonlinear
%   acceptance rule is HYDRUS-style: maximum endpoint total-water-equivalent change where both
%   iterates are unsaturated is below TolTh, maximum head change at saturated or saturation-crossing
%   nodes is below TolH, and the actual-flux water-balance residual is below tol_water_mass_mm. A failed mass
%   check continues the same modified-Picard iteration until all gates pass or NIT is exhausted. Ports the
%   water branch of the legacy Picard loop COMPASS.m:1657-1856. It solves ONE trial substep only: it
%   does NOT tile the hour, does NOT snapshot/restore, does NOT commit the outer step. It writes ONLY
%   trial scratch; never Y/pending/the runtime clock. Freeze-thaw phase storage
%   uses the same scalar water/heat solve; nitrogen remains outside this solver.
%
%   Readable call hierarchy for one trial substep:
%
%     attemptSolveSoilEquations
%     +-- initialize immutable accepted-state anchors
%     +-- for each Modified-Picard iteration
%     |   +-- synchronize h/U/T/P_g trial coordinates
%     |   +-- soil.refreshSoilCoefficients
%     |   |   +-- soil.refreshWaterAux
%     |   |   +-- soil.refreshThermalGasTransport
%     |   |   +-- soil.buildEquationCoefficients       [shared flux owner]
%     |   |       +-- water coefficient bank
%     |   |       +-- one shared element-flux family
%     |   |       +-- dry-air coefficient bank from shared flux
%     |   |       +-- energy coefficient bank from shared flux
%     |   +-- soil.assembleEquationMatrices
%     |   +-- ordinary or block Modified-Picard map    [globalized solve]
%     |   +-- soil.coupling.evaluateResidualLedger     [exact ledger owner]
%     |   |   +-- exact water/dry-air/energy residual blocks
%     |   |   +-- one column ledger per conserved quantity
%     |   |   +-- exact top state and exact top water flux
%     |   +-- soil.water.boundary.evaluateSurfaceTrial
%     |       +-- surface storage ledger
%     |       +-- atmospheric complementarity decision
%     +-- return a trial outcome; the outer step alone may commit
%
%   Frozen coefficients describe the Picard operator. The exact ledger owner
%   evaluates the updated candidate produced by that operator. They therefore
%   use two different nonlinear states by design; one candidate state is
%   evaluated exactly once and that record is the sole acceptance authority.
%
%   Inputs:
%     model     [inout] model container. Reads/updates model.built (trial scratch: workspace.step.
%                       trial/.derived/.solver via the called units). Reads model.runtime.substep.
%                       Delt_old [s] (the A3 previous-substep dt; read-only).
%     kt        [in]    [-] outer-step index (KT).
%     dt        [in]    [s] trial substep length (Delt_t).
%     solveCtl  [in]    struct with .NIT [-], .tol_theta [-], .tol_head_saturated_cm [cm], and
%                       .tol_water_mass_mm [mm]. Optional .tol_surface_water_cm
%                       [cm] gates the atmospheric surface-reservoir identity;
%                       its default is tol_water_mass_mm converted from mm to cm.
%                       Optional .water_boundary_flux contains the current
%                       pending EVAP/Throughfall, and .water_boundary_request contains the immutable
%                       boundary request materialized for this exact trial interval.
%                       Optional .tol_air_mass_g_cm2 is the accepted dry-air
%                       column conservation tolerance [g cm^-2].
%                       Optional .tol_air_mass_relative is the accepted
%                       dry-air column residual divided by the accepted
%                       column dry-air mass [-]. The production default is
%                       1e-4, with tol_air_mass_g_cm2 retained as an absolute
%                       floor for small-mass columns.
%                       Optional .tol_energy_balance_J_cm2 is the accepted
%                       internal-energy conservation tolerance [J cm^-2].
%
%   Returns:
%     model    with model.built trial scratch advanced by the water solve (trial.h etc.); Y unchanged.
%     outcome  struct. converged [-] means all HYDRUS nonlinear and mass gates passed; residual [cm]
%              retains max(CHKWATER) as diagnostic telemetry only; iters [-] is the actual iteration
%              count. final_theta_increment [-], final_saturated_head_increment_cm [cm],
%              final_active_set_changed [-], and final_mass_residual_mm [mm] expose every acceptance
%              component. surface_converged [-], surface_switch_count [-],
%              and final_surface_mass_residual_cm [cm] expose the independent
%              atmospheric complementarity gate. termination_reason is
%              'accepted', 'nit_exhausted', 'surface_active_set_cycle',
%              'temperature_inversion_failed', or
%              'dry_air_admissibility_failed'.
%
%   References:
%     HYDRUS benchmark controls: tests/benchmarks/sv62/code/prepare_hydrus_sv62.m
%       (MaxIt, TolTh, and TolH input record).
%     Legacy Anderson map wiring: COMPASS.m:1665-1675 and 1883-1887.
%     Gatti et al. (2024), arXiv:2407.07865, Sect. 2.2 Eq. (2):
%       no-ponding seepage is a pressure/flux Signorini condition. COMPASS
%       retains its finite-depth DSTOR reservoir and uses the same
%       complementarity structure without adopting Nitsche or mixed FEM.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-08 Jizixing Zhu -- WBC-2: return the converged tile waterFlux in outcome for the hourly accumulator.
%   Change: 2026-07-09 Jizixing Zhu -- I4: wire Anderson apply(KIT>1)/commit into the water Picard (legacy COMPASS.m:1665-1675,1883-1887) via the A4 setup; marshal trial.{h,T} through pack/unpack; outcome.anderson stats; lAnderson default true.
%   Change: 2026-07-09 Jizixing Zhu -- add candidate actual-mass residual telemetry and gating.
%   Change: 2026-07-09 Jizixing Zhu -- on accept, refresh derived from the accepted trial.h so the committed Y.Theta_L (via carryThetaLayer) is VG-consistent with Y.h (closes the pre-solve/post-solve theta ruler gap; telemetry == committed WBC).
%   Change: 2026-07-13 Jizixing Zhu -- merge the former water-only engine into the
%           switch-controlled water/heat/dry-air solver as its sole implementation.
%   Change: 2026-07-13 Jizixing Zhu -- require relative dry-air pressure
%           convergence and admit the existing freeze-thaw phase closure.
%   Change: 2026-07-13 Jizixing Zhu -- include p_gas in the scaled Anderson
%           vector whenever the dry-air equation is active.
%   Change: 2026-07-09 Jizixing Zhu -- consume this-step water-boundary forcing solveCtl.water_boundary_flux (EVAP/Throughfall; default built.workspace.flux) in BOTH the solve and the mass scorer, so the solve/telemetry use the current-step pending forcing (not the last-committed, one-step-stale built.workspace.flux).
%   Change: 2026-07-10 Jizixing Zhu -- compute accepted theta_eq_n once per
%           substep and thread it into the conservative water time term.
%   Change: 2026-07-12 Jizixing Zhu -- remove water_time_form; conservative
%           mixed storage is the sole water-equation time formulation.
%   Change: 2026-07-12 Jizixing Zhu -- keep accepted endpoint water-equivalent
%           conversion local after consolidating the mixed storage helpers.
%   Change: 2026-07-10 Jizixing Zhu -- anchor mixed theta_eq_n on current accepted Y.h/Y.T/Y.p_gas, matching legacy h/T; workspace.step.previous is the hOLD/TOLD extrapolation anchor.
%   Change: 2026-07-10 Jizixing Zhu -- resolve solveCtl.water_root_sink (this-step root sink Srt; default committed built.workspace.flux.Srt) ONCE and pass it to every solveWaterTrial call (fixed across all substeps); backward-compatible when absent.
%   Change: 2026-07-10 Jizixing Zhu -- score candidate mass with FE-endpoint liquid-plus-vapor thetaEq storage, matching the mixed-form time term instead of the liquid/node-collapse WBC ruler.
%   Change: 2026-07-11 Jizixing Zhu -- extend the same FE storage ruler with ice
%     liquid-water-equivalent mass for freeze-thaw states.
%   Change: 2026-07-11 Jizixing Zhu -- W1b: pass the immutable per-trial BoundaryRequest to every Picard kernel call and return it as audit-only outcome metadata; no request is committed to model state.
%   Change: 2026-07-11 Jizixing Zhu -- W1b: return the active mathematical boundary kind separately so atmospheric head/flux switching is visible in audit results.
%   Change: 2026-07-11 Jizixing Zhu -- add an opt-in HYDRUS convergence policy:
%     max dtheta in unsaturated endpoints and max dh in saturated nodes.
%   Change: 2026-07-16 Jizixing Zhu -- expose heat/dry-air Applied boundary
%           kinds, recover each accepted coupled boundary reaction once, and
%           require the exact dry-air column mass identity before acceptance.
%   Change: 2026-07-12 Jizixing Zhu -- retain the staged freeze latent heat in every
%     accepted-state, Picard, and mass-scoring refresh.
%   Change: 2026-07-12 Jizixing Zhu -- make HYDRUS TolTh/TolH plus finite
%     actual-mass closure the sole acceptance policy; remove legacy raw-head
%     selection, adaptive head-tolerance tightening, and their telemetry.
%   Change: 2026-07-13 Jizixing Zhu -- extend the new scalar Picard engine to
%     optional dry-air and heat equations. All equations share one coefficient
%     refresh, shared FEM assembly, and shared Thomas kernel; process RHS,
%     boundaries, safeguards, and accepted-state commit remain separate.
%   Change: 2026-07-13 Jizixing Zhu -- anchor the heat time term to accepted
%     endpoint internal energy and use an exact mixed energy-storage residual.
%   Change: 2026-07-13 Jizixing Zhu -- anchor dry-air mass to the accepted
%     endpoint storage and use an exact mixed dry-air residual.
%   Change: 2026-07-14 Jizixing Zhu -- stabilize the strict VG-Mualem gas
%     permeability with bounded vector-Aitken relaxation of the dry-air
%     pressure update; convergence remains based on the raw Thomas map defect.
%   Change: 2026-07-14 Jizixing Zhu -- match the HYDRUS variably saturated
%     convergence contract: use total-water-equivalent change only where both
%     iterates are unsaturated, use head change for saturated/crossing nodes,
%     and retain active-set changes as diagnostics rather than rejection gates.
%   Change: 2026-07-14 Jizixing Zhu -- add the staged h/U/P_g correction path:
%     complete residual rows, matrix-level T-to-U column transforms, sequential
%     scalar Thomas solves, protected U-to-T inversion, and U-aware Anderson.
%   Change: 2026-07-15 Jizixing Zhu -- make h/U/P_g the sole production
%     coordinate path after the 365-day coupled validation; remove the
%     temporary T-primary driver while retaining reusable scalar kernels.
%   Change: 2026-07-15 Jizixing Zhu -- document that the energy convergence
%     gate uses the unscaled raw Picard defect, independently of line-search damping.
%   Change: 2026-07-15 Jizixing Zhu -- consume one immutable water/heat/dry-air
%     trial boundary request and return solved boundary fluxes for all enabled equations.
%   Change: 2026-07-15 Jizixing Zhu -- recover accepted dry-air boundary
%     fluxes from the complete h/T/P_g correction and exact storage change.
%   Change: 2026-07-16 Jizixing Zhu -- require the exact dry-air column mass
%     residual, in addition to the pressure-map defect, before acceptance.
%   Change: 2026-07-17 Jizixing Zhu -- pass surface-water capacity into gas
%     coefficients so a truly ponded top receives no artificial mobility floor.
%   Change: 2026-07-17 Jizixing Zhu -- resolve atmospheric complementarity
%     inside Picard and require an independently converged surface-water balance.
%   Change: 2026-07-17 Jizixing Zhu -- require both finite-ponding KKT pairs,
%     in addition to the surface reservoir mass identity, before acceptance.
%   Change: 2026-07-17 Jizixing Zhu -- require the exact internal-energy
%     column identity in addition to the U-primary Picard-map defect.
%   Change: 2026-07-17 Jizixing Zhu -- feed the post-damping hcrit projection
%     signal to the typed atmospheric active set while retaining that
%     projection as a nonlinear safeguard.
%   Change: 2026-07-17 Jizixing Zhu -- hold each surface active set fixed
%     until the coupled subsurface equations and conservation ledgers close;
%     only then use its boundary reaction to retain or switch the set.
%   Change: 2026-07-17 Jizixing Zhu -- allow a direct head/supply constraint
%     hit to enter a feasible active set before branch convergence, while
%     still deferring reaction-flux exits until the branch has converged.
%   Change: 2026-07-17 Jizixing Zhu -- recover water boundary reactions from the
%     complete h/T/P_g correction and exact total-water storage change.
%   Change: 2026-07-17 Jizixing Zhu -- restart dry-head transitions from the common
%     substep predictor so a projected dry trial cannot seed its conjugate set.
%   Change: 2026-07-17 Jizixing Zhu -- let typed boundary rows own hcrit; disable the
%     legacy post-solve dry projection that overwrote Neumann solutions.
%   Change: 2026-07-17 Jizixing Zhu -- activate the dry head set only from the
%     converged residual-water endpoint and retain its state-consistent head.
%   Change: 2026-07-17 Jizixing Zhu -- make unconverged surface feasibility probes
%     transactional so a branch switch cannot inherit rejected h/U/P_g state.
%   Change: 2026-07-19 Jizixing Zhu -- replace potential bare-soil latent
%     cooling by the dry active set's realized evaporation inside the same
%     h/U/P_g Picard map, without rerunning the canopy energy balance.
%   Change: 2026-07-20 Jizixing Zhu — share Richards-domain storage and root
%     integration with the accepted outer ledger for MODFLOW boundaries.
%   Change: 2026-08-08 Jizixing Zhu -- when a converged zero-head reaction
%     releases to supply limitation, retain that warm start but defer the
%     reverse head test until the prescribed-flux branch has converged.
%   Change: 2026-08-12 Jizixing Zhu -- keep dry-air density and storage in
%     their nonnegative physical domain with a coupled fraction-to-boundary
%     line search and a final accepted-state gate.
%   Change: 2026-08-13 Jizixing Zhu -- damp a dry-air mass-coordinate release
%     from zero and expose iteration-level audit data for nonsmooth diagnosis.
%   Change: 2026-08-15 Jizixing Zhu -- add a frozen h/U/P_g block
%     Modified-Picard solve with a common residual line search after the
%     ordinary sequential Picard map stalls.
%   Change: 2026-08-15 Jizixing Zhu -- expose the nested solver hierarchy and
%     make one candidate evaluation own all exact ledgers, physical
%     admissibility, and the state/flux consumed by surface complementarity.
%   Change: 2026-08-15 Jizixing Zhu -- scale the dry-air conservation gate by
%     accepted column dry-air mass at a one-part-per-million production
%     default while retaining the original absolute tolerance as a floor.
%   Change: 2026-08-15 Jizixing Zhu -- publish the exact dry-air and energy
%     boundary fluxes used by acceptance, eliminating a frozen-map/outer-
%     ledger mismatch after a candidate has converged.
%   Change: 2026-08-16 Jizixing Zhu -- use the explicit SITE-owned 1e-4
%     dry-air inventory-relative conservation tolerance while retaining the
%     absolute small-inventory floor and all physical admissibility gates.
    built   = model.built;
    wbFlux  = getNumDef(solveCtl, 'water_boundary_flux', built.workspace.flux);   % this-step pending water-boundary forcing (EVAP/Throughfall); default = built.workspace.flux (backward-compat)
    soilBoundaryRequest = getNumDef(solveCtl, 'soil_boundary_request', []);
    boundaryRequest = getNumDef(solveCtl, 'water_boundary_request', []);           % compatibility for direct water-only callers
    heatBoundaryRequest = [];
    airBoundaryRequest = [];
    if ~isempty(soilBoundaryRequest)
        boundaryRequest = soilBoundaryRequest.water;
        heatBoundaryRequest = soilBoundaryRequest.heat;
        airBoundaryRequest = soilBoundaryRequest.dry_air;
    end
    waterRootSink = getNumDef(solveCtl, 'water_root_sink', built.workspace.flux.Srt);   % this-step root sink Srt; default = committed built.workspace.flux.Srt (backward-compat); fixed for all substeps of this call
    prev_dt = model.runtime.substep.Delt_old;                         % [s] A3 previous-substep dt (used only kt>1)

    built = soil.water.extrapolateSubstepGuess(built, kt, dt, prev_dt);   % A3: variable-dt Picard guess
    built = soil.refreshDualPorositySink(built, dt);                      % A5: dual-porosity sink (no-op if lDualPor=0)
    activeSetStart = built;                                                % immutable trial predictor for same-dt branch restarts
    setup       = solver.setupPicardAnderson(built);   % A4 SETUP (I4: now consumed)
    lAnderson   = setup.lAnderson;
    aa          = setup.aaPicard;                      % [] -> filled by the loop (Anderson history)
    aaStateOpts = setup.aaStateOpts;
    aaOpts      = setup.aaOpts;
    nApply      = 0;
    nCommit     = 0;

    converged = false;
    iters     = solveCtl.NIT;
    residual  = NaN;
    waterFlux = struct();                                            % WBC-2: overwritten by each solve iterate (NIT>=1)
    heatFlux = struct('QET', 0, 'QEB', 0);
    airFlux = struct();
    boundaryApplied = struct();
    heatBoundaryApplied = struct();
    heatBoundaryRequestApplied = heatBoundaryRequest;
    airBoundaryApplied = struct();
    CHKHEAT = zeros(built.domain.NN, 1);
    temperatureIncrement = zeros(built.domain.NN, 1);
    CHKAIR = zeros(built.domain.NN, 1);
    heatEnabled = switchEnabled(built.control.switches, 'lThmrlefc');
    airEnabled = switchEnabled(built.control.switches, 'lSoilairefc');
    tolEnergy = getNumDef(solveCtl, 'tol_energy', ...
        getNumDef(built.control.numerics, 'tol_energy', 0.1));          % [degC]
    tolEnergyBalance = getNumDef(solveCtl, 'tol_energy_balance_J_cm2', ...
        getNumDef(built.control.numerics, ...
            'tol_energy_balance_J_cm2', 1e-2));                         % [J cm^-2]
    tolAirRelative = getNumDef(solveCtl, 'tol_air_relative', ...
        getNumDef(built.control.numerics, 'tol_air_relative', 1e-3));   % [-]
    tolAirMass = getNumDef(solveCtl, 'tol_air_mass_g_cm2', ...
        getNumDef(built.control.numerics, 'tol_air_mass_g_cm2', 1e-8)); % [g cm^-2]
    assert(isscalar(tolAirRelative) && isfinite(tolAirRelative) ...
        && tolAirRelative > 0, ...
        'attemptSolveSoilEquations:badTolAirRelative', ...
        'tol_air_relative (%g) must be a finite positive scalar.', ...
        tolAirRelative);
    assert(isscalar(tolAirMass) && isfinite(tolAirMass) ...
            && tolAirMass >= 0, ...
        'attemptSolveSoilEquations:badTolAirMass', ...
        'tol_air_mass_g_cm2 (%g) must be finite and nonnegative.', ...
        tolAirMass);
    assert(isscalar(tolEnergyBalance) && isfinite(tolEnergyBalance) ...
            && tolEnergyBalance >= 0, ...
        'attemptSolveSoilEquations:badTolEnergyBalance', ...
        'tol_energy_balance_J_cm2 (%g) must be finite and nonnegative.', ...
        tolEnergyBalance);
    L_f = getNumDef(built.workspace.step, 'L_f', 0);
    refreshOpts = struct('L_f', L_f, 'dt_s', dt, 'KT', kt, ...
        'NoTime_KT', kt, ...
        'DSTMAX', built.workspace.surface.DSTMAX, ...
        'refreshThermalBoundary', isempty(heatBoundaryRequest));
    % Candidate mass telemetry uses the mixed equation's native storage ruler:
    % FE endpoint condensed plus vapor-equivalent water. The same accepted
    % storage is computed below together with the required thetaEqN array.
    tolMass = getNumDef(solveCtl, 'tol_water_mass_mm', ...
        getNumDef(built.control.numerics, 'tol_water_mass_mm', 1e-4));                   % [mm]
    tolSurfaceWaterCm = getNumDef(solveCtl, 'tol_surface_water_cm', ...
        max(1e-10, tolMass / 10));                                                       % [cm]
    firstNonlinearPassIter = NaN;
    firstNonlinearPassCHK = NaN;
    firstNonlinearPassMass = NaN;
    finalIter = NaN;
    finalCHK = NaN;
    finalMass = NaN;
    tolTheta = getNumDef(solveCtl, 'tol_theta', 1e-4);                    % [-] HYDRUS TolTh
    tolHeadSaturated = getNumDef(solveCtl, 'tol_head_saturated_cm', 0.1); % [cm] HYDRUS TolH
    assert(isscalar(tolTheta) && isfinite(tolTheta) && tolTheta > 0 ...
        && isscalar(tolHeadSaturated) && isfinite(tolHeadSaturated) ...
        && tolHeadSaturated > 0, ...
        'attemptSolveSoilEquations:badHydrusTolerance', ...
        'HYDRUS tolerances must be finite and > 0.');
    assert(isscalar(tolMass) && isfinite(tolMass) && tolMass >= 0, ...
        'attemptSolveSoilEquations:badTolWaterMass', ...
        'tol_water_mass_mm (%g) must be a finite nonnegative scalar.', tolMass);
    assert(isscalar(tolSurfaceWaterCm) && isfinite(tolSurfaceWaterCm) ...
            && tolSurfaceWaterCm >= 0, ...
        'attemptSolveSoilEquations:badTolSurfaceWater', ...
        'tol_surface_water_cm (%g) must be finite and nonnegative.', ...
        tolSurfaceWaterCm);
    terminationReason = 'nit_exhausted';     % default: loop ran to NIT without an early break
    finalThetaIncrement = NaN;                % [-] HYDRUS fully-unsaturated endpoint convergence metric
    finalSaturatedHeadIncrement = NaN;        % [cm] HYDRUS saturated-or-crossing node convergence metric
    finalActiveSetChanged = false;             % [-] h<0 versus h>=0 classification changed
    finalAirRelativeIncrement = 0;              % [-] dry-air Picard increment
    finalAirMassCoordinateIncrement = 0;        % [g cm^-2] conservative nodal map defect
    finalAirMassResidual = NaN;                 % [g cm^-2] exact accepted-state closure
    finalAirAdmissibility = disabledDryAirAdmissibility();
    finalEnergyBalanceResidual = NaN;            % [J cm^-2] exact accepted-state closure
    airRelaxation = initializeAirRelaxation();  % numerical state; reset for every trial/retry
    correctionDiagnostics = struct('max_delta_h_cm', NaN, ...
        'max_delta_p', NaN, 'max_delta_U', NaN, ...
        'max_delta_dry_air_mass_g_cm2', NaN, ...
        'correction_scale', NaN, 'temperature_failure', '', ...
        'dry_air_failure', '');
    hasSurfaceComplementarity = hasAtmosphericSurface(boundaryRequest);
    surfaceDecision = struct();
    finalSurfaceState = struct();
    surfaceConverged = ~hasSurfaceComplementarity;
    surfaceSwitchCount = 0;
    surfaceMassResidualCm = NaN;
    maxSurfaceSwitches = 6;
    surfaceInitialRegime = surfaceRegime(boundaryRequest);
    surfaceTrace = strings(0, 1);
    iterationTrace = repmat(struct( ...
        'iteration', 0, 'air_coordinate', '', ...
        'dry_air_active_node_count', 0, ...
        'max_dry_air_active_reaction', 0, ...
        'max_air_relative_increment', NaN, ...
        'max_air_relative_node', NaN, ...
        'saturated_head_increment_cm', NaN, ...
        'energy_increment_C', NaN), solveCtl.NIT, 1);
    blockPicardRequested = false;
    % Evaluate the accepted substep-start storage exactly once. The accepted
    % Y state, not the extrapolated Picard trial, is the backward-Euler time
    % anchor. Reusing thetaEqN through all Picard iterations is both physically
    % required and avoids an accepted-state refresh inside the iteration loop.
    exactContext = prepareExactResidualContext( ...
        built, dt, kt, wbFlux, waterRootSink, soilBoundaryRequest, ...
        boundaryRequest, heatBoundaryRequest, airBoundaryRequest, ...
        refreshOpts);
    thetaEqN = exactContext.theta_equivalent_accepted;
    dryAirStorageN = exactContext.dry_air_storage_accepted;
    energyStorageN = exactContext.energy_storage_accepted;
    energyCapacityN = exactContext.energy_capacity_accepted;
    % Dry-air residuals scale with the conserved column inventory. Keep the
    % original absolute gate as the small-inventory floor and use a 100 ppm
    % relative production default. This prevents a normal 0.29 g cm^-2 air
    % column from being forced to an inventory-relative error near 3e-8,
    % which otherwise creates nonlinear retries and sub-hour shrinkage even
    % on the stable orchard case. Every outcome reports both terms and the
    % effective gate; callers may explicitly set zero for strict diagnostics.
    tolAirMassRelative = getNumDef(solveCtl, ...
        'tol_air_mass_relative', getNumDef(built.control.numerics, ...
            'tol_air_mass_relative', 1e-4));
    assert(isscalar(tolAirMassRelative) && isfinite(tolAirMassRelative) ...
            && tolAirMassRelative >= 0, ...
        'attemptSolveSoilEquations:badTolAirMassRelative', ...
        'tol_air_mass_relative (%g) must be finite and nonnegative.', ...
        tolAirMassRelative);
    if airEnabled
        dryAirColumnMassScale = sum(lumpEndpointState( ...
            dryAirStorageN, built.domain.dz));
        tolAirMassGate = max(tolAirMass, ...
            tolAirMassRelative * max(dryAirColumnMassScale, realmin));
    else
        dryAirColumnMassScale = NaN;
        tolAirMassGate = tolAirMass;
    end
    exactContext.scales.water_cm_s = max(tolMass / (10 * dt), realmin);
    exactContext.scales.dry_air_g_cm2_s = ...
        max(tolAirMassGate / dt, realmin);
    exactContext.scales.energy_J_cm2_s = ...
        max(tolEnergyBalance / dt, realmin);
    exactContext.scale_id = ...
        'compass.soil.coupled_acceptance_scaled_residual.v1';
    aaStateOpts.includeThermal = heatEnabled;
    aaStateOpts.includeEnergy = heatEnabled;
    aaStateOpts.includePressure = airEnabled;
    if heatEnabled
        aaStateOpts.scale_UU = max(abs(energyCapacityN) * 10, 1e-8);
    end
    assert(all(isfinite(thetaEqN(:))) ...
        && isequal(size(thetaEqN), size(built.workspace.derived.Theta_LL)) && dt > 0, ...
        'attemptSolveSoilEquations:badMixedInputs', ...
        'thetaEqN must be finite NL x 2 and dt must be positive.');
    for KIT = 1:solveCtl.NIT
        if heatEnabled && airEnabled && KIT == 20 ...
                && ~blockPicardRequested
            blockPicardRequested = true;
            aa = [];
        end
        % I4 Anderson APPLY (legacy COMPASS.m:1665-1675): KIT>1 accelerated guess BEFORE refresh.
        if lAnderson
            stateBeforeAnderson = built;
            [xIn, meta] = solver.packSoilPicardState( ...
                picardStateFromTrial(built, heatEnabled), ...
                aaStateOpts);
            if KIT > 1
                [xGuess, aa] = solver.andersonApplyGuess(xIn, aa, aaOpts);
                if isfield(aa, 'lastApplied') && aa.lastApplied, nApply = nApply + 1; end
                built = writeTrialFromPicardState(built, ...
                    solver.unpackSoilPicardState( ...
                        picardStateFromTrial(built, heatEnabled), ...
                        xGuess, meta));
                if heatEnabled
                    [built.workspace.step.trial, inversionGuess] = recoverTrialTemperature( ...
                        built, refreshOpts);
                    if ~inversionGuess.allConverged
                        built = stateBeforeAnderson;
                        aa = [];
                    end
                end
                if airEnabled
                    [airGuess, airGuessEquations] = ...
                        soil.refreshSoilCoefficients(built, refreshOpts);
                    guessAdmissibility = soil.gas.evaluateDryAirAdmissibility( ...
                        airGuess.workspace.derived.RHODA, ...
                        airGuessEquations.dry_air.dryAirStorage);
                    if guessAdmissibility.admissible
                        built = airGuess;
                    else
                        built = stateBeforeAnderson;
                        aa = [];
                    end
                end
                [xIn, ~] = solver.packSoilPicardState( ...
                    picardStateFromTrial(built, heatEnabled), ...
                    aaStateOpts);   % re-pack guessed state (commit residual = solved - guessed)
            end
        end

        if heatEnabled
            [built.workspace.step.trial, inversionAtIteration] = recoverTrialTemperature( ...
                built, refreshOpts);
            if ~inversionAtIteration.allConverged
                terminationReason = 'temperature_inversion_failed';
                break;
            end
        end

        % Layer 1-2: refresh constitutive state and build the frozen shared
        % coefficient/flux bank for the current Picard base state.
        [built, equations] = soil.refreshSoilCoefficients(built, refreshOpts);
        if airEnabled
            iterationAirStatus = soil.gas.evaluateDryAirAdmissibility( ...
                built.workspace.derived.RHODA, ...
                equations.dry_air.dryAirStorage);
            if ~iterationAirStatus.admissible
                finalAirAdmissibility = iterationAirStatus;
                terminationReason = 'dry_air_admissibility_failed';
                break;
            end
        end
        matrixWaterSink = effectiveWaterEquationSink(built, waterRootSink);
        % Layer 3: assemble discrete rows without solving or accepting.
        matrices = soil.assembleEquationMatrices( ...
            equations, built.domain.dz, matrixWaterSink);
        hIteration = built.workspace.step.trial.h;
        dIteration = built.workspace.derived;
        thetaIteration = waterEquivalentEndpoints( ...
            dIteration.Theta_UU, dIteration.Theta_V, dIteration.RHOV);
        % Layer 4: apply the ordinary or block Modified-Picard map.
        try
            [built, mapResult, airRelaxation] = solveInternalEnergyPrimaryMap( ...
                built, equations, matrices, thetaEqN, thetaIteration, dryAirStorageN, ...
                energyStorageN, dt, kt, wbFlux, waterRootSink, ...
                boundaryRequest, heatBoundaryRequest, airBoundaryRequest, ...
                refreshOpts, airEnabled, heatEnabled, ...
                tolAirRelative, airRelaxation, exactContext, ...
                blockPicardRequested);
        catch err
            if strcmp(err.identifier, ...
                    'solveLaggedPicardBlock:singularBlockPivot')
                % A singular frozen Picard pivot invalidates only this trial.
                % Return a normal rejection so the outer scheduler can roll
                % back the snapshot and retry with a smaller substep.
                terminationReason = 'singular_block_pivot';
                break;
            end
            rethrow(err);
        end
        CHKWATER = mapResult.CHKWATER;
        CHKAIR = mapResult.CHKAIR;
        CHKHEAT = mapResult.CHKHEAT;
        temperatureIncrement = mapResult.temperatureIncrement;
        waterFlux = mapResult.waterFlux;
        heatFlux = mapResult.heatFlux;
        airFlux = mapResult.airFlux;
        boundaryApplied = mapResult.boundaryApplied;
        heatBoundaryApplied = mapResult.heatBoundaryApplied;
        heatBoundaryRequestApplied = mapResult.heatBoundaryRequest;
        airBoundaryApplied = mapResult.airBoundaryApplied;
        finalAirRelativeIncrement = mapResult.airRelativeIncrement;
        finalAirMassCoordinateIncrement = ...
            mapResult.airMassCoordinateIncrement;
        correctionDiagnostics = mapResult.diagnostics;
        iterationTrace(KIT) = struct( ...
            'iteration', KIT, ...
            'air_coordinate', mapResult.airCoordinate, ...
            'dry_air_active_node_count', ...
                correctionDiagnostics.dry_air_active_node_count, ...
            'max_dry_air_active_reaction', ...
                correctionDiagnostics.max_dry_air_active_reaction, ...
            'max_air_relative_increment', ...
                mapResult.airRelativeIncrement, ...
            'max_air_relative_node', ...
                correctionDiagnostics.max_air_relative_node, ...
            'saturated_head_increment_cm', NaN, ...
            'energy_increment_C', max(mapResult.CHKHEAT(:)));
        if ~mapResult.correctionAccepted
            if ~mapResult.temperatureInversion.allConverged
                terminationReason = 'temperature_inversion_failed';
            else
                terminationReason = 'dry_air_admissibility_failed';
            end
            finalAirAdmissibility = mapResult.airAdmissibility;
            break;
        end
        residual = max(CHKWATER(:));                                     % [cm] legacy Error_Water = max(CHKWATER)
        [nonlinearPass, thetaIncrement, saturatedHeadIncrement, activeSetChanged] = ...
            evaluateHydrusConvergence(built, hIteration, thetaIteration, ...
                tolTheta, tolHeadSaturated, refreshOpts, ...
                modflowFirstActiveNode(boundaryRequest, built.domain.NN));
        finalThetaIncrement = thetaIncrement;
        finalSaturatedHeadIncrement = saturatedHeadIncrement;
        finalActiveSetChanged = activeSetChanged;
        iterationTrace(KIT).saturated_head_increment_cm = ...
            saturatedHeadIncrement;
        heatPass = ~heatEnabled || max(CHKHEAT(:)) < tolEnergy;
        airPass = ~airEnabled || finalAirRelativeIncrement < tolAirRelative;
        subsurfaceConverged = false;
        if nonlinearPass && heatPass && airPass
            % Layer 5: one exact refresh supplies every acceptance ledger
            % and the surface interface state. No legacy per-equation scorer
            % is allowed to refresh this same candidate again.
            exactCandidate = reuseExactCandidateForState( ...
                built, mapResult, exactContext, boundaryRequest, ...
                mapResult.heatBoundaryRequest);
            finalSurfaceState = exactCandidate.surface_state;
            massR = exactCandidate.ledgers.water.residual_mm;
            waterFlux = exactCandidate.ledgers.water_boundary_flux;
            if airEnabled
                airFlux = exactCandidate.ledgers.dry_air_boundary_flux;
            end
            if heatEnabled
                heatFlux = exactCandidate.ledgers.energy_boundary_flux;
            end
            if isnan(firstNonlinearPassIter)
                firstNonlinearPassIter = KIT;
                firstNonlinearPassCHK = residual;
                firstNonlinearPassMass = massR;
            end
            finalIter = KIT;
            finalCHK = residual;
            finalMass = massR;
            if airEnabled
                finalAirMassResidual = ...
                    exactCandidate.ledgers.dry_air.residual_g_cm2;
                finalAirAdmissibility = ...
                    exactCandidate.feasibility.dry_air;
                airMassPass = abs(finalAirMassResidual) < tolAirMassGate;
                airPhysicalPass = finalAirAdmissibility.admissible;
            else
                airMassPass = true;
                airPhysicalPass = true;
            end
            if heatEnabled
                finalEnergyBalanceResidual = ...
                    exactCandidate.ledgers.energy.residual_J_cm2;
                energyBalancePass = isfinite(finalEnergyBalanceResidual) ...
                    && abs(finalEnergyBalanceResidual) ...
                    < tolEnergyBalance;
            else
                energyBalancePass = true;
            end
            subsurfaceConverged = abs(massR) < tolMass ...
                && airMassPass && airPhysicalPass && energyBalancePass;
            % A failed exact ledger is not, by itself, evidence that the
            % ordinary Picard map is stalled. Keep the inexpensive scalar
            % map through iteration 19; the deterministic KIT=20 trigger at
            % the top of the loop then resets Anderson history and activates
            % the block map. This preserves the exact ledger as acceptance
            % authority without paying for a block line search on every
            % ordinary step.
        end

        % A boundary reaction is physically meaningful only after the
        % currently selected water/heat/dry-air branch and its conservation
        % ledgers have converged.  Holding the active set fixed until then is
        % the 1-D nonlinear active-set analogue of a complementarity solve;
        % it prevents an intermediate Picard reaction from toggling dry/head
        % and flux rows.  Abrupt wetting is still handled before this loop by
        % the wet-feasibility predictor, so an extreme rain pulse need not
        % first converge on an infeasible dry flux branch.
        if hasSurfaceComplementarity
            surfaceConverged = false;
            immediateConstraintHit = shouldCheckSurfaceConstraintNow( ...
                boundaryRequest, built.workspace.step.trial.h(end), ...
                mapResult.surfaceDryHeadProjected, tolSurfaceWaterCm);
            if subsurfaceConverged || immediateConstraintHit
                exactCandidate = reuseExactCandidateForState( ...
                    built, mapResult, exactContext, boundaryRequest, ...
                    mapResult.heatBoundaryRequest);
                finalSurfaceState = exactCandidate.surface_state;
                waterFlux = exactCandidate.ledgers.water_boundary_flux;
                if airEnabled
                    airFlux = exactCandidate.ledgers.dry_air_boundary_flux;
                end
                if heatEnabled
                    heatFlux = exactCandidate.ledgers.energy_boundary_flux;
                end
                surfaceDecision = soil.water.boundary.evaluateSurfaceTrial( ...
                    boundaryRequest, exactCandidate.surface_state, ...
                    tolSurfaceWaterCm, struct( ...
                        'dry_head_projected', ...
                        mapResult.surfaceDryHeadProjected, ...
                        'dry_head_target_cm', ...
                        built.workspace.step.trial.h(end)));
                surfaceTrace(end+1, 1) = formatSurfaceTrace( ...
                    surfaceDecision, KIT); %#ok<AGROW>
                if ~surfaceDecision.valid
                    surfaceSwitchCount = surfaceSwitchCount + 1;
                    surfaceMassResidualCm = NaN;
                    if surfaceSwitchCount > maxSurfaceSwitches
                        terminationReason = 'surface_active_set_cycle';
                        break;
                    end

                    % The nonlinear operator changes with the top boundary
                    % row. An unconverged feasibility probe is disposable: a
                    % switch must restart from the immutable same-dt predictor
                    % rather than seed the conjugate set with rejected h/U/P_g
                    % state. A converged wet branch remains a warm start, but
                    % its released supply branch must converge before the
                    % reverse head test is physically meaningful. Every switch
                    % discards histories and reassembles all h/U/P_g rows.
                    previousRegime = char( ...
                        boundaryRequest.top.surface.regime);
                    nextRegime = char(surfaceDecision.next_regime);
                    deferSupplyHeadTest = subsurfaceConverged ...
                        && strcmp(previousRegime, 'capacity_limited') ...
                        && strcmp(nextRegime, 'supply_limited');
                    boundaryRequest.top.surface.regime = ...
                        nextRegime;
                    if strcmp(nextRegime, 'dry_evaporation_limited') ...
                            && isfinite( ...
                                surfaceDecision.next_dry_head_limit_cm)
                        boundaryRequest.top.surface.dry_head_limit_cm = ...
                            surfaceDecision.next_dry_head_limit_cm;
                    end
                    if ~subsurfaceConverged ...
                            || strcmp(previousRegime, ...
                                'dry_evaporation_limited') ...
                            || strcmp(nextRegime, ...
                                'dry_evaporation_limited')
                        built = activeSetStart;
                    end
                    if deferSupplyHeadTest
                        boundaryRequest.top.surface.initialization = ...
                            'converged_reaction_release';
                    else
                        boundaryRequest.top.surface.initialization = ...
                            'explicit_active_set';
                    end
                    aa = [];
                    airRelaxation = initializeAirRelaxation();
                    continue;
                end
                if subsurfaceConverged
                    surfaceMassResidualCm = surfaceDecision.residual_cm;
                    surfaceConverged = surfaceDecision.complementarity_valid ...
                        && abs(surfaceMassResidualCm) <= tolSurfaceWaterCm;
                end
            end
        end
        if subsurfaceConverged && surfaceConverged
            converged = true;
            iters = KIT;
            terminationReason = 'accepted';
        end

        % I4 Anderson COMMIT (legacy COMPASS.m:1883-1887): store the solved map EVERY iterate (incl. the
        % converged one); the break is AFTER commit so the converged map is never skipped.
        if lAnderson
            [xMap, ~] = solver.packSoilPicardState( ...
                picardStateFromTrial(built, heatEnabled), ...
                aaStateOpts);
            aa = solver.andersonCommitMap(xIn, xMap, aa, aaOpts);        % no nonlinear boundary reset in this scalar path
            nCommit = nCommit + 1;
        end

        if converged, break; end
    end

    if converged
        % Consistency: leave derived.Theta_LL == Theta_LL(accepted trial.h) so carryThetaLayer (mid-tile
        % accept + outer commit) writes a Y.Theta_L that is VG-consistent with the committed Y.h. The
        % loop's last refresh was PRE-solve (top of the accepting iterate); this final refresh is
        % POST-solve. Only on accept -- a rejected trial must not carry a refreshed candidate.
        [built, eqsAccept] = soil.refreshSoilCoefficients(built, refreshOpts);
        if isfield(eqsAccept, 'heat') && isfield(eqsAccept.heat, 'energyStorage') ...
                && ~isempty(eqsAccept.heat.energyStorage)
            acceptedEnergyStorage = eqsAccept.heat.energyStorage;
            acceptedRootEnergySink = eqsAccept.heat.rootEnergySink;
        else
            acceptedEnergyStorage = [];
            acceptedRootEnergySink = [];
        end
        if airEnabled
            assert(isfield(eqsAccept, 'dry_air') ...
                    && isfield(eqsAccept.dry_air, 'dryAirStorage') ...
                    && ~isempty(eqsAccept.dry_air.dryAirStorage) ...
                    && isfield(mapResult, 'airFluxContext') ...
                    && ~isempty(fieldnames(mapResult.airFluxContext)), ...
                'attemptSolveSoilEquations:missingAirFluxClosure', ...
                ['An accepted dry-air solve requires final storage and the ', ...
                 'last complete boundary-reaction context.']);
            airContext = mapResult.airFluxContext;
            [QAT, QAB] = soil.gas.calculateCoupledBoundaryFlux( ...
                airContext.SAVE, airContext.matrix, ...
                airContext.stateBefore, airContext.stateAfter, ...
                airContext.storageBefore, ...
                eqsAccept.dry_air.dryAirStorage, built.domain.dz, dt);
            [QAT, QAB] = soil.gas.applyActiveSetBoundaryReaction( ...
                QAT, QAB, airContext.active_reaction, ...
                airContext.active_nodes, airBoundaryApplied);
            airFlux = struct('QAT', QAT, 'QAB', QAB);
            acceptedDryAirStorage = eqsAccept.dry_air.dryAirStorage;
            finalAirAdmissibility = soil.gas.evaluateDryAirAdmissibility( ...
                built.workspace.derived.RHODA, acceptedDryAirStorage);
            assert(finalAirAdmissibility.admissible, ...
                'attemptSolveSoilEquations:inadmissibleAcceptedDryAir', ...
                ['A converged dry-air candidate crossed the nonnegative ', ...
                 'density/storage boundary before acceptance.']);
        else
            acceptedDryAirStorage = [];
        end
        acceptedSolverDiagnostics = ...
            collectSolverDiagnostics(eqsAccept, built);
    else
        acceptedEnergyStorage = [];
        acceptedRootEnergySink = [];
        acceptedDryAirStorage = [];
        acceptedSolverDiagnostics = struct();
    end
    model.built = built;
    outcome = struct('converged', converged, 'residual', residual, 'iters', iters, ...
        'energyStorage', acceptedEnergyStorage, ...
        'rootEnergySink', acceptedRootEnergySink, ...
        'dryAirStorageBefore', dryAirStorageN, ...
        'dryAirStorageAfter', acceptedDryAirStorage);
    outcome.waterFlux = waterFlux;   % WBC-2: the converged tile's water-boundary flux (for the accumulator)
    outcome.waterBoundaryRequest = boundaryRequest;   % W1b audit metadata only; never committed to Y/workspace
    outcome.waterBoundaryApplied = boundaryApplied;   % W1b active mathematical BC kind; especially atmospheric head/flux switch
    outcome.heatBoundaryApplied = heatBoundaryApplied;
    outcome.heatBoundaryRequest = heatBoundaryRequestApplied;
    outcome.airBoundaryApplied = airBoundaryApplied;
    outcome.anderson = struct('lAnderson', lAnderson, 'applied', nApply, 'committed', nCommit);   % I4 stats
    outcome.first_nonlinear_pass_iter = firstNonlinearPassIter;                  % [-]
    outcome.first_nonlinear_pass_CHKWATER = firstNonlinearPassCHK;               % [cm], diagnostic only
    outcome.first_nonlinear_pass_mass_residual_mm = firstNonlinearPassMass;       % [mm]
    outcome.final_iters = finalIter;                                               % [-]
    outcome.final_CHKWATER = finalCHK;                                             % [cm], diagnostic only
    outcome.final_mass_residual_mm = finalMass;                                   % [mm]
    outcome.nonlinear_converged = ~isnan(firstNonlinearPassIter);                  % [-]
    outcome.mass_converged = ~isnan(finalMass) && abs(finalMass) < tolMass;        % [-]
    outcome.termination_reason = terminationReason;                               % [-]
    outcome.final_theta_increment = finalThetaIncrement;                          % [-]
    outcome.final_saturated_head_increment_cm = finalSaturatedHeadIncrement;     % [cm]
    outcome.final_active_set_changed = finalActiveSetChanged;                    % [-]
    outcome.surfaceWater = surfaceDecision;
    outcome.surfaceState = finalSurfaceState;
    outcome.surface_converged = surfaceConverged;                                % [-]
    outcome.surface_switch_count = surfaceSwitchCount;                           % [-]
    outcome.surface_initial_regime = surfaceInitialRegime;                       % [-]
    outcome.surface_final_regime = surfaceRegime(boundaryRequest);               % [-]
    outcome.surface_trace = strjoin(surfaceTrace, '|');                          % audit-only
    outcome.final_surface_mass_residual_cm = surfaceMassResidualCm;              % [cm]
    outcome.tol_surface_water_cm = tolSurfaceWaterCm;                            % [cm]
    outcome.heatFlux = heatFlux;                                                  % [J cm-2 s-1]
    outcome.airFlux = airFlux;                                                     % [g cm-2 s-1]
    outcome.final_heat_increment_C = max(temperatureIncrement(:));                % [degC] actual T change
    outcome.final_energy_equivalent_increment_C = max(CHKHEAT(:));                % [degC] unscaled raw |dU|/U_T Picard defect
    outcome.final_air_pressure_increment = max(CHKAIR(:));                        % [internal pressure]
    outcome.final_air_mass_coordinate_increment_g_cm2 = ...
        finalAirMassCoordinateIncrement;                                          % [g cm^-2]
    outcome.final_air_relative_increment = finalAirRelativeIncrement;             % [-]
    outcome.final_air_mass_residual_g_cm2 = finalAirMassResidual;                 % [g cm^-2]
    outcome.tol_air_mass_absolute_g_cm2 = tolAirMass;                             % [g cm^-2]
    outcome.tol_air_mass_relative = tolAirMassRelative;                           % [-]
    outcome.tol_air_mass_effective_g_cm2 = tolAirMassGate;                        % [g cm^-2]
    outcome.dry_air_column_mass_scale_g_cm2 = dryAirColumnMassScale;               % [g cm^-2]
    outcome.dry_air_admissible = finalAirAdmissibility.admissible;                % [-]
    outcome.minimum_dry_air_density_g_cm3 = ...
        finalAirAdmissibility.minimum_density_g_cm3;                              % [g cm^-3]
    outcome.minimum_dry_air_storage_g_cm3 = ...
        finalAirAdmissibility.minimum_storage_g_cm3;                              % [g cm^-3]
    outcome.dry_air_admissibility_failure = ...
        finalAirAdmissibility.failure_reason;                                     % [-]
    outcome.final_energy_balance_residual_J_cm2 = ...
        finalEnergyBalanceResidual;                                               % [J cm^-2]
    outcome.tol_energy_balance_J_cm2 = tolEnergyBalance;                          % [J cm^-2]
    outcome.energy_balance_converged = ~heatEnabled ...
        || (isfinite(outcome.final_energy_balance_residual_J_cm2) ...
            && abs(outcome.final_energy_balance_residual_J_cm2) ...
            < tolEnergyBalance);
    outcome.heat_converged = ~heatEnabled ...
        || (outcome.final_energy_equivalent_increment_C < tolEnergy ...
            && outcome.energy_balance_converged);
    outcome.air_converged = ~airEnabled ...
        || (outcome.final_air_relative_increment < tolAirRelative ...
            && abs(outcome.final_air_mass_residual_g_cm2) < tolAirMassGate ...
            && outcome.dry_air_admissible);
    outcome.air_relaxation = struct( ...
        'final_omega', airRelaxation.omega, ...
        'minimum_omega', airRelaxation.minimumOmega, ...
        'applied', airRelaxation.applied);
    outcome.internal_energy_correction = correctionDiagnostics;
    outcome.iteration_trace = iterationTrace(1:max(1, min(solveCtl.NIT, KIT)));
    outcome.solverDiagnostics = acceptedSolverDiagnostics;
end

function diagnosticsOut = collectSolverDiagnostics(equations, built)
% Retain only accepted per-layer flux diagnostics needed by model output.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    diagnosticsOut = struct();
    source = struct();
    if isfield(equations, 'heat') && isstruct(equations.heat) ...
            && isfield(equations.heat, 'QL')
        source = equations.heat;
    elseif isfield(equations, 'dry_air') ...
            && isstruct(equations.dry_air) ...
            && isfield(equations.dry_air, 'QL')
        source = equations.dry_air;
    end
    if isempty(fieldnames(source))
        return;
    end
    copyNames = {'QL','QL_h','QL_T','QL_a','QVH','QVT','QVa','Qa', ...
        'DhDZ','DTDZ','DPgDZ','KLhBAR','KLTBAR','DTDBAR'};
    for i = 1:numel(copyNames)
        name = copyNames{i};
        if isfield(source, name)
            diagnosticsOut.(name) = source.(name);
        end
    end
    derived = built.workspace.derived;
    if isfield(derived, 'Se')
        diagnosticsOut.Se = mean(derived.Se, 2);
    end
    if isfield(derived, 'Sa')
        diagnosticsOut.Sa = mean(derived.Sa, 2);
    end
end

% -----------------------------------------------------------------------------
function state = initializeAirRelaxation()
% Initialize one trial-local vector-Aitken dry-air relaxation state.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    state = struct('omega', 0.5, 'previousCorrection', [], ...
        'minimumOmega', 1, 'applied', 0, 'coordinate', '', ...
        'activeNodes', []);
end

% -----------------------------------------------------------------------------
function [coordinateNext, state] = relaxDryAirCorrection( ...
        coordinateOld, correction, state, rawConverged, coordinateName)
% Relax the dry-air fixed-point correction without changing its fixed point.
%
%   The raw Thomas map is xMap=xOld+correction. For an
%   oscillatory distributed-Picard mode, vector Aitken delta-squared chooses
%
%     omega_m = -omega_(m-1) * r_(m-1)'*(r_m-r_(m-1))
%               / ||r_m-r_(m-1)||^2.
%
%   Bounding omega to [0.05,1] prevents a noisy multi-mode correction from
%   reversing or excessively shrinking the update. Pressure and conservative
%   dry-air mass use separate histories because their units and conditioning
%   differ. Once the raw map defect satisfies its coordinate tolerance, the
%   full Thomas solution is written so boundary values are not left partially
%   relaxed.
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-08-12 Jizixing Zhu -- support pressure/mass hybrid coordinates
%           without carrying Aitken history across a coordinate switch.
    correction = correction(:);
    coordinateOld = coordinateOld(:);
    if ~strcmp(state.coordinate, coordinateName)
        state.omega = 0.5;
        state.previousCorrection = [];
        state.coordinate = coordinateName;
        state.activeNodes = [];
    end
    omega = state.omega;

    if rawConverged
        omega = 1;
    elseif ~isempty(state.previousCorrection)
        delta = correction - state.previousCorrection;
        denominator = dot(delta, delta);
        scale = max([dot(correction, correction), ...
            dot(state.previousCorrection, state.previousCorrection), 1]);
        if denominator > eps * scale
            candidate = -state.omega * ...
                dot(state.previousCorrection, delta) / denominator;
            if isfinite(candidate) && candidate > 0
                omega = min(1, max(0.05, candidate));
            else
                omega = max(0.05, 0.5 * state.omega);
            end
        end
    end

    coordinateNext = coordinateOld + omega * correction;
    state.previousCorrection = correction;
    state.omega = omega;
    state.minimumOmega = min(state.minimumOmega, omega);
    state.applied = state.applied + double(omega < 1);
end

% -----------------------------------------------------------------------------
function [built, result, airRelaxation] = solveInternalEnergyPrimaryMap( ...
        built, equations, matrices, thetaEqN, waterStorageM, dryAirStorageN, ...
        energyStorageN, dt_s, kt, wbFlux, waterRootSink, ...
        boundaryRequest, heatBoundaryRequest, airBoundaryRequest, ...
        refreshOpts, airEnabled, heatEnabled, ...
        tolAirRelative, airRelaxation, residualContext, ...
        blockPicardRequested)
% Solve one distributed frozen-coefficient h/U/P_g Picard map.
%
%   All residuals and all h/T/P_g Jacobian columns are assembled at the same
%   current Picard state.  The T columns are then changed to U coordinates at
%   matrix level. Normal states retain the validated P_g coordinate and the
%   stalled path solves their three rows as one block-tridiagonal
%   Modified-Picard map. If the
%   frozen pressure map predicts P_g below vapor partial pressure, all three
%   rows are changed to conservative nodal dry-air mass M_a before the actual
%   correction sweep. After the sequential gas/U updates, the frozen water row is
%   back-substituted once so J_wP*deltaP and J_wU*deltaU cannot leak into the
%   proposed water residual. The dry-air-mass chart uses the same frozen
%   block-Picard contract.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-08-12 Jizixing Zhu -- add one frozen-Jacobian water-row
%           back substitution after the P_g/U corrections.
%   Change: 2026-08-12 Jizixing Zhu -- solve the coupled linear map in
%           conservative nodal dry-air-mass coordinates when a pressure-map
%           predictor would cross the dry-air component boundary; retain the
%           pressure coordinate elsewhere.
%   Change: 2026-08-15 Jizixing Zhu -- solve the stalled ordinary P_g chart
%           with one frozen block Modified-Picard map and one common
%           exact-residual globalization factor.

    NN = built.domain.NN;
    tr = built.workspace.step.trial;
    hCurrent = tr.h(:);
    TCurrent = tr.T(:);
    PCurrent = tr.p_gas(:);

    % Assemble every equation at q^m without accepting any absolute-state map.
    assemblyBuilt = built;
    [assemblyBuilt, ~, waterFlux, boundaryApplied, systemWater] = ...
        soil.water.solveWaterTrial(assemblyBuilt, dt_s, kt, wbFlux, ...
            thetaEqN, waterRootSink, boundaryRequest, equations.water, ...
            matrices.water, true);
    systems = struct('water', systemWater);
    airBoundaryApplied = struct();
    if airEnabled
        [assemblyBuilt, ~, ~, systemAir, airBoundaryApplied] = ...
            soil.gas.solveDryAirTrial( ...
            assemblyBuilt, equations.dry_air, matrices.dry_air, ...
            dryAirStorageN, dt_s, kt, true, airBoundaryRequest);
        systems.air = systemAir;
    else
        systems.air = struct();
    end
    deltaU = zeros(NN,1);
    heatBoundaryApplied = struct();
    heatRequestForMap = heatBoundaryRequest;
    if heatEnabled
        heatRequestForMap = applyRealizedSurfaceEvaporation( ...
            heatBoundaryRequest, boundaryRequest, waterFlux);
        [assemblyBuilt, ~, ~, systemHeat, heatBoundaryApplied] = ...
            soil.thermal.solveThermalTrial( ...
            assemblyBuilt, equations.heat, matrices.heat, energyStorageN, ...
            dt_s, kt, waterFlux, wbFlux, true, heatRequestForMap);
        systems.heat = systemHeat;
    else
        systems.heat = struct();
    end
    built.workspace.solver.SAVE = assemblyBuilt.workspace.solver.SAVE;

    state = struct('h', hCurrent, 'T', TCurrent, 'P', PCurrent);
    rows = solver.buildSoilCorrectionRows( ...
        matrices, systems, state, dt_s, heatEnabled, airEnabled);

    dryMassCurrent = zeros(NN,1);
    dryMass_h = zeros(NN,1);
    dryMass_T = zeros(NN,1);
    dryMass_P = ones(NN,1);
    if airEnabled
        dryMassCurrent = lumpEndpointState( ...
            equations.dry_air.dryAirStorage, built.domain.dz);
        dryMass_h = matrices.dry_air.C1(:,1);
        dryMass_T = matrices.dry_air.C2(:,1);
        dryMass_P = matrices.dry_air.C3(:,1);
        assert(all(isfinite(dryMassCurrent)) && all(dryMassCurrent >= 0) ...
                && all(isfinite(dryMass_P)) && all(dryMass_P > 0), ...
            'attemptSolveSoilEquations:badDryAirMassCoordinate', ...
            ['Nodal dry-air mass must be nonnegative and its pressure ', ...
             'tangent must be finite and positive.']);
        waterMassRow = solver.transformDryAirStorageCoordinate( ...
            rows.water, dryMass_h, dryMass_T, dryMass_P);
        airMassRow = solver.transformDryAirStorageCoordinate( ...
            rows.air, dryMass_h, dryMass_T, dryMass_P);
        if heatEnabled
            heatMassRow = solver.transformDryAirStorageCoordinate( ...
                rows.heat, dryMass_h, dryMass_T, dryMass_P);
        end
    end

    if heatEnabled
        [~, U_T_P, U_h_P, U_P] = soil.thermal.assembleNodalInternalEnergy( ...
            equations.heat.energyStorage, equations.heat.CTh, ...
            equations.heat.CTT, equations.heat.CTa, built.domain.dz);
        waterPressureRow = solver.transformTemperatureEquationColumns( ...
            rows.water, U_h_P, U_T_P, U_P);
        heatPressureRow = solver.transformTemperatureEquationColumns( ...
            rows.heat, U_h_P, U_T_P, U_P);
        if airEnabled
            airPressureRow = solver.transformTemperatureEquationColumns( ...
                rows.air, U_h_P, U_T_P, U_P);
            p_h_at_mass = -dryMass_h ./ dryMass_P;
            p_T_at_mass = -dryMass_T ./ dryMass_P;
            p_M = 1 ./ dryMass_P;
            U_h = U_h_P + U_P .* p_h_at_mass;
            U_T = U_T_P + U_P .* p_T_at_mass;
            U_M = U_P .* p_M;
            waterRow = solver.transformTemperatureMassEquationColumns( ...
                waterMassRow, U_h, U_T, U_M);
            airRow = solver.transformTemperatureMassEquationColumns( ...
                airMassRow, U_h, U_T, U_M);
            heatRow = solver.transformTemperatureMassEquationColumns( ...
                heatMassRow, U_h, U_T, U_M);
        else
            waterRow = waterPressureRow;
            heatRow = heatPressureRow;
            U_h = U_h_P;
            U_T = U_T_P;
            U_M = zeros(NN,1);
        end
    else
        U_T = ones(NN,1);
        U_h = zeros(NN,1);
        U_M = zeros(NN,1);
        if airEnabled
            waterPressureRow = rows.water;
            airPressureRow = rows.air;
            waterRow = waterMassRow;
            airRow = airMassRow;
        else
            waterRow = rows.water;
        end
    end

    useDryMassCoordinate = false;
    dryMassSeedActive = false(NN,1);
    if airEnabled
        assert(all(isfinite(built.workspace.derived.RHODA(:))) ...
                && all(isfinite(built.workspace.derived.Xaa(:))) ...
                && all(built.workspace.derived.Xaa(:) > 0), ...
            'attemptSolveSoilEquations:badDryAirCoordinateProbe', ...
            ['Dry-air density and its pressure derivative must be finite, ', ...
             'with a positive pressure derivative.']);
        dryPartialPressure = built.workspace.derived.RHODA(:) ...
            ./ built.workspace.derived.Xaa(:);
        vaporPartialPressure = PCurrent - dryPartialPressure;
        dryMassSeedActive = dryPartialPressure ...
            <= 1e-10 * max(abs(PCurrent), 1);
        useDryMassCoordinate = strcmp( ...
            airRelaxation.coordinate, 'dry_air_mass') ...
            || any(dryPartialPressure ...
                <= 0.25 * max(PCurrent, realmin));
        if ~useDryMassCoordinate
            probeHRaw = solveBank( ...
                waterPressureRow.h, waterPressureRow.rhs, ...
                [modflowFirstActiveNode(boundaryRequest, NN), NN]);
            probeH = soil.water.applyHeadUpdateSafeguards( ...
                hCurrent, hCurrent + probeHRaw, ...
                built.workspace.surface.hcrit, ...
                isfield(waterFlux, 'Throughfall') ...
                    && waterFlux.Throughfall > 1e-10, ...
                true, modflowFirstActiveNode(boundaryRequest, NN), ...
                isempty(boundaryRequest));
            probeDeltaH = probeH - hCurrent;
            probeRhsAir = airPressureRow.rhs ...
                - solver.applyTridiagonalBank( ...
                    airPressureRow.h, probeDeltaH);
            probeDeltaP = solveBank(airPressureRow.P, probeRhsAir);
            useDryMassCoordinate = any( ...
                PCurrent + probeDeltaP <= vaporPartialPressure);
        end
        if useDryMassCoordinate
            coordinateName = 'dry_air_mass';
        else
            coordinateName = 'gas_pressure';
            waterRow = waterPressureRow;
            airRow = airPressureRow;
            if heatEnabled
                heatRow = heatPressureRow;
                U_h = U_h_P;
                U_T = U_T_P;
                U_M = zeros(NN,1);
            end
        end
    else
        coordinateName = 'disabled';
    end

    waterFirstActiveNode = modflowFirstActiveNode(boundaryRequest, NN);
    if heatEnabled && airEnabled && useDryMassCoordinate ...
            && waterFirstActiveNode == 1
        if ~any(dryMassSeedActive)
            % A sequential pressure probe can predict a component crossing
            % even when the coupled block direction has a feasible damped
            % pressure step. Try that complete block first; switch charts only
            % when no admissible residual-decreasing pressure step exists.
            [pressureBuilt, pressureResult, pressureRelaxation] = ...
                solvePressureBlockPicardMap( ...
                built, equations, matrices, systems, ...
                waterPressureRow, heatPressureRow, airPressureRow, ...
                hCurrent, TCurrent, PCurrent, U_T_P, dt_s, ...
                waterStorageM, boundaryRequest, heatRequestForMap, ...
                boundaryApplied, heatBoundaryApplied, airBoundaryApplied, ...
                waterFlux, refreshOpts, tolAirRelative, airRelaxation, ...
                residualContext);
            if pressureResult.correctionAccepted
                built = pressureBuilt;
                result = pressureResult;
                airRelaxation = pressureRelaxation;
                return;
            end
        end
        [built, result, airRelaxation] = solveMassBlockPicardMap( ...
            built, equations, matrices, systems, ...
            waterRow, heatRow, airRow, hCurrent, TCurrent, PCurrent, ...
            U_h, U_T, U_M, dryMassCurrent, dryMass_h, dryMass_T, ...
            dryMass_P, dryMassSeedActive, dt_s, waterStorageM, ...
            boundaryRequest, heatRequestForMap, boundaryApplied, ...
            heatBoundaryApplied, airBoundaryApplied, waterFlux, ...
            refreshOpts, tolAirRelative, airRelaxation, residualContext);
        return;
    end
    if blockPicardRequested && heatEnabled && airEnabled ...
            && ~useDryMassCoordinate ...
            && waterFirstActiveNode == 1
        [blockBuilt, blockResult, blockAirRelaxation] = ...
            solvePressureBlockPicardMap( ...
            built, equations, matrices, systems, ...
            waterPressureRow, heatPressureRow, airPressureRow, ...
            hCurrent, TCurrent, PCurrent, U_T_P, dt_s, waterStorageM, ...
            boundaryRequest, heatRequestForMap, boundaryApplied, ...
            heatBoundaryApplied, airBoundaryApplied, waterFlux, ...
            refreshOpts, tolAirRelative, airRelaxation, residualContext);
        if blockResult.correctionAccepted
            built = blockBuilt;
            result = blockResult;
            airRelaxation = blockAirRelaxation;
            return;
        end
    end

    % 1. Water correction.  Other coordinate corrections are zero in this sweep.
    deltaHRaw = solveBank(waterRow.h, waterRow.rhs, ...
        [waterFirstActiveNode, NN]);
    hRaw = hCurrent + deltaHRaw;
    isRaining = isfield(waterFlux, 'Throughfall') ...
        && waterFlux.Throughfall > 1e-10;
    hcrit = built.workspace.surface.hcrit;
    [hNext, ~, ~] = ...
        soil.water.applyHeadUpdateSafeguards( ...
        hCurrent, hRaw, hcrit, isRaining, true, waterFirstActiveNode, ...
        isempty(boundaryRequest));
    deltaH = hNext - hCurrent;

    % 2. Dry-air correction. Normal states use total gas pressure. A predicted
    % component-boundary crossing uses conservative nodal dry-air mass, with
    % pressure recovered only after h/U/M_a are all known.
    CHKAIR = zeros(NN,1);
    airRelativeIncrement = 0;
    airRelativeByNode = zeros(NN,1);
    deltaM = zeros(NN,1);
    dryMassActive = false(NN,1);
    dryMassReaction = zeros(NN,1);
    releasedDryMassNode = false(NN,1);
    deltaP = zeros(NN,1);
    if airEnabled
        rhsAir = airRow.rhs ...
            - solver.applyTridiagonalBank(airRow.h, deltaH);
        if useDryMassCoordinate
            initialDryMassActive = dryMassSeedActive;
            if ~isempty(airRelaxation.activeNodes)
                initialDryMassActive = initialDryMassActive ...
                    | airRelaxation.activeNodes;
            end
            [deltaMRaw, dryMassActive, dryMassReaction] = ...
                solver.solveLowerBoundCorrection( ...
                    airRow.M, rhsAir, dryMassCurrent, ...
                    initialDryMassActive);
            airRelaxation.activeNodes = dryMassActive;
            massDefect = abs(deltaMRaw);
            massDefect(dryMassActive) = 0;
            massScale = max(max(abs(dryMassCurrent)), realmin);
            airRelativeByNode = massDefect ./ ...
                max(abs(dryMassCurrent), 1e-6 * massScale);
            airRelativeIncrement = max(airRelativeByNode);
            rawAirPass = airRelativeIncrement < tolAirRelative;
            releasedDryMassNode = initialDryMassActive & ~dryMassActive;
            if any(dryMassActive)
                % A bound-active coordinate must land exactly on M_a=0.
                % Relaxing -M_a would leave a positive linear target while
                % the constitutive projection below sets rho_a=0, creating
                % an inconsistent mass/pressure pair.
                deltaM = deltaMRaw;
                airRelaxation.coordinate = coordinateName;
                airRelaxation.previousCorrection = [];
                airRelaxation.omega = 1;
            else
                if any(releasedDryMassNode)
                    % A primal-dual release has the correct sign but an
                    % unbounded relative step from M_a=0. Start that one
                    % nonsmooth transition at the existing Aitken lower
                    % bound; subsequent free iterations adapt normally.
                    airRelaxation.coordinate = coordinateName;
                    airRelaxation.previousCorrection = [];
                    airRelaxation.omega = 0.05;
                    airRelaxation.activeNodes = dryMassActive;
                end
                [massNext, airRelaxation] = relaxDryAirCorrection( ...
                    dryMassCurrent, deltaMRaw, airRelaxation, rawAirPass, ...
                    coordinateName);
                deltaM = massNext - dryMassCurrent;
            end
        else
            deltaPRaw = solveBank(airRow.P, rhsAir);
            CHKAIR = abs(deltaPRaw);
            airRelativeByNode = CHKAIR ./ max(abs(PCurrent), 1);
            airRelativeIncrement = max(airRelativeByNode);
            rawAirPass = airRelativeIncrement < tolAirRelative;
            [pressureNext, airRelaxation] = relaxDryAirCorrection( ...
                PCurrent, deltaPRaw, airRelaxation, rawAirPass, ...
                coordinateName);
            deltaP = pressureNext - PCurrent;
        end
    end

    % 3. Internal-energy correction at the forward-sweep h/P_g corrections.
    CHKHEAT = zeros(NN,1);
    temperatureIncrement = zeros(NN,1);
    heatFlux = struct('QET', 0, 'QEB', 0);
    airFlux = struct();
    airFluxContext = struct();
    waterFluxContext = struct();
    inversion = struct('allConverged', true, 'convergedMask', true(NN,1));
    correctionScale = 1;
    mapEnergyResidual = NaN;
    mapEnergyStorageBeforeMismatch = NaN;
    if heatEnabled
        rhsHeat = heatRow.rhs ...
            - solver.applyTridiagonalBank(heatRow.h, deltaH);
        if airEnabled && useDryMassCoordinate
            rhsHeat = rhsHeat ...
                - solver.applyTridiagonalBank(heatRow.M, deltaM);
        else
            rhsHeat = rhsHeat ...
                - solver.applyTridiagonalBank(heatRow.P, deltaP);
        end
        deltaU = solveBank(heatRow.U, rhsHeat);
        % Use the unscaled Picard-map defect as the conservative convergence
        % gate. Line-search damping changes the applied correction, not this
        % fixed-point defect; temperatureIncrement records the applied T change.
        CHKHEAT = abs(deltaU) ./ U_T;
    end

    % 4. Back-substitute the complete gas/U corrections into the same frozen
    % water row. This returns a TOTAL head correction, not an increment added
    % to the forward-sweep head. Reapply the existing Richards safeguards to
    % the back-substituted proposal before installing any trial coordinate.
    if useDryMassCoordinate
        deltaGasCoordinate = deltaM;
    else
        deltaGasCoordinate = deltaP;
    end
    [deltaHBackRaw, waterCoupledRhs, waterBackResidualRaw] = ...
        solver.backSubstituteWaterCorrection( ...
            waterRow, deltaU, deltaGasCoordinate, ...
            [waterFirstActiveNode, NN]);
    hBackRaw = hCurrent + deltaHBackRaw;
    [hNext, CHKWATER, dryHeadProjected] = ...
        soil.water.applyHeadUpdateSafeguards( ...
            hCurrent, hBackRaw, hcrit, isRaining, true, ...
            waterFirstActiveNode, isempty(boundaryRequest));
    deltaH = hNext - hCurrent;
    if useDryMassCoordinate && heatEnabled
        % The dry-air disappearance limit is strongly coupled: sequential
        % scalar back substitution is not a Newton step and can cycle between
        % h, U, and M_a active sets. Solve the same frozen three rows once as
        % their native block-tridiagonal Jacobian.
        [deltaHBlockRaw, deltaUBlockRaw, deltaMBlockRaw, dryMassActive, ...
                dryMassReaction] = ...
            solver.solveCoupledDryAirMassCorrection( ...
                waterRow, heatRow, airRow, waterRow.rhs, ...
                heatRow.rhs, airRow.rhs, dryMassCurrent, ...
                initialDryMassActive);
        [hBlock, ~, ~] = soil.water.applyHeadUpdateSafeguards( ...
            hCurrent, hCurrent + deltaHBlockRaw, hcrit, isRaining, ...
            true, waterFirstActiveNode, isempty(boundaryRequest));
        % A fixed fraction keeps the semismooth bound transition inside the
        % local basin while preserving the exact M_a=0 value on active nodes.
        % The raw block defect remains the convergence gate, so this changes
        % globalization only and cannot manufacture convergence.
        blockScale = 0.1;
        deltaHBlockSafeguarded = hBlock - hCurrent;
        deltaH = blockScale * deltaHBlockSafeguarded;
        deltaU = blockScale * deltaUBlockRaw;
        deltaM = deltaMBlockRaw;
        deltaM(~dryMassActive) = ...
            blockScale * deltaMBlockRaw(~dryMassActive);
        releasedDryMassNode = initialDryMassActive & ~dryMassActive;
        airRelaxation.activeNodes = dryMassActive;
        airRelaxation.previousCorrection = [];
        airRelaxation.omega = 1;
        % Convergence is judged with the undamped block correction. Using
        % the applied 0.1 step here would reduce every reported defect by a
        % factor of ten and could manufacture a false nonlinear PASS.
        CHKHEAT = abs(deltaUBlockRaw) ./ U_T;
        massDefect = abs(deltaMBlockRaw);
        massDefect(dryMassActive) = 0;
        massScale = max(max(abs(dryMassCurrent)), realmin);
        airRelativeByNode = massDefect ./ ...
            max(abs(dryMassCurrent), 1e-6 * massScale);
        airRelativeIncrement = max(airRelativeByNode);
        hNext = hCurrent + deltaHBlockSafeguarded;
        [~, CHKWATER, dryHeadProjected] = ...
            soil.water.applyHeadUpdateSafeguards( ...
                hCurrent, hNext, hcrit, isRaining, true, ...
                waterFirstActiveNode, isempty(boundaryRequest));
    end
    if airEnabled && useDryMassCoordinate
        if heatEnabled
            deltaTCoordinate = (deltaU - U_h .* deltaH ...
                - U_M .* deltaM) ./ U_T;
        else
            deltaTCoordinate = zeros(NN,1);
        end
        deltaP = (deltaM - dryMass_h .* deltaH ...
            - dryMass_T .* deltaTCoordinate) ./ dryMass_P;
        CHKAIR = abs(deltaP);
    end
    surfaceDryHeadProjected = hasAtmosphericSurface(boundaryRequest) ...
        && dryHeadProjected;
    waterBackResidualSafeguarded = waterCoupledRhs ...
        - solver.applyTridiagonalBank(waterRow.h, deltaH);
    activeWater = waterFirstActiveNode:NN;
    maxWaterBackResidualRaw = max(abs( ...
        waterBackResidualRaw(activeWater)));
    maxWaterBackResidualSafeguarded = max(abs( ...
        waterBackResidualSafeguarded(activeWater)));

    % 5. Install the correction through one coupled physical line search.
    % U inversion and dry-air feasibility both scale h/U/P_g together; no
    % state component is clipped independently.
    UCurrent = built.workspace.step.trial.U(:);
    airAdmissibility = disabledDryAirAdmissibility();
    correctionAccepted = true;
    if heatEnabled || airEnabled
        referenceAirDensity = [];
        referenceAirStorage = [];
        if airEnabled
            referenceAirDensity = built.workspace.derived.RHODA(:);
            referenceAirStorage = equations.dry_air.dryAirStorage;
        end
        [built, inversion, correctionScale, airAdmissibility, ...
            correctionAccepted] = applyCoupledCorrectionWithLineSearch( ...
                built, hCurrent, PCurrent, UCurrent, ...
                deltaH, deltaP, deltaU, refreshOpts, heatEnabled, ...
                airEnabled, referenceAirDensity, referenceAirStorage, ...
                dryMassCurrent, deltaM, dryMassActive);
        if heatEnabled && correctionAccepted
            TNext = built.workspace.step.trial.T(:);
            temperatureIncrement = abs(TNext - TCurrent);
            deltaHApplied = built.workspace.step.trial.h(:) - hCurrent;
            deltaUApplied = built.workspace.step.trial.U(:) - UCurrent;
            deltaPApplied = built.workspace.step.trial.p_gas(:) - PCurrent;
            [QET, QEB] = soil.thermal.calculateCoupledBoundaryFlux( ...
                systemHeat.SAVE, matrices.heat, TNext, TNext-TCurrent, ...
                deltaHApplied, deltaUApplied, deltaPApplied, ...
                built.domain.dz, dt_s);
            heatFlux = struct('QET', QET, 'QEB', QEB);
            UBeforeMap = soil.thermal.assembleNodalInternalEnergy( ...
                energyStorageN, zeros(NN-1,2), zeros(NN-1,2), ...
                zeros(NN-1,2), built.domain.dz);
            mapBalance = soil.thermal.computeColumnEnergyBalance( ...
                UBeforeMap, built.workspace.step.trial.U(:), ...
                equations.heat.rootEnergySink, built.domain.dz, ...
                heatFlux, dt_s);
            mapEnergyResidual = mapBalance.residual_J_cm2;
            mapEnergyStorageBeforeMismatch = sum( ...
                mapBalance.storage_before_J_cm2 ...
                - sum(controlVolumeWeights(built.domain.dz) ...
                    .* built.Y.U(:)));
        end
    else
        built.workspace.step.trial.h = reshape( ...
            hCurrent + deltaH, size(built.workspace.step.trial.h));
    end

    % The h/U/P_g correction map combines the assembled equation rows after
    % the water-only Thomas solve. For a MODFLOW fixed-head request, nodes at
    % and below the water table are constraints, not correction coordinates.
    % Reapply the fixed interface and its hydrostatic continuation after the
    % final coupled map (including heat line-search damping) before any flux,
    % storage, or accepted-state calculation reads trial.h.
    built = enforceModflowHydrostaticConstraint(built, boundaryRequest);
    hAccepted = built.workspace.step.trial.h(:);
    % Keep the h-only value as a finite predictor for direct kinematic
    % active-set exits. Once a candidate reaches the nonlinear gate,
    % The exact ledger replaces it with the h/T/P_g/storage reaction.
    hBoundary = soil.water.calculatesSoilWaterFluxes( ...
        systemWater.SAVE, struct('hh', hAccepted), waterFirstActiveNode);
    waterFlux.QMT = hBoundary.QMT;
    waterFlux.QMB = hBoundary.QMB;
    if correctionAccepted
        waterFluxContext = struct( ...
            'SAVE', systemWater.SAVE, ...
            'matrix', matrices.water, ...
            'stateBefore', struct('h', hCurrent, ...
                'T', TCurrent, 'P', PCurrent), ...
            'stateAfter', struct( ...
                'h', built.workspace.step.trial.h(:), ...
                'T', built.workspace.step.trial.T(:), ...
                'P', built.workspace.step.trial.p_gas(:)), ...
            'storageBefore', waterStorageM, ...
            'bottom_node', waterFirstActiveNode);
    end
    if airEnabled
        % No P_g-only flux is published here. The exact dry-air reaction is
        %
        %   q_a = q_a(h,T,P_g,B_a),
        %
        % where B_a is the endpoint dry-air storage. It is recovered once,
        % after the accepted state has been refreshed, by
        % calculateCoupledBoundaryFlux. Keeping airFlux empty prevents a
        % linearized pressure-only value from entering the public contract.
        if correctionAccepted
            airFluxContext = struct( ...
                'SAVE', systemAir.SAVE, ...
                'matrix', matrices.dry_air, ...
                'stateBefore', struct('h', hCurrent, ...
                    'T', TCurrent, 'P', PCurrent), ...
                'stateAfter', struct( ...
                    'h', built.workspace.step.trial.h(:), ...
                    'T', built.workspace.step.trial.T(:), ...
                    'P', built.workspace.step.trial.p_gas(:)), ...
                'storageBefore', equations.dry_air.dryAirStorage, ...
                'active_reaction', dryMassReaction, ...
                'active_nodes', dryMassActive, ...
                'boundary_applied', airBoundaryApplied);
        end
    end

    failure = '';
    if ~inversion.allConverged
        failed = find(~inversion.convergedMask, 1, 'first');
        if ~isempty(failed) && isfield(inversion, 'failureReason')
            failure = inversion.failureReason{failed};
        end
    end

    result = struct('CHKWATER', CHKWATER, 'CHKAIR', CHKAIR, ...
        'CHKHEAT', CHKHEAT, 'temperatureIncrement', temperatureIncrement, ...
        'waterFlux', waterFlux, 'heatFlux', heatFlux, 'airFlux', airFlux, ...
        'waterFluxContext', waterFluxContext, ...
        'airFluxContext', airFluxContext, ...
        'boundaryApplied', boundaryApplied, ...
        'heatBoundaryApplied', heatBoundaryApplied, ...
        'heatBoundaryRequest', heatRequestForMap, ...
        'airBoundaryApplied', airBoundaryApplied, ...
        'airRelativeIncrement', airRelativeIncrement, ...
        'airMassCoordinateIncrement', max(abs(deltaM)), ...
        'airCoordinate', coordinateName, ...
        'surfaceDryHeadProjected', surfaceDryHeadProjected, ...
        'temperatureInversion', inversion, ...
        'airAdmissibility', airAdmissibility, ...
        'correctionAccepted', correctionAccepted, ...
        'diagnostics', struct( ...
            'max_delta_h_cm', max(abs(deltaH)), ...
            'max_delta_p', max(abs(deltaP)), ...
            'max_delta_U', max(abs(deltaU)), ...
            'max_delta_dry_air_mass_g_cm2', max(abs(deltaM)), ...
            'dry_air_active_node_count', nnz(dryMassActive), ...
            'max_dry_air_active_reaction', ...
                max([0; abs(dryMassReaction(dryMassActive))]), ...
            'minimum_dry_air_active_reaction', ...
                min([0; dryMassReaction(dryMassActive)]), ...
            'dry_air_release_transition', any(releasedDryMassNode), ...
            'air_relative_by_node', airRelativeByNode, ...
            'max_air_relative_node', find( ...
                airRelativeByNode == max(airRelativeByNode), 1, 'first'), ...
            'correction_scale', correctionScale, ...
            'map_energy_residual_J_cm2', mapEnergyResidual, ...
            'map_energy_storage_before_mismatch_J_cm2', ...
                mapEnergyStorageBeforeMismatch, ...
            'water_backsolve_raw_linear_residual', ...
                maxWaterBackResidualRaw, ...
            'water_backsolve_safeguarded_linear_residual', ...
                maxWaterBackResidualSafeguarded, ...
            'temperature_failure', failure, ...
            'dry_air_failure', airAdmissibility.failure_reason));
end

% -----------------------------------------------------------------------------
function context = prepareExactResidualContext( ...
        built, dt_s, kt, waterBoundaryFlux, waterRootSink, ...
        soilBoundaryRequest, waterRequest, heatRequest, airRequest, ...
        refreshOpts)
% Freeze the accepted anchors used by every exact/globalized evaluation.
%
% The context is constructed once before the Picard loop. Active-set water
% requests are replaced inside each map, while accepted storage anchors and
% forcing remain immutable. This adapter never commits state.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    context = soil.coupling.prepareResidualContext(built, dt_s, kt);
    context.water_boundary_flux = waterBoundaryFlux;
    context.water_root_sink = waterRootSink;
    context.refresh_opts = refreshOpts;
    if ~isempty(soilBoundaryRequest)
        context.water_boundary_request = soilBoundaryRequest.water;
        context.heat_boundary_request = soilBoundaryRequest.heat;
        context.air_boundary_request = soilBoundaryRequest.dry_air;
    else
        context.water_boundary_request = waterRequest;
        context.heat_boundary_request = heatRequest;
        context.air_boundary_request = airRequest;
    end
end

% -----------------------------------------------------------------------------
function exact = reuseExactCandidateForState( ...
        built, mapResult, context, waterRequest, heatRequest)
% Reuse a block line-search ledger only when it belongs to this exact state.
%
% A block line search already evaluates the exact residual at each distinct
% alpha candidate.
% Re-evaluating its accepted alpha would duplicate a complete constitutive
% refresh. Ordinary Picard has no such record, so the acceptance gate creates
% one here. The trial-state hash prevents stale ledger reuse after any active-
% set restart or subsequent map update.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    exact = struct();
    if isfield(mapResult, 'exactResidual') ...
            && isstruct(mapResult.exactResidual) ...
            && ~isempty(fieldnames(mapResult.exactResidual)) ...
            && isfield(mapResult.exactResidual, 'state_hash') ...
            && strcmp(mapResult.exactResidual.state_hash, ...
                soil.coupling.hashTrialState(built))
        exact = mapResult.exactResidual;
        return;
    end
    context.water_boundary_request = waterRequest;
    context.heat_boundary_request = heatRequest;
    exact = soil.coupling.evaluateResidualLedger(built, context);
end

% -----------------------------------------------------------------------------
function [built, result, airRelaxation] = solvePressureBlockPicardMap( ...
        built, equations, matrices, systems, ...
        waterRow, heatRow, airRow, hCurrent, TCurrent, PCurrent, U_T, ...
        dt_s, waterStorageM, boundaryRequest, heatRequest, ...
        boundaryApplied, heatBoundaryApplied, airBoundaryApplied, ...
        waterFlux, refreshOpts, tolAirRelative, airRelaxation, context)
% Solve and globalize one ordinary-pressure block Modified-Picard map.
%
% All three rows were assembled from the same frozen coefficients. The
% block Thomas solve changes no coefficient and is therefore Picard, not
% Newton. One alpha scales h/U/P_g together. Every alpha candidate is detached,
% temperature is recovered from U, and the exact residual owner decides whether
% the nonlinear merit decreases. No candidate is committed here.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    NN = built.domain.NN;
    linearization = struct('available', true, ...
        'coordinate_order', {{'h','U','P_g'}}, ...
        'equation_order', {{'water','energy','dry_air'}}, ...
        'rows', struct('water', waterRow, 'energy', heatRow, ...
            'dry_air', airRow));
    [delta, linearInfo] = solver.solveLaggedPicardBlock(linearization);

    waterFirstActiveNode = modflowFirstActiveNode(boundaryRequest, NN);
    isRaining = isfield(waterFlux, 'Throughfall') ...
        && waterFlux.Throughfall > 1e-10;
    [hTarget, ~, dryHeadProjected] = ...
        soil.water.applyHeadUpdateSafeguards( ...
            hCurrent, hCurrent + delta.h, ...
            built.workspace.surface.hcrit, isRaining, true, ...
            waterFirstActiveNode, isempty(boundaryRequest));
    directionH = hTarget - hCurrent;
    UCurrent = built.workspace.step.trial.U(:);

    context.water_boundary_request = boundaryRequest;
    context.heat_boundary_request = heatRequest;
    phiBefore = blockRowMerit(waterRow, heatRow, airRow, context.scales);
    alphaAccepted = 0;
    exactAfter = [];
    inversion = struct('allConverged', true, ...
        'convergedMask', true(NN,1), ...
        'failureReason', {repmat({''}, NN, 1)});
    airAdmissibility = disabledDryAirAdmissibility();
    candidateAccepted = false;
    for alpha = 0.5 .^ (0:12)
        candidate = built;
        candidate.workspace.step.trial.h = reshape( ...
            hCurrent + alpha * directionH, ...
            size(candidate.workspace.step.trial.h));
        candidate.workspace.step.trial.U = reshape( ...
            UCurrent + alpha * delta.U, ...
            size(candidate.workspace.step.trial.U));
        candidate.workspace.step.trial.p_gas = reshape( ...
            PCurrent + alpha * delta.P, ...
            size(candidate.workspace.step.trial.p_gas));
        [candidate.workspace.step.trial, trialInversion] = recoverTrialTemperature( ...
            candidate, refreshOpts);
        inversion = trialInversion;
        if ~trialInversion.allConverged
            continue;
        end
        [candidateCheck, candidateEquations] = ...
            soil.refreshSoilCoefficients(candidate, refreshOpts);
        trialAir = soil.gas.evaluateDryAirAdmissibility( ...
            candidateCheck.workspace.derived.RHODA, ...
            candidateEquations.dry_air.dryAirStorage);
        airAdmissibility = trialAir;
        if ~trialAir.admissible
            continue;
        end
        trialExact = soil.coupling.evaluateResidualLedger( ...
            candidate, context);
        armijoSlack = 1e-12 * max(1, phiBefore);
        sufficientDecrease = trialExact.merit.phi ...
            <= (1 - 1e-4 * alpha) * phiBefore + armijoSlack;
        if trialExact.feasibility.admissible && sufficientDecrease
            built = candidate;
            exactAfter = trialExact;
            alphaAccepted = alpha;
            candidateAccepted = true;
            break;
        end
    end

    rawAirRelativeByNode = abs(delta.P) ./ max(abs(PCurrent), 1);
    rawAirRelativeIncrement = max(rawAirRelativeByNode);
    CHKWATER = abs(delta.h);
    CHKAIR = abs(delta.P);
    CHKHEAT = abs(delta.U) ./ U_T;
    temperatureIncrement = zeros(NN,1);
    heatFlux = struct('QET', 0, 'QEB', 0);
    airFlux = struct();
    waterFluxContext = struct();
    airFluxContext = struct();
    mapEnergyResidual = NaN;
    mapEnergyStorageBeforeMismatch = NaN;
    phiAfter = NaN;
    if candidateAccepted
        trial = built.workspace.step.trial;
        temperatureIncrement = abs(trial.T(:) - TCurrent);
        deltaHApplied = trial.h(:) - hCurrent;
        deltaUApplied = trial.U(:) - UCurrent;
        deltaPApplied = trial.p_gas(:) - PCurrent;
        [QET, QEB] = soil.thermal.calculateCoupledBoundaryFlux( ...
            systems.heat.SAVE, matrices.heat, trial.T(:), ...
            trial.T(:) - TCurrent, deltaHApplied, deltaUApplied, ...
            deltaPApplied, built.domain.dz, dt_s);
        heatFlux = struct('QET', QET, 'QEB', QEB);
        waterFlux = exactAfter.ledgers.water_boundary_flux;
        mapEnergyResidual = exactAfter.ledgers.energy.residual_J_cm2;
        mapEnergyStorageBeforeMismatch = ...
            exactAfter.ledgers.energy_coordinate_mismatch_J_cm3;
        phiAfter = exactAfter.merit.phi;
        waterFluxContext = struct( ...
            'SAVE', systems.water.SAVE, ...
            'matrix', matrices.water, ...
            'stateBefore', struct('h', hCurrent, ...
                'T', TCurrent, 'P', PCurrent), ...
            'stateAfter', struct('h', trial.h(:), ...
                'T', trial.T(:), 'P', trial.p_gas(:)), ...
            'storageBefore', waterStorageM, ...
            'bottom_node', waterFirstActiveNode);
        airFluxContext = struct( ...
            'SAVE', systems.air.SAVE, ...
            'matrix', matrices.dry_air, ...
            'stateBefore', struct('h', hCurrent, ...
                'T', TCurrent, 'P', PCurrent), ...
            'stateAfter', struct('h', trial.h(:), ...
                'T', trial.T(:), 'P', trial.p_gas(:)), ...
            'storageBefore', equations.dry_air.dryAirStorage, ...
            'active_reaction', zeros(NN,1), ...
            'active_nodes', false(NN,1), ...
            'boundary_applied', airBoundaryApplied);
        airRelaxation.coordinate = 'gas_pressure';
        airRelaxation.previousCorrection = delta.P;
        airRelaxation.omega = alphaAccepted;
        airRelaxation.minimumOmega = min( ...
            airRelaxation.minimumOmega, alphaAccepted);
        airRelaxation.applied = airRelaxation.applied ...
            + double(alphaAccepted < 1);
    end

    maxAirNode = find(rawAirRelativeByNode ...
        == max(rawAirRelativeByNode), 1, 'first');
    result = struct('CHKWATER', CHKWATER, 'CHKAIR', CHKAIR, ...
        'CHKHEAT', CHKHEAT, ...
        'temperatureIncrement', temperatureIncrement, ...
        'waterFlux', waterFlux, 'heatFlux', heatFlux, 'airFlux', airFlux, ...
        'waterFluxContext', waterFluxContext, ...
        'airFluxContext', airFluxContext, ...
        'boundaryApplied', boundaryApplied, ...
        'heatBoundaryApplied', heatBoundaryApplied, ...
        'heatBoundaryRequest', heatRequest, ...
        'airBoundaryApplied', airBoundaryApplied, ...
        'airRelativeIncrement', rawAirRelativeIncrement, ...
        'airMassCoordinateIncrement', 0, ...
        'airCoordinate', 'gas_pressure_block_picard', ...
        'surfaceDryHeadProjected', ...
            hasAtmosphericSurface(boundaryRequest) && dryHeadProjected, ...
        'temperatureInversion', inversion, ...
        'airAdmissibility', airAdmissibility, ...
        'correctionAccepted', candidateAccepted, ...
        'exactResidual', exactAfter, ...
        'diagnostics', struct( ...
            'max_delta_h_cm', max(abs(delta.h)), ...
            'max_delta_p', max(abs(delta.P)), ...
            'max_delta_U', max(abs(delta.U)), ...
            'max_delta_dry_air_mass_g_cm2', 0, ...
            'dry_air_active_node_count', 0, ...
            'max_dry_air_active_reaction', 0, ...
            'minimum_dry_air_active_reaction', 0, ...
            'dry_air_release_transition', false, ...
            'air_relative_by_node', rawAirRelativeByNode, ...
            'max_air_relative_node', maxAirNode, ...
            'correction_scale', alphaAccepted, ...
            'map_energy_residual_J_cm2', mapEnergyResidual, ...
            'map_energy_storage_before_mismatch_J_cm2', ...
                mapEnergyStorageBeforeMismatch, ...
            'water_backsolve_raw_linear_residual', ...
                linearInfo.max_abs_linear_residual, ...
            'water_backsolve_safeguarded_linear_residual', ...
                linearInfo.max_abs_linear_residual, ...
            'block_phi_before', phiBefore, ...
            'block_phi_after', phiAfter, ...
            'block_relative_linear_residual', ...
                linearInfo.relative_linear_residual, ...
            'block_minimum_rcond', ...
                linearInfo.minimum_block_rcond, ...
            'block_newton_used', false, ...
            'temperature_failure', inversionFailure(inversion), ...
            'dry_air_failure', airAdmissibility.failure_reason));
end

% -----------------------------------------------------------------------------
function [built, result, airRelaxation] = solveMassBlockPicardMap( ...
        built, equations, matrices, systems, waterRow, heatRow, airRow, ...
        hCurrent, TCurrent, PCurrent, U_h, U_T, U_M, dryMassCurrent, ...
        dryMass_h, dryMass_T, dryMass_P, initialActive, dt_s, ...
        waterStorageM, boundaryRequest, heatRequest, boundaryApplied, ...
        heatBoundaryApplied, airBoundaryApplied, waterFlux, refreshOpts, ...
        tolAirRelative, airRelaxation, context)
% Solve and globalize one frozen h/U/M_a block Modified-Picard map.
%
% This degeneration chart is used only when pressure ceases to be a
% safe component coordinate. Coefficients remain frozen and the active-set
% solve enforces M_a >= 0 without invoking the production Newton-named block.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    NN = built.domain.NN;
    rows = struct('water', waterRow, 'energy', heatRow, ...
        'dry_air', airRow);
    [delta, linearInfo] = solver.solveLaggedPicardMassBlock( ...
        rows, dryMassCurrent, initialActive);
    waterFirstActiveNode = modflowFirstActiveNode(boundaryRequest, NN);
    isRaining = isfield(waterFlux, 'Throughfall') ...
        && waterFlux.Throughfall > 1e-10;
    [hTarget, ~, dryHeadProjected] = ...
        soil.water.applyHeadUpdateSafeguards( ...
            hCurrent, hCurrent + delta.h, ...
            built.workspace.surface.hcrit, isRaining, true, ...
            waterFirstActiveNode, isempty(boundaryRequest));
    directionH = hTarget - hCurrent;
    directionT = (delta.U - U_h .* directionH ...
        - U_M .* delta.M) ./ U_T;
    directionP = (delta.M - dryMass_h .* directionH ...
        - dryMass_T .* directionT) ./ dryMass_P;
    UCurrent = built.workspace.step.trial.U(:);
    context.water_boundary_request = boundaryRequest;
    context.heat_boundary_request = heatRequest;
    phiBefore = massBlockRowMerit( ...
        waterRow, heatRow, airRow, context.scales);

    alphaAccepted = 0;
    exactAfter = [];
    inversion = struct('allConverged', true, ...
        'convergedMask', true(NN,1), ...
        'failureReason', {repmat({''}, NN, 1)});
    airAdmissibility = disabledDryAirAdmissibility();
    candidateAccepted = false;
    for alpha = 0.5 .^ (0:12)
        candidate = built;
        candidate.workspace.step.trial.h = reshape( ...
            hCurrent + alpha * directionH, ...
            size(candidate.workspace.step.trial.h));
        candidate.workspace.step.trial.U = reshape( ...
            UCurrent + alpha * delta.U, ...
            size(candidate.workspace.step.trial.U));
        candidate.workspace.step.trial.p_gas = reshape( ...
            PCurrent + alpha * directionP, ...
            size(candidate.workspace.step.trial.p_gas));
        activeAtBoundary = alpha == 1 && any(delta.active);
        if activeAtBoundary
            candidate = projectDryAirComponentBoundary( ...
                candidate, delta.active, refreshOpts, false);
        end
        [candidate.workspace.step.trial, trialInversion] = recoverTrialTemperature( ...
            candidate, refreshOpts);
        inversion = trialInversion;
        if ~trialInversion.allConverged
            continue;
        end
        if activeAtBoundary
            candidate = projectDryAirComponentBoundary( ...
                candidate, delta.active, refreshOpts, true);
            [candidate.workspace.step.trial, trialInversion] = recoverTrialTemperature( ...
                candidate, refreshOpts);
            inversion = trialInversion;
            if ~trialInversion.allConverged
                continue;
            end
        end
        [candidateCheck, candidateEquations] = ...
            soil.refreshSoilCoefficients(candidate, refreshOpts);
        trialAir = soil.gas.evaluateDryAirAdmissibility( ...
            candidateCheck.workspace.derived.RHODA, ...
            candidateEquations.dry_air.dryAirStorage);
        airAdmissibility = trialAir;
        if ~trialAir.admissible
            continue;
        end
        trialExact = soil.coupling.evaluateResidualLedger(candidate, context);
        armijoSlack = 1e-12 * max(1, phiBefore);
        sufficientDecrease = trialExact.merit.phi ...
            <= (1 - 1e-4 * alpha) * phiBefore + armijoSlack;
        if trialExact.feasibility.admissible && sufficientDecrease
            built = candidate;
            exactAfter = trialExact;
            alphaAccepted = alpha;
            candidateAccepted = true;
            break;
        end
    end

    massDefect = abs(delta.M);
    massDefect(delta.active) = 0;
    massScale = max(max(abs(dryMassCurrent)), realmin);
    airRelativeByNode = massDefect ...
        ./ max(abs(dryMassCurrent), 1e-6 * massScale);
    airRelativeIncrement = max(airRelativeByNode);
    CHKWATER = abs(delta.h);
    CHKAIR = abs(directionP);
    CHKHEAT = abs(delta.U) ./ U_T;
    temperatureIncrement = zeros(NN,1);
    heatFlux = struct('QET', 0, 'QEB', 0);
    airFlux = struct();
    waterFluxContext = struct();
    airFluxContext = struct();
    mapEnergyResidual = NaN;
    mapEnergyStorageBeforeMismatch = NaN;
    phiAfter = NaN;
    if candidateAccepted
        trial = built.workspace.step.trial;
        temperatureIncrement = abs(trial.T(:) - TCurrent);
        deltaHApplied = trial.h(:) - hCurrent;
        deltaUApplied = trial.U(:) - UCurrent;
        deltaPApplied = trial.p_gas(:) - PCurrent;
        [QET, QEB] = soil.thermal.calculateCoupledBoundaryFlux( ...
            systems.heat.SAVE, matrices.heat, trial.T(:), ...
            trial.T(:) - TCurrent, deltaHApplied, deltaUApplied, ...
            deltaPApplied, built.domain.dz, dt_s);
        heatFlux = struct('QET', QET, 'QEB', QEB);
        waterFlux = exactAfter.ledgers.water_boundary_flux;
        mapEnergyResidual = exactAfter.ledgers.energy.residual_J_cm2;
        mapEnergyStorageBeforeMismatch = ...
            exactAfter.ledgers.energy_coordinate_mismatch_J_cm3;
        phiAfter = exactAfter.merit.phi;
        waterFluxContext = struct( ...
            'SAVE', systems.water.SAVE, 'matrix', matrices.water, ...
            'stateBefore', struct('h', hCurrent, ...
                'T', TCurrent, 'P', PCurrent), ...
            'stateAfter', struct('h', trial.h(:), ...
                'T', trial.T(:), 'P', trial.p_gas(:)), ...
            'storageBefore', waterStorageM, ...
            'bottom_node', waterFirstActiveNode);
        airFluxContext = struct( ...
            'SAVE', systems.air.SAVE, 'matrix', matrices.dry_air, ...
            'stateBefore', struct('h', hCurrent, ...
                'T', TCurrent, 'P', PCurrent), ...
            'stateAfter', struct('h', trial.h(:), ...
                'T', trial.T(:), 'P', trial.p_gas(:)), ...
            'storageBefore', equations.dry_air.dryAirStorage, ...
            'active_reaction', delta.reaction, ...
            'active_nodes', delta.active, ...
            'boundary_applied', airBoundaryApplied);
        airRelaxation.coordinate = 'dry_air_mass';
        airRelaxation.previousCorrection = delta.M;
        airRelaxation.activeNodes = delta.active;
        airRelaxation.omega = alphaAccepted;
        airRelaxation.minimumOmega = min( ...
            airRelaxation.minimumOmega, alphaAccepted);
        airRelaxation.applied = airRelaxation.applied ...
            + double(alphaAccepted < 1);
    end
    maxAirNode = find(airRelativeByNode ...
        == max(airRelativeByNode), 1, 'first');
    result = struct('CHKWATER', CHKWATER, 'CHKAIR', CHKAIR, ...
        'CHKHEAT', CHKHEAT, ...
        'temperatureIncrement', temperatureIncrement, ...
        'waterFlux', waterFlux, 'heatFlux', heatFlux, 'airFlux', airFlux, ...
        'waterFluxContext', waterFluxContext, ...
        'airFluxContext', airFluxContext, ...
        'boundaryApplied', boundaryApplied, ...
        'heatBoundaryApplied', heatBoundaryApplied, ...
        'heatBoundaryRequest', heatRequest, ...
        'airBoundaryApplied', airBoundaryApplied, ...
        'airRelativeIncrement', airRelativeIncrement, ...
        'airMassCoordinateIncrement', max(abs(delta.M)), ...
        'airCoordinate', 'dry_air_mass_block_picard', ...
        'surfaceDryHeadProjected', ...
            hasAtmosphericSurface(boundaryRequest) && dryHeadProjected, ...
        'temperatureInversion', inversion, ...
        'airAdmissibility', airAdmissibility, ...
        'correctionAccepted', candidateAccepted, ...
        'exactResidual', exactAfter, ...
        'diagnostics', struct( ...
            'max_delta_h_cm', max(abs(delta.h)), ...
            'max_delta_p', max(abs(directionP)), ...
            'max_delta_U', max(abs(delta.U)), ...
            'max_delta_dry_air_mass_g_cm2', max(abs(delta.M)), ...
            'dry_air_active_node_count', nnz(delta.active), ...
            'max_dry_air_active_reaction', ...
                max([0; abs(delta.reaction(delta.active))]), ...
            'minimum_dry_air_active_reaction', ...
                min([0; delta.reaction(delta.active)]), ...
            'dry_air_release_transition', ...
                any(initialActive & ~delta.active), ...
            'air_relative_by_node', airRelativeByNode, ...
            'max_air_relative_node', maxAirNode, ...
            'correction_scale', alphaAccepted, ...
            'map_energy_residual_J_cm2', mapEnergyResidual, ...
            'map_energy_storage_before_mismatch_J_cm2', ...
                mapEnergyStorageBeforeMismatch, ...
            'water_backsolve_raw_linear_residual', ...
                linearInfo.max_abs_linear_residual, ...
            'water_backsolve_safeguarded_linear_residual', ...
                linearInfo.max_abs_linear_residual, ...
            'block_phi_before', phiBefore, ...
            'block_phi_after', phiAfter, ...
            'block_relative_linear_residual', ...
                linearInfo.relative_linear_residual, ...
            'block_minimum_rcond', ...
                linearInfo.minimum_block_rcond, ...
            'block_newton_used', false, ...
            'temperature_failure', inversionFailure(inversion), ...
            'dry_air_failure', airAdmissibility.failure_reason));
end

% -----------------------------------------------------------------------------
function phi = massBlockRowMerit(waterRow, heatRow, airRow, scales)
% Evaluate the exact scaled norm; coordinate transforms leave each rhs intact.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    rhoWater = norm(waterRow.rhs(:) ./ scales.water_cm_s, 2) ...
        / sqrt(numel(waterRow.rhs));
    rhoEnergy = norm(heatRow.rhs(:) ./ scales.energy_J_cm2_s, 2) ...
        / sqrt(numel(heatRow.rhs));
    rhoAir = norm(airRow.rhs(:) ./ scales.dry_air_g_cm2_s, 2) ...
        / sqrt(numel(airRow.rhs));
    phi = 0.5 * (rhoWater^2 + rhoEnergy^2 + rhoAir^2);
end

% -----------------------------------------------------------------------------
function phi = blockRowMerit(waterRow, heatRow, airRow, scales)
% Evaluate the exact scaled residual norm from the assembled rows.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    rhoWater = norm(waterRow.rhs(:) ./ scales.water_cm_s, 2) ...
        / sqrt(numel(waterRow.rhs));
    rhoEnergy = norm(heatRow.rhs(:) ./ scales.energy_J_cm2_s, 2) ...
        / sqrt(numel(heatRow.rhs));
    rhoAir = norm(airRow.rhs(:) ./ scales.dry_air_g_cm2_s, 2) ...
        / sqrt(numel(airRow.rhs));
    phi = 0.5 * (rhoWater^2 + rhoEnergy^2 + rhoAir^2);
end

% -----------------------------------------------------------------------------
function failure = inversionFailure(inversion)
% Return the first protected U-to-T inversion failure for telemetry.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    failure = '';
    if inversion.allConverged || ~isfield(inversion, 'convergedMask') ...
            || ~isfield(inversion, 'failureReason')
        return;
    end
    first = find(~inversion.convergedMask, 1, 'first');
    if ~isempty(first)
        failure = inversion.failureReason{first};
    end
end

% -----------------------------------------------------------------------------
function built = enforceModflowHydrostaticConstraint(built, boundaryRequest)
% Reinstall an external MODFLOW head after a coupled correction map.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    if isempty(boundaryRequest) || ~isfield(boundaryRequest, 'bottom') ...
            || ~isfield(boundaryRequest.bottom, 'kind') ...
            || ~strcmp(boundaryRequest.bottom.kind, 'modflow_fixed_head')
        return;
    end
    modflow = boundaryRequest.bottom.modflow;
    node = modflow.active_node;
    h = built.workspace.step.trial.h(:);
    assert(node >= 1 && node <= numel(h) && node == floor(node) ...
            && isfinite(boundaryRequest.bottom.head_cm), ...
        'attemptSolveSoilEquations:badModflowConstraint', ...
        'MODFLOW fixed-head request has an invalid active node or head.');
    h(node) = boundaryRequest.bottom.head_cm;
    dz = built.domain.dz(:);
    for i = node:-1:2
        h(i-1) = h(i) + dz(i-1);
    end
    built.workspace.step.trial.h = reshape( ...
        h, size(built.workspace.step.trial.h));
end

function request = applyRealizedSurfaceEvaporation( ...
        request, waterRequest, waterFlux)
% Replace potential latent cooling by the dry branch's realized evaporation.
%
% The surface-energy response initially contains SCOPE's potential bare-soil
% latent loss. Only the residual-water active set can make actual evaporation
% smaller. Its water reaction is available before the heat row is assembled,
% so the unmet fraction is returned to the soil energy flux and its latent
% temperature tangent is removed in the same Picard map.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    if isempty(request) || ~isstruct(request) || ~isfield(request, 'top') ...
            || ~strcmp(request.top.kind, 'surface_energy_response') ...
            || ~hasAtmosphericSurface(waterRequest) ...
            || ~strcmp(waterRequest.top.surface.regime, ...
                'dry_evaporation_limited')
        return;
    end
    demand = requestFieldScalar( ...
        request.top, 'soil_evaporation_demand_cm_s', NaN);
    latent = requestFieldScalar( ...
        request.top, 'soil_evaporation_latent_heat_flux_W_m2', NaN);
    dLatent = requestFieldScalar(request.top, ...
        'soil_evaporation_latent_heat_flux_temperature_derivative_W_m2_K', ...
        NaN);
    if ~(isfinite(demand) && demand > 0 && isfinite(latent) && latent >= 0 ...
            && isfinite(dLatent) && dLatent >= 0 ...
            && isfield(waterFlux, 'QMT') && isfinite(waterFlux.QMT))
        return;
    end
    supply = waterRequest.top.liquid_supply_cm_s;
    actual = min(max(supply + waterFlux.QMT, 0), demand);
    unmetFraction = min(max((demand - actual) / demand, 0), 1);
    request.top.requested_flux_into_domain_J_cm2_s = ...
        request.top.requested_flux_into_domain_J_cm2_s ...
        + unmetFraction * latent * 1e-4;
    request.top.flux_temperature_derivative_J_cm2_s_K = ...
        request.top.flux_temperature_derivative_J_cm2_s_K ...
        + unmetFraction * dLatent * 1e-4;
    assert(request.top.flux_temperature_derivative_J_cm2_s_K <= 1e-12, ...
        'attemptSolveSoilEquations:destabilizingDryEnergyResponse', ...
        'The dry-surface heat response must retain a nonpositive tangent.');
    request.top.value_source = ...
        'coupled_surface_energy_response_with_realized_evaporation';
end

function value = requestFieldScalar(container, name, fallback)
% Read optional finite scalar boundary metadata.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    value = fallback;
    if isfield(container, name) && isscalar(container.(name)) ...
            && isfinite(container.(name))
        value = double(container.(name));
    end
end

function value = formatSurfaceTrace(decision, iteration)
% Compact audit string for diagnosing active-set cycling on rejected trials.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = sprintf(['k%d:%s>%s,v%d,proj%d,h=%.9g,', ...
        'b=%.9g,qin=%.9g'], iteration, decision.regime, ...
        decision.next_regime, decision.valid, ...
        decision.trial_dry_head_projected, decision.psi_top_cm, ...
        decision.trial_balance_rate_cm_s, ...
        decision.q_in_actual_cm_s);
end

function volume = controlVolumeWeights(dz)
% Return the bottom-up lumped nodal control-volume widths.
% Editor: Jizixing Zhu; Date: 2026-07-15.
    dz = dz(:);
    volume = [dz(1)/2; (dz(1:end-1)+dz(2:end))/2; dz(end)/2];
end

% -----------------------------------------------------------------------------
function nodal = lumpEndpointState(endpoint, dz)
% Integrate an NL-by-2 endpoint density into conservative nodal mass.
% Editor: Jizixing Zhu; Date: 2026-08-12.
    dz = dz(:);
    NL = numel(dz);
    assert(isequal(size(endpoint), [NL,2]) ...
            && all(isfinite(endpoint(:))) && all(isfinite(dz)) ...
            && all(dz > 0), ...
        'attemptSolveSoilEquations:badEndpointLumping', ...
        'Endpoint state must be finite NL-by-2 and dz finite positive.');
    nodal = zeros(NL+1,1);
    nodal(1:NL) = nodal(1:NL) + endpoint(:,1) .* dz / 2;
    nodal(2:NL+1) = nodal(2:NL+1) + endpoint(:,2) .* dz / 2;
end

function x = solveBank(bank, rhs, activeRange)
% Solve one compact tridiagonal correction block.
% Editor: Jizixing Zhu; Date: 2026-07-14.
% Change: 2026-07-20 Jizixing Zhu -- return zero corrections outside a constrained
%         MODFLOW Richards range instead of leaking untouched RHS entries.
    n = numel(rhs);
    if nargin < 3 || isempty(activeRange)
        activeRange = [1,n];
    end
    solved = solver.solveTridiagonal(bank.lower, bank.matrix(:,1), ...
        bank.matrix(:,2), rhs(:), activeRange);
    x = zeros(n,1);
    x(activeRange(1):activeRange(2)) = ...
        solved(activeRange(1):activeRange(2));
end

% -----------------------------------------------------------------------------
function first = modflowFirstActiveNode(boundaryRequest, NN)
% Return the first Richards unknown; MODFLOW owns the saturated continuation.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    first = 1;
    if isempty(boundaryRequest) || ~isfield(boundaryRequest, 'bottom') ...
            || ~isfield(boundaryRequest.bottom, 'kind') ...
            || ~strcmp(boundaryRequest.bottom.kind, 'modflow_fixed_head')
        return;
    end
    first = boundaryRequest.bottom.modflow.active_node;
    assert(isscalar(first) && isfinite(first) && first == floor(first) ...
            && first >= 1 && first <= NN, ...
        'attemptSolveSoilEquations:badModflowActiveNode', ...
        'MODFLOW active Richards node must be an integer in 1:NN.');
end

% -----------------------------------------------------------------------------
function [trial, inversion] = recoverTrialTemperature(built, refreshOpts)
% Recover the compatibility temperature mirror from prognostic (h,U,P_g).
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-17 Jizixing Zhu -- enable the safeguarded local Newton
%           fast path; difficult states retain the full bracketed fallback.
%   Change: 2026-08-20 Jizixing Zhu -- return only the small trial-state
%           container so temperature recovery does not copy the full built
%           model on every nonlinear iteration.
%   Change: 2026-08-20 Codex -- capture only the fields used by the
%           thermodynamic refresh closure. This prevents every U-to-T
%           inversion from retaining/traversing the multi-year forcing and
%           diagnostics trees carried by the full model container.

    trial = built.workspace.step.trial;
    tr = trial;
    dt_s = refreshOpts.dt_s;
    L_f = refreshOpts.L_f;
    Tmin = getNumDef(refreshOpts, 'Tmin', 0);
    thermalTemplate = compactInternalEnergyTemplate(built);
    refreshState = @(h,T,P,Lf) soil.thermal.refreshInternalEnergyState( ...
        thermalTemplate, h, T, P, dt_s, Lf, Tmin);
    [T, inversion] = soil.thermal.invertTemperatureFromInternalEnergy( ...
        tr.U(:), tr.h(:), tr.p_gas(:), tr.T(:), refreshState, ...
        built.domain.dz(:), L_f, struct('T_min', -60, 'T_max', 90, ...
            'use_local_fast_path', true));
    if inversion.allConverged
        trial.T = reshape(T, size(trial.T));
    end
end

% -----------------------------------------------------------------------------
function template = compactInternalEnergyTemplate(built)
% Retain the exact production-kernel inputs, excluding unrelated run services.
% Editor: Jizixing Zhu; Date: 2026-08-20.
    template = struct();
    template.domain = built.domain;
    template.control = struct('switches', built.control.switches);
    template.parameters = struct('soil', built.parameters.soil);
    template.Y = struct( ...
        'h', built.Y.h, ...
        'T', built.Y.T, ...
        'p_gas', built.Y.p_gas, ...
        'Theta_L', built.Y.Theta_L, ...
        'Theta_I', built.Y.Theta_I);
    template.workspace = struct( ...
        'step', struct('trial', built.workspace.step.trial), ...
        'derived', built.workspace.derived);
end

% -----------------------------------------------------------------------------
function [accepted, inversion, scale, airStatus, correctionAccepted] = ...
        applyCoupledCorrectionWithLineSearch( ...
            built, h0, P0, U0, deltaH, deltaP, deltaU, refreshOpts, ...
            heatEnabled, airEnabled, referenceAirDensity, referenceAirStorage, ...
            dryMass0, deltaDryMass, dryMassActive)
% Find the largest coupled correction inside the physical feasible domain.
%
%   Scaling the complete (dh,dU,dP_g) correction preserves its fixed point
%   and prevents a strongly nonlinear Richards update from pushing the local
%   thermodynamic state outside either the documented [-60,90] degC
%   inversion bracket or the nonnegative dry-air density/storage domain.
%   The conservative mass coordinate supplies the fraction-to-boundary
%   condition. The constitutive check admits the closed rho_a=0 boundary but
%   rejects every negative density/storage state. The raw, unscaled defects
%   remain the convergence metrics.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-08-12 Jizixing Zhu -- extend the coupled line search to
%           enforce dry-air density/storage admissibility with or without heat.

    scale = 1;
    minimumScale = 2^-12;
    remainingFraction = 0;
    accepted = built;
    inversion = struct('allConverged', ~heatEnabled, ...
        'convergedMask', repmat(~heatEnabled, numel(h0), 1), ...
        'failureReason', {repmat({''}, numel(h0), 1)});
    airStatus = disabledDryAirAdmissibility();
    correctionAccepted = false;
    while scale >= minimumScale
        dryMassScale = max(max(abs(dryMass0)), 1);
        targetDryMass = dryMass0 + scale * deltaDryMass;
        massCoordinateFeasible = ~airEnabled ...
            || all(targetDryMass >= -128 * eps(dryMassScale));
        if ~massCoordinateFeasible
            scale = scale / 2;
            continue;
        end
        candidate = built;
        candidate.workspace.step.trial.h = reshape( ...
            h0 + scale * deltaH, size(candidate.workspace.step.trial.h));
        candidate.workspace.step.trial.p_gas = reshape( ...
            P0 + scale * deltaP, size(candidate.workspace.step.trial.p_gas));
        candidate.workspace.step.trial.U = reshape( ...
            U0 + scale * deltaU, size(candidate.workspace.step.trial.U));
        activeAtBoundary = airEnabled && scale == 1 ...
            && any(dryMassActive);
        if activeAtBoundary
            candidate = projectDryAirComponentBoundary( ...
                candidate, dryMassActive, refreshOpts, false);
        end
        if heatEnabled
            [candidate.workspace.step.trial, trialInversion] = recoverTrialTemperature( ...
                candidate, refreshOpts);
        else
            trialInversion = inversion;
        end
        trialAirStatus = disabledDryAirAdmissibility();
        if trialInversion.allConverged && airEnabled
            if activeAtBoundary
                candidate = projectDryAirComponentBoundary( ...
                    candidate, dryMassActive, refreshOpts, heatEnabled);
                if heatEnabled
                    [candidate, trialInversion] = ...
                        recoverTrialTemperature(candidate, refreshOpts);
                end
            end
            [candidate, candidateEquations] = ...
                soil.refreshSoilCoefficients(candidate, refreshOpts);
            trialAirStatus = soil.gas.evaluateDryAirAdmissibility( ...
                candidate.workspace.derived.RHODA, ...
                candidateEquations.dry_air.dryAirStorage, ...
                referenceAirDensity, referenceAirStorage, ...
                remainingFraction);
        end
        if trialInversion.allConverged && trialAirStatus.admissible
            accepted = candidate;
            inversion = trialInversion;
            airStatus = trialAirStatus;
            correctionAccepted = true;
            return;
        end
        inversion = trialInversion;
        airStatus = trialAirStatus;
        scale = scale / 2;
    end
    scale = 0;
end

% -----------------------------------------------------------------------------
function candidate = projectDryAirComponentBoundary( ...
        candidate, active, refreshOpts, repeatWithEnergy)
% Put active nodes on rho_a=0 using the production constitutive derivative.
%
%   At fixed h,T, Density_DA is affine in total gas pressure and Xaa is its
%   exact positive derivative. Therefore P_g <- P_g-rho_a/Xaa is an exact
%   projection to the vapor-only component boundary. When U is prognostic,
%   the caller alternates this projection with the protected U-to-T inversion.
%   A machine-scale inward offset prevents roundoff from reporting a negative
%   component density without imposing a physical tolerance.
%   Editor: Jizixing Zhu; Date: 2026-08-12.
    if ~any(active)
        return;
    end
    if repeatWithEnergy
        maxProjection = 12;
    else
        maxProjection = 1;
    end
    for iteration = 1:maxProjection
        [candidate, ~] = soil.refreshSoilCoefficients( ...
            candidate, refreshOpts);
        rhoa = candidate.workspace.derived.RHODA(:);
        Xaa = candidate.workspace.derived.Xaa(:);
        pressure = candidate.workspace.step.trial.p_gas(:);
        pressure(active) = pressure(active) ...
            - rhoa(active) ./ Xaa(active);
        inward = 64 * eps(max(max(abs(pressure)), 1));
        pressure(active) = pressure(active) + inward;
        candidate.workspace.step.trial.p_gas = reshape( ...
            pressure, size(candidate.workspace.step.trial.p_gas));
        if ~repeatWithEnergy
            break;
        end
        [candidate.workspace.step.trial, inversion] = recoverTrialTemperature( ...
            candidate, refreshOpts);
        if ~inversion.allConverged
            break;
        end
    end
end

% -----------------------------------------------------------------------------
function [passed, thetaIncrement, saturatedHeadIncrement, activeSetChanged] = ...
        evaluateHydrusConvergence(built, hBefore, thetaBefore, ...
        tolTheta, tolHeadSaturated, refreshOpts, firstActiveNode)
% Apply the HYDRUS unsaturated-storage and saturated/crossing-head test.
%
%   If both successive iterates are unsaturated, convergence is measured by
%   the total-water-equivalent endpoint storage. If either iterate is
%   saturated, head is the informative variable because saturated storage is
%   nearly insensitive to pressure. A legitimate crossing of h=0 is therefore
%   tested with TolH and reported, but it is not rejected merely for crossing.
%   firstActiveNode excludes an externally constrained MODFLOW saturated
%   continuation from the Richards convergence norm.
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-20 Jizixing Zhu -- exclude prescribed MODFLOW rows from the
%           nonlinear Richards convergence test.
    if nargin < 7 || isempty(firstActiveNode)
        firstActiveNode = 1;
    end
    candidate = soil.refreshSoilCoefficients(built, refreshOpts);
    hAfter = built.workspace.step.trial.h(:);
    hBefore = hBefore(:);
    dAfter = candidate.workspace.derived;
    thetaAfter = waterEquivalentEndpoints( ...
        dAfter.Theta_UU, dAfter.Theta_V, dAfter.RHOV);
    hBeforeEndpoint = [hBefore(1:end-1), hBefore(2:end)];
    hAfterEndpoint = [hAfter(1:end-1), hAfter(2:end)];
    activeMask = false(size(hBefore));
    activeMask(firstActiveNode:end) = true;
    activeEndpointMask = [activeMask(1:end-1), activeMask(2:end)];
    thetaChecked = activeEndpointMask ...
        & hBeforeEndpoint < 0 & hAfterEndpoint < 0;
    headChecked = activeMask & ~(hBefore < 0 & hAfter < 0);
    activeSetChanged = any(activeMask & ((hBefore >= 0) ~= (hAfter >= 0)));
    thetaDelta = abs(thetaAfter - thetaBefore);
    if any(thetaChecked(:))
        thetaIncrement = max(thetaDelta(thetaChecked));
    else
        thetaIncrement = 0;
    end
    if any(headChecked)
        saturatedHeadIncrement = max(abs(hAfter(headChecked) - hBefore(headChecked)));
    else
        saturatedHeadIncrement = 0;
    end
    passed = thetaIncrement < tolTheta ...
        && saturatedHeadIncrement < tolHeadSaturated;
end

% -----------------------------------------------------------------------------
function pv = picardStateFromTrial(built, includeEnergy)
% Primary coupled state offered to the switch-aware Anderson pack.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    if nargin < 2, includeEnergy = false; end
    tr = built.workspace.step.trial;
    pv = struct('hh', tr.h(:), 'TT', tr.T(:), ...
        'P_gg', tr.p_gas(:));   % [cm] / [K] / [internal pressure]
    if includeEnergy
        pv.UU = tr.U(:);          % [J cm^-3] prognostic internal energy
    end
end

% -----------------------------------------------------------------------------
function built = writeTrialFromPicardState(built, pv)
% Write the unpacked Anderson guess back to active trial unknowns.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    built.workspace.step.trial.h = reshape(pv.hh, size(built.workspace.step.trial.h));   % [cm]
    built.workspace.step.trial.T = reshape(pv.TT, size(built.workspace.step.trial.T));   % [K]
    if isfield(pv, 'P_gg')
        built.workspace.step.trial.p_gas = reshape( ...
            pv.P_gg, size(built.workspace.step.trial.p_gas)); % [internal pressure]
    end
    if isfield(pv, 'UU')
        built.workspace.step.trial.U = reshape( ...
            pv.UU, size(built.workspace.step.trial.U)); % [J cm^-3]
    end
end

% -----------------------------------------------------------------------------
function status = disabledDryAirAdmissibility()
% Return the neutral physical-gate result when dry-air physics is disabled.
% Editor: Jizixing Zhu; Date: 2026-08-12.
    status = struct('admissible', true, ...
        'minimum_density_g_cm3', NaN, ...
        'minimum_storage_g_cm3', NaN, ...
        'failure_reason', '', ...
        'remaining_fraction', 0);
end

% -----------------------------------------------------------------------------
function thetaEq = waterEquivalentEndpoints(Theta_UU, Theta_V, RHOV)
% Convert condensed and vapor storage to endpoint total-water equivalent.
% Editor: Jizixing Zhu; Date: 2026-07-12.
    NL = size(Theta_UU, 1);
    RHOV = RHOV(:);
    assert(isequal(size(Theta_UU), [NL, 2]) ...
        && isequal(size(Theta_V), [NL, 2]) && numel(RHOV) == NL + 1, ...
        'attemptSolveSoilEquations:badThetaEqShape', ...
        'Theta_UU/Theta_V must be NL-by-2 and RHOV must contain NL+1 nodes.');
    RHOV_endpoints = [RHOV(1:NL), RHOV(2:NL + 1)];
    thetaEq = Theta_UU ...
        + Theta_V .* RHOV_endpoints / equations.getPhysicalConstant('RHOL');
end

% -----------------------------------------------------------------------------
function tf = hasAtmosphericSurface(request)
% Test whether one materialized water request owns the DSTOR active set.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    tf = isstruct(request) && isfield(request, 'top') ...
        && isstruct(request.top) ...
        && isfield(request.top, 'kind') ...
        && strcmp(request.top.kind, 'atmospheric') ...
        && isfield(request.top, 'surface') ...
        && isstruct(request.top.surface);
end

% -----------------------------------------------------------------------------
function tf = shouldCheckSurfaceConstraintNow( ...
        request, surfaceHeadCm, dryHeadProjected, massToleranceCm)
%CHECKSURFACECONSTRAINTNOW  Detect kinematic active-set violations early.
%
%   Head-bound and external-supply violations do not require a converged
%   boundary reaction. Reaction-sign tests normally do; the sole exception
%   is the declared zero-head abrupt-wetting feasibility predictor. The
%   dry-head evaporation reaction is deliberately absent here.
%   Editor: Jizixing Zhu; Date: 2026-07-17.
%   Change: 2026-08-08 Jizixing Zhu -- defer the reverse supply-head test
%     after a converged zero-head reaction releases that branch.
    tf = false;
    if ~hasAtmosphericSurface(request)
        return;
    end
    surface = request.top.surface;
    dmax = surface.capacity_cm;
    headTol = max(1e-10, 1e-8 * max(1, dmax));
    fluxTol = max(1e-14, massToleranceCm / request.dt_s);
    switch char(surface.regime)
        case 'supply_limited'
            releasedByConvergedReaction = ...
                isfield(surface, 'initialization') ...
                && strcmp(surface.initialization, ...
                    'converged_reaction_release');
            tf = ~releasedByConvergedReaction ...
                && surfaceHeadCm > headTol;
        case 'storage'
            tf = surfaceHeadCm < -headTol ...
                || surfaceHeadCm > dmax + headTol;
        case 'dry_evaporation_limited'
            tf = request.top.liquid_supply_cm_s > fluxTol;
        case 'capacity_limited'
            % The abrupt-wetting predictor is explicitly a zero-head Darcy
            % feasibility probe, so its reaction is the predictor output and
            % must be inspected immediately. A prognostic capacity state is
            % instead retained until its coupled branch has converged.
            tf = isfield(surface, 'initialization') ...
                && strcmp(surface.initialization, ...
                    'wet_feasibility_predictor');
        otherwise
            error('attemptSolveSoilEquations:badSurfaceRegime', ...
                'Unknown atmospheric surface regime "%s".', ...
                char(surface.regime));
    end
end

% -----------------------------------------------------------------------------
function regime = surfaceRegime(request)
% Return the active atmospheric surface regime for audit telemetry.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    regime = '';
    if hasAtmosphericSurface(request) ...
            && isfield(request.top.surface, 'regime')
        regime = char(request.top.surface.regime);
    end
end

% -----------------------------------------------------------------------------
function v = getNumDef(nm, name, dflt)
% Read a control.numerics field with a default (keeps the diagnostic self-contained).
% Editor: Jizixing Zhu; Date: 2026-07-14.
    v = dflt;
    if isfield(nm, name) && ~isempty(nm.(name)), v = nm.(name); end
end

% -----------------------------------------------------------------------------
function enabled = switchEnabled(switches, name)
% Read one master equation switch without inventing an enabled default.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    enabled = isfield(switches, name) && switches.(name) == 1;
end

% -----------------------------------------------------------------------------
function sink = effectiveWaterEquationSink(built, rootSink)
% Add the fixed-per-substep dual-porosity exchange to the assembled water load.
% Editor: Jizixing Zhu; Date: 2026-07-14.
% Change: 2026-07-20 Jizixing Zhu -- read the trial-only exchange result instead of
%   accepted Y, so rejected trials cannot advance the immobile domain.
    sink = rootSink;
    if switchEnabled(built.control.switches, 'lDualPor') ...
            && isfield(built.workspace.step, 'trial_aux') ...
            && isfield(built.workspace.step.trial_aux, 'dual_porosity') ...
            && isfield(built.workspace.step.trial_aux.dual_porosity, ...
                'SinkIm_current') ...
            && ~isempty(built.workspace.step.trial_aux.dual_porosity. ...
                SinkIm_current)
        sink = sink + built.workspace.step.trial_aux.dual_porosity. ...
            SinkIm_current;
    end
end
