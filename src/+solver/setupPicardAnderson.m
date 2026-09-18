function setup = setupPicardAnderson(built)
%SETUPPICARDANDERSON  Per-substep Picard diagnostics + Anderson acceleration SETUP (Layer A / A4).
%
%   setup = solver.setupPicardAnderson(built) builds the per-substep convergence-diagnostics seeds
%   and the Anderson-acceleration option bundle that the (future) Picard loop consumes. Ports the
%   SETUP block COMPASS.m:1610-1638 only. It performs NO acceleration: it does not apply an Anderson
%   guess (legacy COMPASS.m:1659-1669), does not commit a Picard map (COMPASS.m:1877-1883), does not
%   open the Picard loop, and wires no runOneStep/attemptSolve/refreshSoilCoefficients/solve. It does
%   not mutate built.
%
%   Inputs:
%     built  [in]  container. Reads control.numerics.NIT [-]; control.numerics.lAnderson [-]
%                  (optional solver toggle); control.switches.lNefc [-]; domain.NN [-]; and the
%                  post-A3 surface-node trial iterate workspace.step.trial.{h [cm], T [K], nh4, no3}
%                  at node NN.
%
%   Returns:
%     setup  struct with:
%       .picardStep    struct  solver.initPicardStepDiagnostics(NIT) output (per-iterate probes)
%       .lAnderson     [-]     logical; TRUE unless control.numerics.lAnderson is present and false
%       .aaPicard      [-]     [] (Anderson state, filled by the Picard loop)
%       .aaStateOpts   struct  scaled stacked-state options; pressure is included only when the
%                      dry-air equation and numerics.includePressureInAnderson are both enabled
%       .aaOpts        struct  Anderson options (enabled=lAnderson; memory/damping/...)
%       .hSAVE/.TSAVE          surface-node trial head [cm] / temperature [K] snapshots [1 x nCell]
%       .NSAVE/.NSAVE2         surface-node trial NH4/NO3 (only when lNefc==1; [] otherwise)
%       .IterUsed      [-]     SEED = NIT (legacy initial value); the future Picard loop overwrites
%                              it with the actual iteration count -- it is NOT a used-iteration count
%       .StepConverged [-]     false
%       .Error_Water/.Error_Energy  NaN
%       .lastNonConvergenceReason   'not evaluated'
%       .P_gg_before_air_debug      []
%
%   Method: pure construction from built; no solve, no acceleration. lAnderson defaults TRUE
%     (legacy COMPASS.m:1622: a missing flag is NOT false). The scaled stacked state
%     [hh/scale_hh; TT/scale_TT; P_gg/scale_P_gg] is used only when dry air is enabled and pressure
%     acceleration is explicitly requested. The default omits pressure for compatibility and to
%     avoid contaminating water/heat acceleration with a weakly constrained gas-pressure block.
%
%   Reference: COMPASS.m:1610-1638 (setup block only). Anderson apply/commit are at COMPASS.m:1659-1669
%     and 1877-1883 (the Picard loop, NOT ported here).
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-13 Jizixing Zhu -- include dry-air pressure in the
%           Anderson state exactly when lSoilairefc is enabled.
%   Change: 2026-07-14 Jizixing Zhu -- honor includePressureInAnderson and
%           add an explicit conditioning gate for Anderson least squares.
    num   = built.control.numerics;
    sw    = built.control.switches;
    NIT   = num.NIT;
    NN    = built.domain.NN;
    trial = built.workspace.step.trial;

    lAnderson = true;   % legacy default (COMPASS.m:1622): a missing flag is NOT false
    if isfield(num, 'lAnderson') && ~isempty(num.lAnderson)
        lAnderson = logical(num.lAnderson);
    end

    setup = struct();
    setup.picardStep  = solver.initPicardStepDiagnostics(NIT);
    setup.lAnderson   = lAnderson;
    setup.aaPicard    = [];
    dryAirEnabled = isfield(sw, 'lSoilairefc') && sw.lSoilairefc == 1;
    pressureRequested = isfield(num, 'includePressureInAnderson') && ...
        ~isempty(num.includePressureInAnderson) && logical(num.includePressureInAnderson);
    includePressure = dryAirEnabled && pressureRequested;
    setup.aaStateOpts = struct('includePressure', includePressure, ...
        'scale_hh', 100, 'scale_TT', 10, 'scale_P_gg', 1e5);
    setup.aaOpts = struct('enabled', lAnderson, ...
        'memory', 3, ...
        'minHistory', 3, ...
        'damping', 0.5, ...
        'regularization', 1e-10, ...
        'minSingularRatio', 1e-8, ...
        'maxGammaNorm', 1e4, ...
        'maxStepNorm', 2.0, ...
        'safeguardFactor', 1.2);

    % Per-substep convergence-diagnostics seeds at the surface node NN (regional-safe row per cell).
    setup.hSAVE = trial.h(NN, :);   % [cm] surface trial head
    setup.TSAVE = trial.T(NN, :);   % [K]  surface trial temperature
    if isfield(sw, 'lNefc') && sw.lNefc == 1
        setup.NSAVE  = trial.nh4(NN, :);   % surface trial NH4 (legacy C_NH4_new)
        setup.NSAVE2 = trial.no3(NN, :);   % surface trial NO3 (legacy C_NO3_new)
    else
        setup.NSAVE  = [];
        setup.NSAVE2 = [];
    end
    setup.IterUsed                 = NIT;    % seed = NIT (legacy initial value); overwritten by loop
    setup.StepConverged            = false;
    setup.Error_Water              = NaN;
    setup.Error_Energy             = NaN;
    setup.lastNonConvergenceReason = 'not evaluated';
    setup.P_gg_before_air_debug    = [];
end
