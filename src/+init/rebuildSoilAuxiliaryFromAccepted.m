function built = rebuildSoilAuxiliaryFromAccepted(built, resumeKt)
%REBUILDSOILAUXILIARYFROMACCEPTED  Rebuild soil aux/cache from accepted Y only.
%
%   built = init.rebuildSoilAuxiliaryFromAccepted(built) reallocates disposable
%   solver workspace and recomputes workspace.derived from accepted built.Y.
%   Unlike init.initializeSoilAux, this path is safe on checkpoint resume:
%     - does NOT call soil.water.initializeFreezingHead (that mutates Y.h);
%     - does NOT overwrite a finite prognostic Y.U with a kernel re-evaluation.
%
%   built = init.rebuildSoilAuxiliaryFromAccepted(built, resumeKt) records the
%   completed checkpoint KT for documentation and future callers. Thermal
%   bottom-boundary refresh is intentionally NOT performed here: that field is
%   time-indexed and is rebuilt at the start of the next runOneStep with the
%   real outer-step index (ckpt.kt+1). Using KT=1 would inject the first-hour
%   harmonic bottom temperature into a winter resume.
%
%   Pure rebuild semantics (audit 2026-07-19):
%     B1 water/phase: dt_s = 0 so Dall'Amico time-relaxation cannot advance
%       ice/liquid relative to accepted Y.Theta_I / SFCC endpoint state.
%     B2 thermal-gas transport: dt_s = Delt_t (>0 required by Density_DA /
%       thermal liquid dispersion denominators). This only sets transport
%       coefficient scale; it is not a model time advance.
%     Thermal boundary: not refreshed (refreshThermalBoundary equivalent false).
%     Equation bank: built with the same transport dt only when missing U must
%       be materialised for old checkpoints.
%
%   Inputs:
%     built    [inout] model container with accepted built.Y and control/domain
%     resumeKt [in, optional] [-] completed checkpoint KT (informational)
%
%   Returns:
%     built [out] same prognostic Y (see invariance rules), with
%                 workspace.trial/previous warm-started from Y and
%                 workspace.derived rebuilt without time-step physics
%
%   Y invariance:
%     - Current finite-U checkpoint: isequaln(built.Y, Y0) after rebuild.
%     - Missing / all-nonfinite U: only Y.U may be filled; rest of Y equals Y0.
%     - Mixed finite/nonfinite U: error (no silent retain).
%
%   Ownership:
%     Y          prognostic (immutable except optional missing-U fill)
%     workspace  disposable aux/cache (rebuilt)
%
%   Editor: Jizixing Zhu; Date: 2026-07-19.
%   Change: 2026-07-19 Jizixing Zhu -- pure B1 dt_s=0, skip Tbh refresh, full Y
%     invariance, refuse partial-NaN U; optional resumeKt is not used to fake KT=1.
%   Change: 2026-07-20 Jizixing Zhu -- rename the local equation bank so it
%     cannot shadow the equations physical-constants package.

    if ~isstruct(built) || ~isfield(built, 'Y')
        error('rebuildSoilAuxiliaryFromAccepted:badBuilt', ...
            'built.Y is required.');
    end
    if nargin < 2 || isempty(resumeKt)
        resumeKt = NaN;
    end
    if ~(isnumeric(resumeKt) && isscalar(resumeKt))
        error('rebuildSoilAuxiliaryFromAccepted:badResumeKt', ...
            'resumeKt must be a numeric scalar when provided.');
    end

    Y0 = built.Y;
    uPolicy = classifyAcceptedU(Y0);

    built = init.allocateWorkspace(built);

    L_f = 0;
    if isfield(built.control, 'switches') ...
            && isfield(built.control.switches, 'lFreezing') ...
            && built.control.switches.lFreezing == 1
        L_f = equations.getPhysicalConstant('lambda_fus');
        override = str2double(io.getEnv('LFUSION_OVERRIDE'));
        if isfinite(override)
            L_f = override;
        end
    end
    built.workspace.step.L_f = L_f;

    % B1: instantaneous water/phase cache from accepted trial=Y (dt_s=0).
    built = soil.refreshWaterAux(built, 0, L_f, 0);

    % B2: transport coefficients require a positive scale dt.
    dtTransport = built.control.time.Delt_t;
    assert(isscalar(dtTransport) && isfinite(dtTransport) && dtTransport > 0, ...
        'rebuildSoilAuxiliaryFromAccepted:badDt', ...
        'built.control.time.Delt_t must be a finite positive outer step [s].');
    built = soil.refreshThermalGasTransport(built, dtTransport);
    % Deliberately omit soil.refreshThermalBoundaryAuxiliary: next outer step
    % owns the time-indexed Tbh refresh with the true KT.

    if isfinite(resumeKt)
        % Stash only for diagnostics/tests; not read by physics.
        built.workspace.resume_completed_kt = resumeKt;
    end

    if uPolicy.needsFill
        opts = struct( ...
            'dt_s', dtTransport, ...
            'L_f', L_f, ...
            'Tmin', 0, ...
            'DSTMAX', 0, ...
            'KT', 1, ...
            'NoTime_KT', 1, ...
            'Srt', built.workspace.flux.Srt, ...
            'airHeatImpedance', 0, ...
            'refreshThermalBoundary', false);
        equationBank = soil.buildEquationCoefficients(built, opts);
        if ~isfield(equationBank, 'heat') ...
                || ~isfield(equationBank.heat, 'energyStorage') ...
                || isempty(equationBank.heat.energyStorage)
            energyBuilt = built;
            energyBuilt.control.switches.lThmrlefc = 1;
            energyEquations = soil.buildEquationCoefficients(energyBuilt, opts);
            equationBank.heat = energyEquations.heat;
        end
        assert(isfield(equationBank.heat, 'energyStorage') ...
            && ~isempty(equationBank.heat.energyStorage), ...
            'rebuildSoilAuxiliaryFromAccepted:missingInitialEnergy', ...
            'Cannot materialize missing Y.U without energyStorage.');
        built = init.setPrognosticInternalEnergy( ...
            built, equationBank.heat.energyStorage);
        assertYEqualExceptU(Y0, built.Y);
    else
        built.workspace.step.trial.U = built.Y.U;
        built.workspace.step.previous.U = built.Y.U;
        assert(isequaln(built.Y, Y0), ...
            'rebuildSoilAuxiliaryFromAccepted:mutatedY', ...
            'Pure aux rebuild must not mutate any field of accepted Y.');
    end
end

function policy = classifyAcceptedU(Y)
% Decide whether Y.U may be filled, must be retained, or is illegal mixed.
    policy = struct('needsFill', false);
    if ~isfield(Y, 'U') || isempty(Y.U)
        policy.needsFill = true;
        return;
    end
    u = double(Y.U(:));
    isFin = isfinite(u);
    if all(~isFin)
        policy.needsFill = true;
        return;
    end
    if any(~isFin)
        error('rebuildSoilAuxiliaryFromAccepted:partialU', ...
            ['Y.U has mixed finite and non-finite values (%d finite / %d total). ', ...
             'Refuse silent retain or partial fill on resume.'], ...
            nnz(isFin), numel(u));
    end
    policy.needsFill = false;
end

function assertYEqualExceptU(Y0, Y1)
% Require whole-tree equality except the allowed missing-U fill of Y.U.
    f0 = fieldnames(Y0);
    f1 = fieldnames(Y1);
    extra = setdiff(f1, f0);
    if ~isempty(extra) && ~(isscalar(extra) && strcmp(extra{1}, 'U'))
        error('rebuildSoilAuxiliaryFromAccepted:mutatedY', ...
            'Pure aux rebuild added unexpected Y fields: %s', strjoin(extra, ', '));
    end
    missing = setdiff(f0, f1);
    if ~isempty(missing)
        error('rebuildSoilAuxiliaryFromAccepted:mutatedY', ...
            'Pure aux rebuild dropped Y fields: %s', strjoin(missing, ', '));
    end
    for i = 1:numel(f0)
        name = f0{i};
        if strcmp(name, 'U')
            continue;
        end
        if ~isequaln(Y0.(name), Y1.(name))
            error('rebuildSoilAuxiliaryFromAccepted:mutatedY', ...
                'Pure aux rebuild mutated Y.%s while filling missing U.', name);
        end
    end
end
