function snapshot = snapshotSubstepState(built)
%SNAPSHOTSUBSTEPSTATE  Capture the mutable state of one solver trial.
%
%   The snapshot covers accepted Y, step anchors/trial scratch, surface and
%   SCOPE scratch, and the shared boundary-flux SAVE array. Derived
%   coefficients are excluded because they are a recomputable Picard cache.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.

    snapshot = struct();
    snapshot.Y = built.Y;
    snapshot.step = built.workspace.step;
    snapshot.surface = built.workspace.surface;
    snapshot.scope_step = built.workspace.scope_step;
    snapshot.solver_SAVE = built.workspace.solver.SAVE;
end
