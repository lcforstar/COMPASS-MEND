function built = restoreSubstepState(built, snapshot)
%RESTORESUBSTEPSTATE  Restore a snapshot after a rejected solver trial.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.

    built.Y = snapshot.Y;
    built.workspace.step = snapshot.step;
    built.workspace.surface = snapshot.surface;
    built.workspace.scope_step = snapshot.scope_step;
    built.workspace.solver.SAVE = snapshot.solver_SAVE;
end
