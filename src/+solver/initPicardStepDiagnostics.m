function diag = initPicardStepDiagnostics(nit)
%INITPICARDSTEPDIAGNOSTICS Allocate diagnostics for one outer time step.

diag = struct();
diag.water_error = NaN(nit, 1);
diag.energy_error = NaN(nit, 1);
diag.hh_surface = NaN(nit, 1);
diag.TT_surface = NaN(nit, 1);
diag.aa_applied = false(nit, 1);
diag.bc_switch = false(nit, 1);
diag.stage = repmat({''}, nit, 1);
end
