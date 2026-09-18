function [summary, history, reason, aaUsed] = finalizePicardDiagnostics( ...
    summary, history, kt, stepDiag, saveHistory, converged, iterations, ...
    waterError, energyError, cfg, deltT, fallbackReason)
%FINALIZEPICARDDIAGNOSTICS Store and print compact diagnostics for one time step.

idx = find(~isnan(stepDiag.water_error) | ~isnan(stepDiag.energy_error) | ~cellfun('isempty', stepDiag.stage));
if isempty(idx)
    idx = 1:min(iterations, cfg.NIT);
end
lastIdx = idx(end);

reason = stepDiag.stage{lastIdx};
if isempty(reason)
    reason = fallbackReason;
end
aaUsed = nnz(stepDiag.aa_applied(idx));

summary(kt).converged = converged;
summary(kt).iterations = iterations;
summary(kt).water_error = waterError;
summary(kt).energy_error = energyError;
summary(kt).reason = reason;
summary(kt).aa_used = aaUsed;

if saveHistory
    history(kt).water_error = stepDiag.water_error(idx);
    history(kt).energy_error = stepDiag.energy_error(idx);
    history(kt).hh_surface = stepDiag.hh_surface(idx);
    history(kt).TT_surface = stepDiag.TT_surface(idx);
    history(kt).aa_applied = stepDiag.aa_applied(idx);
    history(kt).bc_switch = stepDiag.bc_switch(idx);
    history(kt).stage = stepDiag.stage(idx);
end

% 注：'NOT-CONV' = Picard 内部不收敛（dt 缩半重试，可救回），不是程序崩溃。
% 真正的 FATAL = 最小步 dt 仍不收敛（src/COMPASS.m 里的 ROLLBACK-MIN warning）。
if converged
    status = 'OK';
else
    status = 'NOT-CONV';
end
fprintf('[Picard] KT=%d %s Iter=%d/%d Water=%9.3e Energy=%9.3e AA=%d Reason=%s\n', ...
    kt, status, iterations, cfg.NIT, waterError, energyError, aaUsed, reason);
end
