function report = finalizeReport(report, runtime, artifacts)
%FINALIZEREPORT  Attach accepted runtime and file-artifact status to report.
%
%   Inputs:
%     report    [inout] Fixed outer-step report.
%     runtime   [in] Accepted runtime clock state.
%     artifacts [in] Result from io.writePeriodicStepArtifacts.
%
%   Return:
%     report    [out] Report synchronized with accepted calendar and I/O.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    report.current_day = runtime.step.current_day;
    report.current_hour = runtime.step.current_hour;
    report.checkpoint_file = artifacts.checkpoint_file;
    report.partial_output_written = artifacts.partial_output_written;
end
