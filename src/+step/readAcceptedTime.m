function elapsed_s = readAcceptedTime(runtime)
%READACCEPTEDTIME  Read elapsed time at the last accepted outer step.
%
%   Input:
%     runtime   [in] Runtime service state.
%
%   Return:
%     elapsed_s [s] Accepted cumulative elapsed time; zero before step one.
%
%   This accessor does not infer time from kt and never reads trial state.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    elapsed_s = 0;
    if isfield(runtime, 'step') && isfield(runtime.step, 'TIME')
        elapsed_s = runtime.step.TIME;
    end
end

