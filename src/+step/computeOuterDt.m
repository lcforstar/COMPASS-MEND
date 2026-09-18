function dt_s = computeOuterDt(dt_config_s, time_control, elapsed_s)
%COMPUTEOUTERDT  Clip the final outer step to the configured end time.
%
%   Inputs:
%     dt_config_s  [s] Configured outer-step duration.
%     time_control [in] Time configuration, optionally containing TEND [s].
%     elapsed_s    [s] Accepted cumulative elapsed time before this step.
%
%   Return:
%     dt_s         [s] Actual duration proposed for this outer step.
%
%   Interior steps retain dt_config_s. Only the final positive remainder is
%   shorter, matching the legacy COMPASS final-step clipping contract.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    dt_s = dt_config_s;
    if isfield(time_control, 'TEND')
        remaining_s = time_control.TEND - elapsed_s;
        if remaining_s > 0 && remaining_s < dt_config_s
            dt_s = remaining_s;
        end
    end
end

