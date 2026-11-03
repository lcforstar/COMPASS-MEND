function hourly = finalizeWaterFluxAccum(accum)
%FINALIZEWATERFLUXACCUM  Weight-normalize the hourly waterFlux accumulator (WBC-2).
%
%   hourly = solver.finalizeWaterFluxAccum(accum) divides each accumulated flux by weight_sum to give
%   the time-weighted average rate, mirroring solver.finalizeHourlySubstepOutputs; when nothing was
%   accepted (weight_sum <= 0) every field is 0, so a fully-rejected hour yields zeros (no failed-tile
%   flux leaks into the hourly output). The field names are IDENTICAL to the waterFlux struct (no
%   alias), so WBC-3 reads QMB/Infiltration_act/... directly.
%
%   Inputs:
%     accum  [in]  accumulator from solver.initWaterFluxAccum / solver.accumulateWaterFluxAccum.
%
%   Returns:
%     hourly struct with the eight waterFlux fields [cm/s] as time-weighted average rates (or 0).
%
%   Reference: legacy solver.finalizeHourlySubstepOutputs (src/+solver/finalizeHourlySubstepOutputs.m:5-21).
%   Editor: Jizixing Zhu; Date: 2026-07-08.
    fn = setdiff(fieldnames(accum), {'weight_sum'});
    hourly = struct();
    for i = 1:numel(fn)
        if accum.weight_sum > 0
            hourly.(fn{i}) = accum.(fn{i}) / accum.weight_sum;   % [cm/s] weighted-average rate
        else
            hourly.(fn{i}) = 0;
        end
    end
end
