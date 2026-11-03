function accum = initWaterFluxAccum()
%INITWATERFLUXACCUM  Zeroed hour-local accumulator for the WBC-1 waterFlux struct (WBC-2).
%
%   accum = solver.initWaterFluxAccum() returns the per-outer-step accumulator for the water-boundary
%   fluxes, mirroring solver.initHourlySubstepAccum but over the waterFlux struct fields. weight_sum is
%   the running sum of accepted-tile weights; every flux field starts at 0 [cm/s].
%
%   Returns:
%     accum  struct with weight_sum [-] and the eight waterFlux fields (QMT/QMB/AVAIL0/Throughfall/
%            Infiltration/Infiltration_act/TopFlux_potential/runoffHort) [cm/s], all 0. For production
%            atmospheric requests runoffHort is the accepted capacity-overflow rate. Hour-local
%            scratch: never scattered to Y/pending/checkpoint.
%
%   Reference: legacy solver.initHourlySubstepAccum (src/+solver/initHourlySubstepAccum.m).
%   Editor: Jizixing Zhu; Date: 2026-07-08.
%   Change: 2026-07-17 Jizixing Zhu -- document realized surface-overflow semantics.
%   Change: 2026-07-17 Jizixing Zhu -- mark potential/actual infiltration fields as
%           private legacy-kernel adapter values; acceptance and diagnostics
%           use QMT/QMB plus the surface-reservoir ledger.
    accum = struct('weight_sum', 0, ...
        'QMT', 0, 'QMB', 0, 'AVAIL0', 0, 'Throughfall', 0, ...
        'Infiltration', 0, 'Infiltration_act', 0, 'TopFlux_potential', 0, 'runoffHort', 0);
end
