function accum = accumulateWaterFluxAccum(accum, waterFlux, weight)
%ACCUMULATEWATERFLUXACCUM  Add one accepted substep's waterFlux to the hourly accumulator (WBC-2).
%
%   accum = solver.accumulateWaterFluxAccum(accum, waterFlux, weight) time-weights one ACCEPTED tile's
%   water-boundary flux into accum, mirroring solver.accumulateAcceptedSubstep's rule
%   (accum.f += flux.f * weight) for the waterFlux struct.
%
%   Inputs:
%     accum     [inout] accumulator from solver.initWaterFluxAccum.
%     waterFlux [in]    the tile's waterFlux struct (soil.water.solveWaterTrial 3rd output) [cm/s].
%     weight    [in]    [-] the accepted tile weight = dt / dtOuterActual (A1-clipped accepted dt).
%
%   Returns:
%     accum  with weight_sum and every flux field incremented by weight * waterFlux.field.
%
%   Reference: legacy solver.accumulateAcceptedSubstep (src/+solver/accumulateAcceptedSubstep.m:5-13).
%   Editor: Jizixing Zhu; Date: 2026-07-08.
    accum.weight_sum = accum.weight_sum + weight;
    fn = fieldnames(waterFlux);
    for i = 1:numel(fn)
        accum.(fn{i}) = accum.(fn{i}) + waterFlux.(fn{i}) * weight;   % [cm/s] weighted flux
    end
end
