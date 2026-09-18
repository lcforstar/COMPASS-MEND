function accum = initHourlySubstepAccum()
% initHourlySubstepAccum Initialize weighted hourly accumulators for substeps.

accum.weight_sum = 0;
accum.QMT = 0;
accum.QMB = 0;
accum.QET = 0;
accum.QEB = 0;
accum.RS = 0;
accum.NITT = 0;
end
