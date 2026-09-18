function accum = accumulateAcceptedSubstep(accum, Delt_t, Delt_t_outer, KT, ...
    QMT, QMB, QET, QEB, RS_trial, IterUsed, cfg)
% accumulateAcceptedSubstep Accumulate accepted internal substep fluxes.

sub_weight = Delt_t / Delt_t_outer;
accum.weight_sum = accum.weight_sum + sub_weight;
accum.QMT = accum.QMT + QMT(KT) * sub_weight;
accum.QMB = accum.QMB + QMB(KT) * sub_weight;
if cfg.lThmrlefc == 1
    accum.QET = accum.QET + QET(KT) * sub_weight;
    accum.QEB = accum.QEB + QEB(KT) * sub_weight;
end
accum.RS = accum.RS + RS_trial * sub_weight;
accum.NITT = accum.NITT + IterUsed;
end
