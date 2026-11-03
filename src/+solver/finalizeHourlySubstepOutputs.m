function [QMT, QMB, QET, QEB, RS] = finalizeHourlySubstepOutputs( ...
    accum, KT, cfg, QMT, QMB, QET, QEB, RS)
% finalizeHourlySubstepOutputs Write accepted substep averages to hourly arrays.

if accum.weight_sum > 0
    QMT(KT) = accum.QMT / accum.weight_sum;
    QMB(KT) = accum.QMB / accum.weight_sum;
    if cfg.lThmrlefc == 1
        QET(KT) = accum.QET / accum.weight_sum;
        QEB(KT) = accum.QEB / accum.weight_sum;
    end
    RS(KT) = accum.RS / accum.weight_sum;
else
    QMT(KT) = 0;
    QMB(KT) = 0;
    if cfg.lThmrlefc == 1
        QET(KT) = 0;
        QEB(KT) = 0;
    end
    RS(KT) = 0;
end
end
