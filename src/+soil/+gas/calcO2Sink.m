function [R_O2_sink, nitrificationSink] = calcO2Sink(Snit, respirationOrNL, NL)
%CALCO2SINK Convert nitrification rate to stoichiometric oxygen demand.
%
%   Snit is the positive NH4-N consumption rate produced by
%   soil.nitrogen.chemistry.nitrification [ug N cm^-3 bulk s^-1, NL x 2].
%   Complete nitrification consumes two moles of O2 per mole of N:
%
%     NH4+ + 2 O2 -> NO3- + H2O + 2 H+.
%
%   The returned rate is [mol O2 cm^-3 bulk s^-1, NL x 1]. Soil organic-N
%   mineralization Smin is deliberately excluded: Smin is a calibrated
%   organic-N-to-NH4 source, not a carbon-respiration rate, so converting it
%   with a carbon molar mass has no valid unit or stoichiometric basis.
%
%   Editor: Jizixing Zhu; Date: 2026-05-16.
%   Change: 2026-08-16 Jizixing Zhu -- remove the dimensionally invalid Smin
%   oxygen demand and reject negative nitrification instead of taking abs().

% Change: 2026-10-23 Wenjie Liu -- retain the original two-input nitrogen
% conversion; optional (Snit, respiratoryC, NL) adds a named carbon source
% [ug C cm^-3 bulk s^-1, NL-by-1 or NL-by-2], at one mol O2 per mol C.
% Optional nitrificationSink is the same layer component before respiratory demand is added.
    respiratoryC = [];
    if nargin == 2
        NL = respirationOrNL;
    else
        respiratoryC = respirationOrNL;
    end
    assert(isequal(size(Snit), [NL, 2]) && all(isfinite(Snit(:))), ...
        'calcO2Sink:badNitrificationShape', ...
        'Snit must be one finite NL-by-2 endpoint rate field.');
    negativeTolerance = 1e-18;
    assert(min(Snit(:)) >= -negativeTolerance, ...
        'calcO2Sink:negativeNitrification', ...
        'Nitrification oxygen demand requires nonnegative Snit.');
    Snit = max(Snit, 0);
    Snit_avg = mean(Snit, 2);                   % [ug N cm^-3 s^-1]
    nitrificationSink = Snit_avg * 1e-6 / 14 * 2;      % [mol O2 cm^-3 s^-1]
    R_O2_sink = nitrificationSink;
    if ~isempty(respiratoryC)
        assert((isequal(size(respiratoryC),[NL,1]) ...
                || isequal(size(respiratoryC),[NL,2])) ...
                && all(isfinite(respiratoryC(:))) && all(respiratoryC(:)>=0), ...
            'calcO2Sink:badRespiration','Respiratory C source must be finite and nonnegative.');
        R_O2_sink = R_O2_sink + mean(respiratoryC,2)*1e-6/12;
    end
end

% Change: 2026-10-23 Wenjie Liu -- restore the COMPASS-COMPASS_MATLAB
% oxygen transport/source boundary implementation; retain current carrier physics.
