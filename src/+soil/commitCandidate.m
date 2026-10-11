function built = commitCandidate(built, candidate)
%COMMITCANDIDATE Install the accepted coupled-soil endpoint into legacy state.
%
%   Inputs:
%     built     [struct] Legacy model shell before soil acceptance.
%     candidate [struct] Soil-owned t^(n+1) trial and derived endpoint.
%
%   Returns:
%     built     [struct] Shell whose accepted soil state equals candidate.
%
%   Physical contract:
%     `trial` is the converged discrete h/U/P_g solution.  `derived` supplies
%     theta liquid/ice layers evaluated at that same endpoint.  Installing
%     both together preserves the soil state-vector consistency required by
%     water and internal-energy conservation; this routine has no boundary
%     calculation and no diagnostic side effect.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
% Change: 2026-07-20 Jizixing Zhu -- commit soil-N pools, O2, and N solver memory atomically with soil state.
% Change: 2026-07-20 Jizixing Zhu -- merge only soil-owned trial fields so surface
%   and vegetation scratch cannot enter through the soil candidate.
% Change: 2026-08-29 Codex -- commit native AOM/DOM/MEMC/MEND state in the
%   same accepted transaction as water, heat, gas, and mineral N.
% Change: 2026-09-12 Codex -- retain accepted post-root DIN availability as
%   diagnostic-only runtime memory alongside the root reservation ledger.
% Change: 2026-09-14 Codex -- commit O2 capacity history with its concentration.
% Change: 2026-09-21 Wenjie Liu -- commit the passive fertilizer-N source
%   label atomically with the total soil/organic candidate.
    assert(isstruct(candidate) && isfield(candidate, 'trial') ...
        && isfield(candidate, 'derived'), ...
        'soil:commitCandidate:invalidCandidate', ...
        'Soil candidate must contain trial and derived endpoint fields.');
    trialNames = fieldnames(candidate.trial);
    for i = 1:numel(trialNames)
        name = trialNames{i};
        built.workspace.step.trial.(name) = candidate.trial.(name);
    end
    built.workspace.derived = candidate.derived;
    built = solver.carrySubstepState(built);
    if isfield(candidate, 'nitrogen') && ~isempty(candidate.nitrogen)
        nitrogen = candidate.nitrogen;
        required = {'state','trial','ice','oxygen','diagnostics'};
        for i = 1:numel(required)
            assert(isfield(nitrogen, required{i}), ...
                'soil:commitCandidate:incompleteNitrogenCandidate', ...
                'Accepted soil-N candidate lacks %s.', required{i});
        end
        built.Y.soil_nitrogen = nitrogen.state;
        built.Y.C_NH4_ice = nitrogen.ice.nh4;
        built.Y.C_NO3_ice = nitrogen.ice.no3;
        built.Y.soil_gas.C_O2 = nitrogen.oxygen;
        built.Y.soil_gas.oxygen_phase=nitrogen.diagnostics.oxygen.phase_new;
        built.workspace.step.trial.nh4 = nitrogen.trial.nh4;
        built.workspace.step.trial.no3 = nitrogen.trial.no3;
        built.workspace.step.trial.urea = nitrogen.trial.urea;
        built.workspace.step.trial.C_NH4_granule = ...
            nitrogen.trial.C_NH4_granule;
        built.Y.nh4 = nitrogen.trial.nh4;
        built.Y.no3 = nitrogen.trial.no3;
        built.Y.urea = nitrogen.trial.urea;
        built.Y.C_NH4_granule = nitrogen.trial.C_NH4_granule;
        built.workspace.nitrogen = nitrogen.diagnostics;
        if isfield(nitrogen,'root_reservation')
            built.workspace.nitrogen.root_reservation = nitrogen.root_reservation;
        end
        if isfield(nitrogen,'root_remaining')
            built.workspace.nitrogen.root_remaining = nitrogen.root_remaining;
        end
        if isfield(nitrogen,'root_actual')
            built.workspace.nitrogen.root_actual = nitrogen.root_actual;
        end
        if isfield(nitrogen,'organic') && ~isempty(nitrogen.organic)
            built.Y.soil_organic = nitrogen.organic;
            % Keep the DOM sole-owner provenance in continuous and resumed
            % trajectories identical; checkpoint normalization writes the
            % same contract on legacy files.
            built.Y.soil_organic.qon_ownership = 'dom_total_sole_owner';
            built.workspace.soil_organic = struct( ...
                'exchange',nitrogen.organic_exchange, ...
                'diagnostics',nitrogen.organic_diagnostics);
        end
        if isfield(nitrogen,'dom') && ~isempty(nitrogen.dom)
            built.Y.soil_dom=nitrogen.dom.state;
            built.workspace.soil_dom=nitrogen.dom;
        end
        if isfield(nitrogen,'microbes') && ~isempty(nitrogen.microbes)
            built.Y.soil_microbes=nitrogen.microbes.state;
            built.workspace.soil_microbes=nitrogen.microbes;
        end
        if isfield(nitrogen,'fertilizer_tracer') ...
                &&~isempty(nitrogen.fertilizer_tracer)
            tracer=nitrogen.fertilizer_tracer;
            assert(isfield(tracer,'state')&&isfield(tracer,'diagnostics'), ...
                'soil:commitCandidate:fertilizerTracer', ...
                'Accepted fertilizer tracer lacks state or diagnostics.');
            built.Y.soil_fertilizer_n_tracer=tracer.state;
            built.workspace.fertilizer_n_tracer=tracer.diagnostics;
        end
    end
end
