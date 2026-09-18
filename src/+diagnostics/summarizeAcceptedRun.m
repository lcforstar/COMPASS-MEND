function summary = summarizeAcceptedRun(registry, last_kt)
%SUMMARIZEACCEPTEDRUN  Integrate accepted water-ledger diagnostics at run end.
%
%   summary = diagnostics.summarizeAcceptedRun(registry, last_kt) returns
%   the cumulative surface--soil liquid-water budget over accepted rows
%   1:last_kt.  Each input is already an accepted-step amount [mm], so this
%   function only sums the discrete ledger; it neither recomputes fluxes nor
%   changes model state.
%
%   Inputs:
%     registry [in] diagnostics.Registry Accepted-step output registry.
%     last_kt  [in] [-] Last accepted outer-step index.
%
%   Returns:
%     summary [out] struct Water amounts [mm] and closure ratio [%].
%
%   Discrete closure:
%     residual = supply - evaporation - runoff - bottom_loss - root_sink
%                - storage_change.
%   This is exactly the outer-step surface--soil ledger written by
%   diagnostics.writeStep; NaN means an unavailable diagnostic and is not
%   silently treated as a physical zero.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.

    assert(isa(registry, 'diagnostics.Registry'), ...
        'diagnostics:summarizeAcceptedRun:badRegistry', ...
        'registry must be diagnostics.Registry.');
    assert(isscalar(last_kt) && isfinite(last_kt) && last_kt >= 1, ...
        'diagnostics:summarizeAcceptedRun:badLastKt', ...
        'last_kt must be a positive accepted index.');

    summary = struct('has_water_ledger', false, ...
        'liquid_supply_mm', NaN, 'evaporation_mm', NaN, ...
        'runoff_mm', NaN, 'bottom_loss_mm', NaN, 'root_sink_mm', NaN, ...
        'storage_change_mm', NaN, 'residual_mm', NaN, ...
        'closure_percent', NaN);
    names = registry.list();
    required = {'surface_soil_liquid_input', ...
        'surface_soil_evaporation_loss', 'surface_soil_runoff_loss', ...
        'surface_soil_bottom_loss', 'surface_soil_root_sink', ...
        'surface_soil_storage_change', 'surface_soil_water_balance_residual'};
    if ~all(ismember(required, names))
        return;
    end

    values = cell(size(required));
    for i = 1:numel(required)
        series = registry.get(required{i});
        values{i} = series(1:last_kt, 1);
        if any(~isfinite(values{i}))
            return;
        end
    end
    totals = cellfun(@sum, values);
    summary.has_water_ledger = true;
    summary.liquid_supply_mm = totals(1);
    summary.evaporation_mm = totals(2);
    summary.runoff_mm = totals(3);
    summary.bottom_loss_mm = totals(4);
    summary.root_sink_mm = totals(5);
    summary.storage_change_mm = totals(6);
    summary.residual_mm = totals(7);
    summary.closure_percent = 100 * abs(summary.residual_mm) / ...
        max(abs(summary.liquid_supply_mm), 1);
end
