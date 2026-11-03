function accumulator = accumulateOxygenSubstep(accumulator, tile, dt_s)
%ACCUMULATEOXYGENSUBSTEP Add one accepted adaptive O2 tile exactly once.
% Editor: Jizixing Zhu; Date: 2026-08-21.
    assert(tile.converged && tile.admissible, ...
        'accumulateOxygenSubstep:unacceptedTile', ...
        'Only converged and physically admissible O2 tiles may accumulate.');
    if ~accumulator.has_data
        accumulator.storage_before_mol_cm2 = ...
            tile.ledger.storage_before_mol_cm2;
    end
    accumulator.has_data = true;
    accumulator.elapsed_s = accumulator.elapsed_s + dt_s;
    accumulator.accepted_tiles = accumulator.accepted_tiles + 1;
    accumulator.final_y_O2 = tile.y_O2;
    for name = {'Snit','Snit_N2O','f_aerobic','f_anaerobic','C_O2_liq'}
        field = [name{1} '_integral'];
        value = tile.(name{1});
        accumulator.(field) = accumulator.(field) + value * dt_s;
    end
    for name = {'QL','QL_h','QL_T','QL_a'}
        field = [name{1} '_integral'];
        value = tile.(name{1});
        accumulator.(field) = accumulator.(field) + value(:) * dt_s;
    end
    accumulator.QL_abs_integral = accumulator.QL_abs_integral ...
        + abs(tile.QL(:)) * dt_s;
    accumulator.storage_after_mol_cm2 = tile.ledger.storage_after_mol_cm2;
    for name = {'top_input_mol_cm2','bottom_input_mol_cm2', ...
            'reaction_loss_mol_cm2','positivity_correction_mol_cm2'}
        accumulator.(name{1}) = accumulator.(name{1}) + tile.ledger.(name{1});
    end
    accumulator.max_abs_balance_residual_mol_cm2 = max( ...
        accumulator.max_abs_balance_residual_mol_cm2, ...
        abs(tile.ledger.balance_residual_mol_cm2));
    accumulator.max_omitted_carrier_feedback_ratio = max( ...
        accumulator.max_omitted_carrier_feedback_ratio, ...
        tile.carrier_consistency.omitted_reaction_feedback_ratio);
    accumulator.minimum_background_dry_air_storage_g_cm3 = min( ...
        accumulator.minimum_background_dry_air_storage_g_cm3, ...
        tile.carrier_consistency.minimum_background_storage_g_cm3);
    accumulator.last_boundary_request = tile.boundary_request;
    accumulator.last_boundary_applied = tile.boundary_applied;
    accumulator.iterations_max = max(accumulator.iterations_max, tile.iterations);
end

% Change: 2026-10-23 Wenjie Liu -- restore the COMPASS-COMPASS_MATLAB
% oxygen transport/source boundary implementation; retain current carrier physics.
