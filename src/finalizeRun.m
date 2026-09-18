function model = finalizeRun(model)
%FINALIZERUN Flush accepted run diagnostics once after the outer-step loop.
%
% model = finalizeRun(model) writes one CF-style NetCDF file per registered
% diagnostic group. Only accepted rows 1:last_completed_kt are written, so a
% clean early stop never serializes uncomputed future rows. Optional legacy CSV
% output follows runtime.output.legacy_csv.
%
% This function performs output only. It does not update prognostic state,
% constitutive auxiliaries, boundary fluxes, or initialization values.
%
% Editor: Jizixing Zhu; Date: 2026-06-29.
% Change: 2026-07-16 Jizixing Zhu — implement the final Registry-to-NetCDF
% flush and optional legacy CSV output for the migrated run chain.
% Change: 2026-07-20 Jizixing Zhu -- migrate the annual water summary and opt-in three-way comparison hook.

validateFinalizeRun(model);

runtime = model.runtime;
output = runtime.output;
registry = runtime.diagnostics.registry;
lastKt = runtime.run_control.last_completed_kt;
if lastKt < 1
    warning('finalizeRun:noAcceptedSteps', ...
        'No accepted outer step is available; diagnostic files were not written.');
    model.runtime.finalization = struct( ...
        'flushed', false, 'last_completed_kt', 0, ...
        'completed_at', datetime('now'));
    return;
end

time = model.built.control.time;
switches = model.built.control.switches;
depthNode = flipud(model.built.domain.z_node(:));
runTag = sprintf('complete_run_cutoff_KT%05d', lastKt);
runMetadata = struct( ...
    'tag', runTag, ...
    'run_tag', runTag, ...
    'last_completed_index', lastKt, ...
    'start_time', char(time.start_time), ...
    'end_time', char(time.end_time), ...
    'delt_t_seconds', time.Delt_t, ...
    'depth_node', depthNode, ...
    'cfg_lFreezing', double(readSwitch(switches, 'lFreezing')), ...
    'cfg_lNefc', double(readSwitch(switches, 'lNefc')), ...
    'cfg_lSoilairefc', double(readSwitch(switches, 'lSoilairefc')));

diagnostics.flushToNetCDF(registry, output.output_dir, runMetadata, lastKt);
if output.legacy_csv
    diagnostics.flushToCsv(registry, output.output_dir);
end
runOptionalThreeWayComparison(output.output_dir, time.start_time);

% This is the migrated equivalent of the legacy end-of-loop water printout.
% It integrates only registry rows that were already accepted and written by
% diagnostics.writeStep; output must never become an alternate water ledger.
closure = diagnostics.summarizeAcceptedRun(registry, lastKt);
if closure.has_water_ledger
    fprintf(['[finalizeRun] Surface-soil water: supply=%+.2f, evap=%+.2f, ', ...
        'runoff=%+.2f, bottom=%+.2f, root=%+.2f, dS=%+.2f, ', ...
        'residual=%+.4f mm (%.3f%% of supply).\n'], ...
        closure.liquid_supply_mm, closure.evaporation_mm, closure.runoff_mm, ...
        closure.bottom_loss_mm, closure.root_sink_mm, ...
        closure.storage_change_mm, closure.residual_mm, closure.closure_percent);
end

model.runtime.finalization = struct( ...
    'flushed', true, ...
    'last_completed_kt', lastKt, ...
    'netcdf_directory', output.nc_dir, ...
    'legacy_csv_written', logical(output.legacy_csv), ...
    'surface_soil_water_closure', closure, ...
    'completed_at', datetime('now'));
fprintf('[finalizeRun] Wrote diagnostics through KT=%d to %s.\n', ...
    lastKt, output.output_dir);
end

function validateFinalizeRun(model)
% Validate the post-loop services required by the diagnostic flush.
assert(isstruct(model) && isfield(model, 'runtime'), ...
    'finalizeRun:missingRuntime', 'finalizeRun requires model.runtime.');
runtime = model.runtime;
assert(isfield(runtime, 'diagnostics') ...
    && isfield(runtime.diagnostics, 'registry'), ...
    'finalizeRun:missingRegistry', ...
    'finalizeRun requires runtime.diagnostics.registry.');
assert(isfield(runtime, 'output') && isfield(runtime.output, 'output_dir'), ...
    'finalizeRun:missingOutput', ...
    'finalizeRun requires runtime.output.output_dir.');
assert(isfield(runtime, 'run_control') ...
    && isfield(runtime.run_control, 'last_completed_kt'), ...
    'finalizeRun:missingRunControl', ...
    'finalizeRun requires runtime.run_control.last_completed_kt.');
end

function value = readSwitch(switches, name)
% Read one scalar physics switch; an absent optional switch is off.
value = isfield(switches, name) && logical(switches.(name));
end

function runOptionalThreeWayComparison(outputDir, startTime)
%RUNOPTIONALTHREEWAYCOMPARISON  Run the legacy opt-in post-run comparison.
%
%   This diagnostic is deliberately outside the physical time loop: it reads
%   only completed output files and cannot alter accepted state or ledgers.
%   io.getEnv preserves the old STEMMUS_THREEWAY_* environment names.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    if ~strcmpi(strtrim(io.getEnv('THREEWAY_COMPARE')), '1')
        return;
    end
    root = fileparts(fileparts(mfilename('fullpath')));
    comparisonDir = fullfile(root, 'tests', 'integration');
    try
        addpath(comparisonDir);
        days = str2double(io.getEnv('THREEWAY_DAYS'));
        if isfinite(days) && days > 0
            compare_soil_moisture_threeway(outputDir, days, startTime);
        else
            compare_soil_moisture_threeway(outputDir, [], startTime);
        end
    catch comparisonError
        warning('finalizeRun:threeWayComparisonFailed', ...
            'Three-way soil-moisture comparison failed: %s', ...
            comparisonError.message);
    end
end
