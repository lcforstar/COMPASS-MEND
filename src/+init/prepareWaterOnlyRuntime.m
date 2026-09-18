function model = prepareWaterOnlyRuntime(model)
%PREPAREWATERONLYRUNTIME  Allocate runtime state for a laboratory water column.
%
%   model = init.prepareWaterOnlyRuntime(model) creates only the clock,
%   adaptive-substep state, and run-control counters consumed by
%   solver.runWaterOnlyExperimentStep. It deliberately does not register field
%   diagnostics or allocate output, checkpoint, groundwater, phenology, canopy,
%   SCOPE, or management services.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.

    built = model.built;
    dt = built.control.time.Delt_t;

    runtime = struct();
    runtime.substep = struct( ...
        'Delt_t', dt, 'Delt_old', dt, 'retry_count', 0, ...
        'success_count', 0, 'hold_count', 0, 'last_iters', 0);
    runtime.step = struct( ...
        'TIME', 0, 'clock_time', built.control.time.start_time, ...
        'current_day', day(built.control.time.start_time), ...
        'current_hour', hour(built.control.time.start_time));
    runtime.run_control = struct( ...
        'n_steps', built.control.time.tS, 'max_kt', NaN, ...
        'next_kt', 1, 'last_completed_kt', 0);
    runtime.dimensions = struct( ...
        'tS', built.control.time.tS, 'NL', built.domain.NL, ...
        'NN', built.domain.NN, 'nCell', 1);
    runtime.status = "prepared_water_only";
    runtime.created_at = datetime("now");

    model.runtime = runtime;
    model.status = "runtime_prepared";
end
