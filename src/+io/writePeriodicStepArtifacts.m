function artifacts = writePeriodicStepArtifacts(model, kt)
%WRITEPERIODICSTEPARTIFACTS  Write gated checkpoint and partial-step files.
%
%   Inputs:
%     model [in] Accepted model shell after diagnostics are written.
%     kt    [-] One-based accepted outer-step index.
%
%   Return:
%     artifacts [struct] Checkpoint path and partial-output status.
%
%   This function owns disk I/O only. Runtime clock/calendar and report state
%   are committed by the step package before and after this callback.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    checkpoint = model.runtime.checkpoint;
    checkpoint_file = "";
    if ~checkpoint.in_optimizer && checkpoint.checkpoint_interval > 0 ...
            && mod(kt, checkpoint.checkpoint_interval) == 0
        checkpoint_file = string(io.saveCompassCheckpoint(model, kt));
    end

    partial_written = false;
    if checkpoint.write_periodic_partial && ~checkpoint.in_optimizer ...
            && kt > 1 && checkpoint.partial_output_interval > 0 ...
            && mod(kt, checkpoint.partial_output_interval) == 0
        registry = model.runtime.diagnostics.registry;
        config = init.makeLegacyRuntimeCfg(model.built);
        config.OutputPath = checkpoint.output_dir;
        config.Time = model.built.control.time;
        io.saveCompassPeriodicPartialOutput( ...
            registry, config, '', kt, kt + 1);
        partial_written = true;
    end
    artifacts = struct( ...
        'checkpoint_file', checkpoint_file, ...
        'partial_output_written', partial_written);
end

