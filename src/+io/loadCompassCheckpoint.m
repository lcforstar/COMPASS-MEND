function ckpt = loadCompassCheckpoint(file)
%LOADCOMPASSCHECKPOINT  Read a resume checkpoint written by io.saveCompassCheckpoint.
%
%   ckpt = io.loadCompassCheckpoint(file) loads and schema-validates the checkpoint
%   struct. It does NOT apply it to a model — beginRun (layer 3) reconstructs the
%   live state (Y, workspace via rebuildSoilAuxiliaryFromAccepted, LUNA buffer via
%   Climate10dBuffer.fromStruct, runtime memory, next_kt) from this struct.
%
%   Inputs:
%     file [char]  path to a checkpoint .mat (single `ckpt` variable).
%
%   Returns:
%     ckpt [struct]  the validated checkpoint (schema_version 1 or 2).
%
%   Editor: Jizixing Zhu; Date: 2026-07-01.
%   Change: 2026-08-21 Jizixing Zhu -- accept schema 2 compatibility records
%     while retaining schema 1 for strict legacy fallback resumes.
%   Change: 2026-09-12 Jizixing Zhu -- audit DON/QON ownership before resume;
%     reject ambiguous nonzero legacy QON instead of deleting or double counting.

    if ~isfile(file)
        error('loadCompassCheckpoint:notFound', 'checkpoint file not found: %s', file);
    end
    s = load(file, 'ckpt');
    if ~isfield(s, 'ckpt') || ~isstruct(s.ckpt)
        error('loadCompassCheckpoint:badFile', 'file does not contain a ckpt struct: %s', file);
    end
    ckpt = s.ckpt;
    if ~isfield(ckpt, 'schema_version') ...
            || ~ismember(ckpt.schema_version,[1,2])
        error('loadCompassCheckpoint:badSchema', ...
            'unsupported checkpoint schema_version (expected 1 or 2).');
    end
    if ~isfield(ckpt, 'Y') || ~isfield(ckpt, 'runtime') || ~isfield(ckpt, 'kt')
        error('loadCompassCheckpoint:badSchema', 'checkpoint missing Y/runtime/kt.');
    end
    [ckpt.Y,ckpt.don_ownership_audit]=io.normalizeCheckpointDON(ckpt.Y);
end
