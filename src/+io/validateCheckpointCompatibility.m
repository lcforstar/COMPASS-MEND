function validateCheckpointCompatibility(built,ckpt)
%VALIDATECHECKPOINTCOMPATIBILITY Reject incompatible BGC restart topology.
%
%   Schema-1 checkpoints are accepted only when DOM, explicit microbes, and
%   litter are all disabled. Schema-2 signatures must match exactly.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
    expected=io.checkpointCompatibilitySignature(built);
    if ckpt.schema_version==1
        if expected.soil_dom||expected.explicit_soil_microbes||expected.litter_input
            error('validateCheckpointCompatibility:legacyBgcCheckpoint', ...
                'Schema-1 checkpoints cannot initialize active DOM/MEMC/litter states.');
        end
        return;
    end
    assert(isfield(ckpt,'compatibility')&&isstruct(ckpt.compatibility), ...
        'validateCheckpointCompatibility:missingSignature', ...
        'Schema-2 checkpoint lacks a compatibility signature.');
    actual=ckpt.compatibility;
    names=fieldnames(expected);
    for i=1:numel(names)
        name=names{i};
        if ~isfield(actual,name)||~isequal(actual.(name),expected.(name))
            error('validateCheckpointCompatibility:mismatch', ...
                'Checkpoint compatibility field %s does not match the active site.',name);
        end
    end
end
