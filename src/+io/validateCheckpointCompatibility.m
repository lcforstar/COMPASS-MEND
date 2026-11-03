function validateCheckpointCompatibility(built,ckpt)
%VALIDATECHECKPOINTCOMPATIBILITY Reject incompatible BGC restart topology.
%
%   Schema-1 checkpoints are accepted only when DOM, explicit microbes, and
%   litter are all disabled. Schema-2 signatures must match exactly.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-09-30 Wenjie Liu -- schema-1 cannot bypass the new
%   winter-wheat physical-NSC funding restart contract.
% Change: 2026-10-30 Wenjie Liu -- reject removal of a root-residue schedule
%   from an existing root-enabled restart; absent old no-root fields remain valid.
% Change: 2026-11-01 Wenjie Liu -- nonzero NH4 protection requires a matching
%   schema-2 signature and cannot be silently removed from a protected restart.
% Change: 2026-11-03 Wenjie Liu -- reject changes to active DIN transport
% options, including switching from an explicitly routed scheme to zero.
    expected=io.checkpointCompatibilitySignature(built);
    if ckpt.schema_version==1
        if expected.soil_dom||expected.explicit_soil_microbes||expected.litter_input ...
                ||expected.fertilizer_source_tracer ...
                ||isfield(expected,'manure_source_tracer') ...
                ||isfield(expected,'nh4_protection_threshold') ...
                ||isfield(expected,'din_transport_numerics') ...
                ||isfield(expected,'winter_wheat_carbon_funding')
            error('validateCheckpointCompatibility:legacyBgcCheckpoint', ...
                'Schema-1 checkpoint lacks the signature required by the active model processes or DIN transport options.');
        end
        return;
    end
    assert(isfield(ckpt,'compatibility')&&isstruct(ckpt.compatibility), ...
        'validateCheckpointCompatibility:missingSignature', ...
        'Schema-2 checkpoint lacks a compatibility signature.');
    actual=ckpt.compatibility;
    if isfield(actual,'dom_internal_face_flux')&&~isfield(expected,'dom_internal_face_flux')
        error('validateCheckpointCompatibility:mismatch', ...
            'The signed DOM face geometry does not match the active configuration.');
    end
    if isfield(actual,'din_transport_numerics')&&~isfield(expected,'din_transport_numerics')
        error('validateCheckpointCompatibility:mismatch', ...
            'Checkpoint DIN transport numerical options cannot be removed from the active site.');
    end
    if isfield(actual,'manure_source_tracer')&&~isfield(expected,'manure_source_tracer')
        error('validateCheckpointCompatibility:mismatch', ...
            'Checkpoint manure-source labels cannot be silently removed.');
    end
    if isfield(actual,'root_residue_management')&&~isfield(expected,'root_residue_management')
        error('validateCheckpointCompatibility:mismatch', ...
            'Checkpoint root-residue management cannot be removed from the active site.');
    end
    if isfield(actual,'nh4_protection_threshold')&&~isfield(expected,'nh4_protection_threshold')
        error('validateCheckpointCompatibility:mismatch', ...
            'Checkpoint NH4 protection cannot be removed from the active site.');
    end
    names=fieldnames(expected);
    for i=1:numel(names)
        name=names{i};
        if strcmp(name,'fertilizer_source_tracer') ...
                &&~expected.fertilizer_source_tracer ...
                &&~isfield(actual,name)
            % Compatibility for schema-2 checkpoints written before the
            % optional tracer existed.  Enabled tracers still fail fast.
            continue
        end
        if ~isfield(actual,name)||~isequal(actual.(name),expected.(name))
            error('validateCheckpointCompatibility:mismatch', ...
                'Checkpoint compatibility field %s does not match the active site.',name);
        end
    end
end
