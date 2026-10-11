function identity = buildBiogeochemistryCycleIdentity(model)
%BUILDBIOGEOCHEMISTRYCYCLEIDENTITY Bind a repeated-forcing state to its model.
% Input model includes raw inputs and the built production container.
% Only immutable domain, parameters,
% forcing, management, control and vegetation drivers are hashed; accepted Y,
% disposable caches, output paths and wall-clock timestamps are excluded.
% The source hash covers the current src .m files. Equality is required, not
% approximate parameter similarity. This is a cycle/restart contract, not a
% claim that a measured-constrained production initial state is equilibrated.
% Editor: Wenjie Liu; Date: 2026-09-29.
    assert(isfield(model,'built')&&isfield(model,'raw'), ...
        'buildBiogeochemistryCycleIdentity:model','Full model and raw inputs are required.');
    built=model.built;
    identity = struct('schema_version',2);
    groups = {'domain','control','parameters','boundary','sources','forcing', ...
        'vegetation_state'};
    for i = 1:numel(groups)
        name = groups{i};
        assert(isfield(built,name), ...
            'buildBiogeochemistryCycleIdentity:missingGroup', ...
            'Cycle identity requires built.%s.',name);
        identity.([name '_sha256']) = digest( ...
            unicode2native(jsonencode(built.(name),ConvertInfAndNaN=false),'UTF-8'));
    end
    raw=model.raw;
    if isfield(raw,'siteDir'),raw=rmfield(raw,'siteDir');end
    identity.raw_sha256=digest(unicode2native( ...
        jsonencode(raw,ConvertInfAndNaN=false),'UTF-8'));
    identity.loaded_model_library_sha256=digest(unicode2native( ...
        jsonencode(built.p.scope_static,ConvertInfAndNaN=false),'UTF-8'));
    sourceRoot = fileparts(fileparts(mfilename('fullpath')));
    listing = dir(fullfile(sourceRoot,'**','*.m'));
    paths = arrayfun(@(f)fullfile(f.folder,f.name),listing,'UniformOutput',false);
    paths = sort(paths);
    entries = cell(numel(paths),1);
    for i = 1:numel(paths)
        file = paths{i}; fid = fopen(file,'rb');
        assert(fid>=0,'buildBiogeochemistryCycleIdentity:sourceRead', ...
            'Cannot read source file %s.',file);
        cleanup = onCleanup(@()fclose(fid));
        bytes = fread(fid,Inf,'*uint8');
        clear cleanup
        relative = strrep(file(numel(sourceRoot)+2:end),'\','/');
        entries{i} = [relative ':' digest(bytes)];
    end
    identity.source_sha256 = digest(unicode2native(strjoin(entries,newline),'UTF-8'));
    identity.signature = digest(unicode2native(jsonencode(identity),'UTF-8'));
end

function value = digest(bytes)
    engine = java.security.MessageDigest.getInstance('SHA-256');
    engine.update(uint8(bytes(:)));
    value = upper(reshape(dec2hex(typecast(engine.digest(),'uint8'),2).',1,[]));
end
