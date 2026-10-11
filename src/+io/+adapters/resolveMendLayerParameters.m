function [values,owner]=resolveMendLayerParameters(soil,microbe)
%RESOLVEMENDLAYERPARAMETERS Validate the sole owner of four kinetic inputs.
% SOIL is a raw table (rows retain public depth order) or the mapped soil
% struct (NL-by-1, bottom-up). MICROBE is config.microbe. Returned VALUES
% preserve soil order; legacy values remain scalars. No unit conversion or
% depth interpolation happens here. OWNER is 'soil' or 'legacy_site'.
% Rejects partial groups, dual ownership, nonpositive/nonfinite values and
% invalid shapes, rather than silently falling back to another source.
% Editor: Wenjie Liu; Date: 2026-09-29.
    spec=io.adapters.mendLayerParameterSpec();keys=spec.names;
    assert((istable(soil)||isstruct(soil))&&isstruct(microbe), ...
        'resolveMendLayerParameters:container','Expected soil table/struct and microbe struct.');
    if istable(soil),columns=soil.Properties.VariableNames;else,columns=fieldnames(soil);end
    inSoil=ismember(keys,columns);inSite=isfield(microbe,keys);
    if any(inSoil)
        assert(~any(inSite),'resolveMendLayerParameters:duplicateOwner', ...
            'MEND layer kinetics cannot be defined in both soil.csv and site.toml [microbe].');
        assert(all(inSoil),'resolveMendLayerParameters:partialProfile', ...
            'soil.csv must supply the complete four-column MEND kinetic group.');
        source=soil;owner='soil';
    else
        assert(all(inSite),'resolveMendLayerParameters:missing', ...
            'Supply all four MEND kinetics in soil.csv, or the complete legacy site scalar group.');
        source=microbe;owner='legacy_site';
    end
    values=struct();count=[];
    for i=1:numel(keys)
        key=keys{i};v=source.(key);
        assert(isnumeric(v)&&isreal(v)&&iscolumn(v)&&~isempty(v) ...
            &&all(isfinite(v))&&all(v>0), ...
            'resolveMendLayerParameters:value', ...
            '%s:%s must be a finite positive scalar/column with its declared units.',owner,key);
        if strcmp(owner,'legacy_site')
            assert(isscalar(v),'resolveMendLayerParameters:legacyScalar', ...
                'site.toml microbe.%s is legacy scalar-only; use soil.csv for a profile.',key);
        elseif isempty(count),count=numel(v);
        else
            assert(numel(v)==count,'resolveMendLayerParameters:shape', ...
                'MEND soil kinetic columns must have matching lengths.');
        end
        values.(key)=double(v);
    end
end
