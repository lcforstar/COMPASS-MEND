function [values,owner,owners]=resolveMendLayerParameters(soil,microbe)
%RESOLVEMENDLAYERPARAMETERS Resolve layer-capable MEND reaction candidates.
% SOIL is a raw table (rows retain public depth order) or the mapped soil
% struct (NL-by-1, bottom-up). MICROBE is config.microbe. Returned VALUES
% preserve soil order; legacy values remain scalars. No unit conversion or
% depth interpolation happens here. OWNER is 'soil' or 'legacy_site'.
% The historical four-column group retains sole atomic ownership. Each new
% physiological, dormancy, DIN uptake or sorption parameter has one soil-column
% or legacy site-scalar owner. OWNERS records each effective owner.
% Editor: Wenjie Liu; Date: 2026-09-29.
% Change: 2026-10-12 Wenjie Liu — resolve the 8+2 candidate interface and
% expose per-parameter ownership without changing the legacy four-column rule.
% Change: 2026-10-12 Wenjie Liu — allow soil-only ownership for all ten
% candidates and reject duplicate extended-parameter definitions.
% Change: 2026-10-27 Wenjie Liu — validate layered dormancy and DIN uptake domains.
% Change: 2026-10-31 Wenjie Liu -- resolve kinetic sorption by sole owner;
% require it in source_hourly and retain optional absence for legacy RK4 sites.
    spec=io.adapters.mendLayerParameterSpec();base=spec.base_atomic_names;
    assert((istable(soil)||isstruct(soil))&&isstruct(microbe), ...
        'resolveMendLayerParameters:container','Expected soil table/struct and microbe struct.');
    if istable(soil),columns=soil.Properties.VariableNames;else,columns=fieldnames(soil);end
    baseInSoil=ismember(base,columns);baseInSite=isfield(microbe,base);
    if any(baseInSoil)
        assert(~any(baseInSite),'resolveMendLayerParameters:duplicateOwner', ...
            'MEND layer kinetics cannot be defined in both soil.csv and site.toml [microbe].');
        assert(all(baseInSoil),'resolveMendLayerParameters:partialProfile', ...
            'soil.csv must supply the complete four-column MEND kinetic group.');
        owner='soil';
    else
        assert(all(baseInSite),'resolveMendLayerParameters:missing', ...
            'Supply all four MEND kinetics in soil.csv, or the complete legacy site scalar group.');
        owner='legacy_site';
    end
    sourceHourly=isfield(microbe,'reaction_method') ...
        &&strcmp(string(microbe.reaction_method),'source_hourly');
    values=struct();owners=struct();profileCount=[];
    for i=1:numel(spec.names)
        key=spec.names{i};
        if ismember(key,base)
            useSoil=strcmp(owner,'soil');
        else
            useSoil=ismember(key,columns);useSite=isfield(microbe,key);
            if ismember(key,spec.sorption_names)&&~sourceHourly&&~useSoil&&~useSite
                continue;
            end
            assert(~(useSoil&&useSite),'resolveMendLayerParameters:duplicateOwner', ...
                'MEND parameter %s cannot be defined in both soil.csv and site.toml [microbe].',key);
            assert(useSoil||useSite,'resolveMendLayerParameters:missing', ...
                'Supply MEND parameter %s in soil.csv or as a legacy site scalar.',key);
        end
        if useSoil,source=soil;owners.(key)='soil';else,source=microbe;owners.(key)='legacy_site';end
        v=source.(key);
        assert(isnumeric(v)&&isreal(v)&&iscolumn(v)&&~isempty(v) ...
            &&all(isfinite(v)), ...
            'resolveMendLayerParameters:value', ...
            '%s:%s must be a finite real scalar/column.',owners.(key),key);
        validateValue(key,v);
        if ~useSoil
            assert(isscalar(v),'resolveMendLayerParameters:legacyScalar', ...
                'site.toml microbe.%s is legacy scalar-only; use soil.csv for a profile.',key);
        elseif isempty(profileCount),profileCount=numel(v);
        else
            assert(numel(v)==profileCount,'resolveMendLayerParameters:shape', ...
                'MEND soil kinetic columns must have matching lengths.');
        end
        values.(key)=double(v);
    end
end

function validateValue(key,value)
    switch key
        case {'mend_alpha'}
            ok=value>0&value<0.5;
        case {'mend_cue_20'}
            ok=value>0&value<1;
        case {'mend_tau'}
            ok=value>0&value<=1;
        case {'mend_mortality_multiplier','mend_vnup_mb_h', ...
                'mend_kba_g_mg','mend_kdes_mg_g_h'}
            ok=value>=0;
        case {'mend_wdorm'}
            ok=value>=1;
        otherwise
            ok=value>0;
    end
    assert(all(ok),'resolveMendLayerParameters:value', ...
        'MEND parameter %s is outside its physical input domain.',key);
end
