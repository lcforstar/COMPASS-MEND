function [resolved, sources] = compilePlantParameters(siteVegetation, pftType)
%COMPILEPLANTPARAMETERS Freeze numeric PFT values and explicit site overrides.
% Inputs are a read-only site vegetation struct and PFT name. Output scalars
% retain their PFT-documented units; NaN is retained only for inactive PFT
% fields. Site overrides must be finite numeric scalars and cannot silently
% fall back. This view owns parameters, never prognostic pools or diagnostics.
% Source: vegetation.getPFTParameterSet; no alternative physical defaults.
% Editor: Wenjie Liu; Date: 2026-09-30.
    table = vegetation.getPFTParameterSet();
    index = find(strcmpi({table.PFT}, pftType), 1);
    assert(~isempty(index), 'compilePlantParameters:unknownPFT', ...
        'Unknown PFT "%s".', char(pftType));
    row = table(index);
    resolved = struct(); sources = struct();
    names = fieldnames(row);
    for k = 1:numel(names)
        name = names{k}; value = row.(name);
        if ~isempty(value)&&(~isnumeric(value)||~isscalar(value)),continue;end
        source='PFT_table';
        if isfield(siteVegetation, name) && ~isempty(siteVegetation.(name))
            value = siteVegetation.(name);
            assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
                'compilePlantParameters:invalidOverride', ...
                'vegetation.%s must be a finite numeric scalar.', name);
            source = 'site.toml:vegetation';
        end
        if isempty(value) || ~isnumeric(value) || ~isscalar(value), continue; end
        sources.(name)=source;
        resolved.(name) = value;
    end
    % Optional legacy fields have named compatibility values. These match
    % pre-existing consumers and are inactive for annual winter wheat.
    defaults = struct('use_fruittree_alloc',0,'ft_GDDlfmat',0,'ft_bfact',0, ...
        'harvest_DS',0.95,'arootf2',0,'aleafstor',0,'grnrp',1100, ...
        'fleafi',0,'arooti',0,'arootf',0,'astemf',0,'grnfill',400, ...
        'hybgdd',2880,'froot_leaf',0,'croot_stem_alloc',0,'stem_leaf',0, ...
        'winter_wheat_cn',0,'resp_scheme',0,'resp_Km',2, ...
        'resp_jones_kf',0.227,'resp_jones_Yg',0.75, ...
        'xsmr_recovery_days',30,'clm_crop_respiration',0, ...
        'clm_br_s',2.525e-6,'clm_br_root_s',2.525e-6,'clm_resp_Q10',1.5, ...
        'root_turnover_rate',0.00274,'spring_release_days',60, ...
        'croot_stem',0,'nsc_frac',0.07,'f_storage_root',0.18, ...
        'labile_resp_days',1,'C_labile_abs_min',0.01,'initial_seed_c',3, ...
        'ndays_stor',30,'prune_fr',0,'r_mort',0, ...
        'fcur',0.95,'storage_cushion',1.2,'cn_flex_offset',10);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(resolved,name)
            resolved.(name) = defaults.(name);
            sources.(name) = 'explicit_legacy_compatibility';
            if isfield(siteVegetation,name) && ~isempty(siteVegetation.(name))
                value=siteVegetation.(name);
                assert(isnumeric(value)&&isscalar(value)&&isfinite(value), ...
                    'compilePlantParameters:invalidOverride', ...
                    'vegetation.%s must be a finite numeric scalar.',name);
                resolved.(name)=value;sources.(name)='site.toml:vegetation';
            end
        end
    end
end
