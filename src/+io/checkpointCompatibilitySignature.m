function signature=checkpointCompatibilitySignature(built)
% Change: 2026-09-30 Wenjie Liu -- selectively restore COMPASS canopy
% behavior; retain winter-wheat phenology and CoLM C/N allocation.
%CHECKPOINTCOMPATIBILITYSIGNATURE Describe restart-critical BGC topology.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-08-28 Codex -- include native organic-C/N, dormancy, and
%     BIOPARTICLE topology in the restart-compatibility contract.
% Change: 2026-09-14 Jizixing Zhu -- bind source-hour restarts to their explicit
%   hourly policy/parameters in addition to the distinct pool schemas.
% Change: 2026-09-16 Jizixing Zhu -- bind source-hour equilibrium restarts
% to the DOC/DON Kd profiles that define their algebraic phase inventories.
% Change: 2026-09-21 Wenjie Liu -- bind restarts to the fertilizer-source
% tracer schema and immutable source/cohort definition.
% Change: 2026-09-29 Wenjie Liu -- bind restarts to canonical MEND kinetic
% profiles. Pre-migration checkpoints lack this evidence and are rejected.
% Change: 2026-09-30 Wenjie Liu -- winter-wheat physical NSC funding cannot
% resume a legacy debt/writeoff trajectory, even when endpoint xsmr is zero.
% Change: 2026-09-30 Wenjie Liu -- separately bind annual reserve exclusion,
% grain-N lifecycle and harvest cleanup; a debt-free older crop checkpoint
% does not establish compatibility with the revised annual crop contracts.
% Change: 2026-10-10 Wenjie Liu -- bind the once-only plant/manure entry
% route, external C:N closure and conservative management placement.
    sw=built.control.switches;
    signature=struct('NL',built.domain.NL,'NN',built.domain.NN, ...
        'soil_organic',readFlag(sw,'lSoilOrganic'), ...
        'soil_dom',readFlag(sw,'lSoilDOM'), ...
        'explicit_soil_microbes',readFlag(sw,'lExplicitSoilMicrobes'), ...
        'litter_input',readFlag(sw,'lLitterInput'), ...
        'microbial_dormancy',readFlag(sw,'lMicrobialDormancy'), ...
        'microbial_transport',readFlag(sw,'lMicrobialTransport'), ...
        'fertilizer_source_tracer',readFlag(sw,'lFertilizerSourceTracer'), ...
        'organic_schema',stateSchema(built.Y,'soil_organic'), ...
        'dom_schema',stateSchema(built.Y,'soil_dom'), ...
        'microbe_schema',stateSchema(built.Y,'soil_microbes'));
    if readFlag(sw,'lPlantCN')&&readFlag(sw,'lNefc') ...
            &&vegetation.getPFTValue(built.parameters.veg.structure.PFT_type, ...
                'winter_wheat_cn',0)==1
        signature.winter_wheat_carbon_funding='physical_nsc_v1';
        signature.winter_wheat_crop_lifecycle='observed_harvest_maturity_separate_v2';
        signature.winter_wheat_crop_scope='colm_allocation_compass_litter_v1';
        % Old checkpoints predate the layer-weighted actual soil-N inflow
        % and candidate replay. Pool schema alone cannot prove compatibility.
        signature.winter_wheat_actual_root_n='accepted_input_replay_layer_weighted_v1';
        signature.winter_wheat_parameters=canonicalPlantParameters(built.parameters.plant_cn.resolved);
        signature.winter_wheat_harvest_rule=built.parameters.veg.phenology.crop_harvest_rule;
        signature.ecosystem_carbon_exchange='accepted_owner_cumulative_v1';
    end
    % Change: 2026-09-29 Wenjie Liu -- reject a pre-fix initial trajectory
    % when resuming the new measured-water/DIN initialization protocol.
    if isfield(built.p,'initialization') ...
        &&isfield(built.p.initialization,'measured_soil') ...
        &&built.p.initialization.measured_soil.applied
        signature.measured_soil_initialization=1;
    end
    % Change: 2026-09-29 Wenjie Liu -- bind DOM boundary and export-ledger
    % contracts; old unsigned/closed checkpoints cannot resume an open run.
    if signature.soil_dom
        if isfield(built.parameters.soil_dom,'boundary')
            signature.dom_boundary=built.parameters.soil_dom.boundary;
        else
            signature.dom_boundary=soil.dom.normalizeBoundary(struct());
        end
    end
    % Change: 2026-09-29 Wenjie Liu -- bind cell-particle export/clock contract.
    if signature.microbial_transport
        p=built.parameters.soil_microbes.bioparticle;
        if isfield(p,'boundary'),signature.bioparticle_boundary=p.boundary;
        else,signature.bioparticle_boundary=soil.microbe.normalizeTransportBoundary(struct());end
    end
    if signature.fertilizer_source_tracer
        assert(isfield(built.Y,'soil_fertilizer_n_tracer'), ...
            'checkpointCompatibilitySignature:missingFertilizerTracer', ...
            'Enabled fertilizer source tracer has no accepted state.');
        tracer=built.Y.soil_fertilizer_n_tracer;
        signature.fertilizer_tracer_schema=tracer.schema_version;
        signature.fertilizer_tracer_source_definition=tracer.source_definition;
    end
    if signature.explicit_soil_microbes
        p=built.parameters.soil_microbes;n=built.domain.NL;
        signature.external_organic_inputs=struct( ...
            'route',p.external_input_route, ...
            'placement',p.management_placement, ...
            'plant_carbon_fraction',p.input_fraction, ...
            'manure_carbon_fraction',p.manure_input_fraction, ...
            'input_cn_pom2',p.CN_POM2, ...
            'input_cn_lignin_to_labile',p.input_cn_lignin_to_labile);
        keys={'VdPOM_ref_h','KsPOM_ref_mg_g','VdMOM_ref_h','KsMOM_ref_mg_g', ...
            'Vg_ref_h','Vm_ref_h','KsDOM_ref_mg_g','rMORT_ref_h','VmD_ref_h'};
        kinetics=struct('schema_version',1,'grid_order','top_down');
        for j=1:numel(keys)
            key=keys{j};value=p.(key);
            if size(value,1)==1,value=repmat(value,n,1);end
            kinetics.(key)=value;
        end
        signature.mend_layer_kinetics=kinetics;
    end
    if isfield(built.parameters,'soil_microbes') ...
        &&isfield(built.parameters.soil_microbes,'hourly')
        signature.mend_hourly=built.parameters.soil_microbes.hourly;
        if strcmp(built.parameters.soil_dom.sorption_model,'linear_equilibrium')
            signature.dom_equilibrium=struct('sorption_model','linear_equilibrium', ...
                'kd_doc_cm3_g',built.parameters.soil_dom.kd_doc_cm3_g, ...
                'kd_don_cm3_g',built.parameters.soil_dom.kd_don_cm3_g);
        end
    end
end

function params=canonicalPlantParameters(params)
% Inactive PFT entries deliberately use NaN. Encode their declared absence
% so the existing strict isequal restart guard can compare identical views;
% retain every finite active value and reject Inf rather than weakening it.
% Change: 2026-09-30 Wenjie Liu -- bind compiled inactive parameter sentinels.
    keys=fieldnames(params);
    for j=1:numel(keys)
        key=keys{j};value=params.(key);
        if isnumeric(value)
            assert(isscalar(value)&&isreal(value)&&~isinf(value), ...
                'checkpointCompatibilitySignature:plantParameter', ...
                'Compiled plant parameter %s must be a real scalar.',key);
            if isnan(value),params.(key)='inactive_pft_nan';end
        end
    end
end

function value=readFlag(s,name)
    value=false; if isfield(s,name), value=logical(s.(name)); end
end

function value=stateSchema(Y,name)
    value=0;
    if isfield(Y,name)&&isfield(Y.(name),'schema_version')
        value=Y.(name).schema_version;
    end
end
