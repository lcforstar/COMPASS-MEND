function signature=checkpointCompatibilitySignature(built)
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- bind new evaporation policy, multiplier and VG shape to restart.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- bind resolved CSV initialization contract.
% Change: 2026-10-17 Wenjie Liu -- bind AOM-first routing, DOM outer mean,
%   pH switch/shape and active microbial/enzyme C:N to prevent mixed restarts.
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
% Change: 2026-10-12 Wenjie Liu -- bind layer-capable CUE, alpha, mortality
% multiplier and dormancy thresholds in addition to their derived rates.
% Change: 2026-10-23 Wenjie Liu -- bind the inorganic-N kinetic pH policy
%   so a response-enabled trajectory cannot silently resume with unit factors.
% Change: 2026-10-23 Wenjie Liu -- bind surface-layer volatilization pH and
%   weather-driven winter-wheat phenology, including observed-LAI PlantCN-off runs.
% Change: 2026-10-30 Wenjie Liu -- root-residue runs bind the complete
%   organic event ordering, paired masses and fixed Jackson allocation so
%   old applied flags cannot silently resume a changed schedule. No-root
%   configurations retain the previous compatibility fields.
% Change: 2026-11-01 Wenjie Liu -- bind nonzero NH4 protection coefficient
%   and the actual mapped bottom-up threshold profile; zero keeps old signatures.
% Change: 2026-11-03 Wenjie Liu -- bind active nonzero DIN ADE numerical
% options. Existing unsigned zero-option trajectories remain compatible.
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
    signature.surface_resistance_policy='field_capacity_8206_4255_v1';
    if isfield(built.parameters,'scope')
        signature.surface_rss_scale=built.parameters.scope.soil_optical.rss_scale;
    end
    signature.hydraulic_vg_n=built.parameters.soil.n;
    if readFlag(sw,'lNefc')
        signature.volatilization_pH_location='surface_layer_NL_v1';
        signature.nitrogen_ph_response=true;
        if isfield(built.parameters,'soil_n') ...
                &&isfield(built.parameters.soil_n,'ph_response_enabled')
            signature.nitrogen_ph_response=logical(built.parameters.soil_n.ph_response_enabled);
        end
    end
    if isfield(built.parameters,'soil_n')
        pn=built.parameters.soil_n;
        if isfield(pn,'nh4_threshold_enabled')&&pn.nh4_threshold_enabled
            assert(isfield(pn,'nh4_protected_per_clay_mg_kg') ...
                &&isnumeric(pn.nh4_protected_per_clay_mg_kg) ...
                &&isreal(pn.nh4_protected_per_clay_mg_kg) ...
                &&isscalar(pn.nh4_protected_per_clay_mg_kg) ...
                &&isfinite(pn.nh4_protected_per_clay_mg_kg) ...
                &&pn.nh4_protected_per_clay_mg_kg>0 ...
                &&isfield(pn,'S_min_nh4_mg_kg') ...
                &&isequal(size(pn.S_min_nh4_mg_kg),[built.domain.NL,1]) ...
                &&isreal(pn.S_min_nh4_mg_kg)&&all(isfinite(pn.S_min_nh4_mg_kg)) ...
                &&all(pn.S_min_nh4_mg_kg>=0), ...
                'checkpointCompatibilitySignature:nh4Protection', ...
                'Active NH4 protection requires its positive scalar and mapped NL-by-1 profile.');
            signature.nh4_protection_threshold=struct( ...
                'per_clay_mg_kg',pn.nh4_protected_per_clay_mg_kg, ...
                'layer_S_min_mg_kg',pn.S_min_nh4_mg_kg);
            % The original shared scope keeps its historical signature. A
            % different active scope must not resume that trajectory silently.
            if isfield(pn,'nh4_protection_scope')&&strcmp(pn.nh4_protection_scope,'nitrification_only')
                signature.nh4_protection_threshold.scope='nitrification_only';
            end
        end
    end
    if isfield(built.parameters,'veg') ...
            &&vegetation.getPFTValue(built.parameters.veg.structure.PFT_type, ...
                'winter_wheat_cn',0)==1
        signature.winter_wheat_phenology_driver='daily_weather_with_observed_LAI_v1';
    end
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
    if readFlag(sw,'lNefc')
        dinOptions=soil.nitrogen.transport.numericalOptions(built.control);
        if dinOptions.lUpW||dinOptions.lArtD
            signature.din_transport_numerics=dinOptions;
            if dinOptions.lUpW
                signature.din_transport_numerics.upwind_weight='signed_coth_small_pe_series_v1';
            end
        end
    end
    if signature.soil_dom
        % The free DOM trajectory depends on the liquid-water time support.
        % Do not resume old mean/end-theta checkpoints as ordinary v2 runs.
        signature.dom_water_driver='accepted_water_tiles_substep_start_theta_v2';
        dz=built.domain.dz(:);
        if any(diff(dz)~=0)
            signature.dom_internal_face_flux=struct( ...
                'operator','linear_at_physical_face_v1','layer_thickness_cm',dz);
        end
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
    % Change: 2026-10-30 Wenjie Liu -- never infer historical manure labels
    % from unlabelled checkpoints. Disabled legacy signatures stay unchanged.
    if readFlag(sw,'lManureSourceTracer')
        assert(isfield(built.Y,'soil_manure_n_tracer'), ...
            'checkpointCompatibilitySignature:missingManureTracer', ...
            'Enabled manure tracer has no accepted state.');
        tags=built.Y.soil_manure_n_tracer;
        signature.manure_source_tracer=true;
        signature.manure_tracer_sources=struct( ...
            'schema_version',tags.schema_version, ...
            'organic_schema',tags.organic.schema_version, ...
            'tan_schema',tags.tan.schema_version, ...
            'organic_definition',tags.organic.source_definition, ...
            'tan_definition',tags.tan.source_definition);
    end
    if isfield(built,'sources')&&isfield(built.sources,'organic_fertilization_schedule') ...
            &&any(strcmp(built.sources.organic_fertilization_schedule.material_type,'root_residue'))
        schedule=built.sources.organic_fertilization_schedule;
        schedule.time=cellstr(string(schedule.time,'yyyy-MM-dd''T''HH:mm:ss'));
        rai=built.parameters.veg.root.RAI(:);
        assert(all(isfinite(rai))&&all(rai>=0)&&sum(rai)>0, ...
            'checkpointCompatibilitySignature:rootResidueWeights', ...
            'Root-residue compatibility requires a finite positive Jackson profile.');
        signature.root_residue_management=struct('schema_version',1, ...
            'organic_schedule',schedule,'allocation','configured_jackson_RAI', ...
            'weights_bottom_up',rai./sum(rai), ...
            'root_fractions',built.parameters.soil_organic.aom.litter_root_fractions(:), ...
            'input_owner','external_root_residue_into_component_litter_AOM');
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
            'Vg_ref_h','Vm_ref_h','KsDOM_ref_mg_g','rMORT_ref_h','VmD_ref_h', ...
            'Yg_ref','alpha','gamma','SWP_A2D_MPa','tau'};
        kinetics=struct('schema_version',2,'grid_order','top_down');
        for j=1:numel(keys)
            key=keys{j};value=p.(key);
            if size(value,1)==1,value=repmat(value,n,1);end
            kinetics.(key)=value;
        end
        signature.mend_layer_kinetics=kinetics;
        signature.mend_ph_response=struct('enabled',p.ph_response_enabled, ...
            'sensitivity_scale',p.ph_sensitivity_scale, ...
            'optimum_adapt_per_ph',p.ph_optimum_adapt_per_ph);
        signature.mend_stoichiometry=struct('CN_MB_min',p.CN_MB_min, ...
            'CN_MB_max',p.CN_MB_max,'CN_ENZP',p.CN_ENZP,'CN_ENZM',p.CN_ENZM);
        signature.aom_particulate_n='pom2_target_then_lignin_remainder_low_n_carbon_shares_v1';
        if isfield(p,'manure_fast_input_enabled')&&p.manure_fast_input_enabled
            signature.external_organic_inputs.manure_fast_input_enabled=true;
            signature.aom_particulate_n='material_manure_partition_fast_dom_remaining_AOM_v1';
        end
    end
    if isfield(built.parameters,'soil_organic')&&isfield(built.parameters.soil_organic,'initial_dom')
        signature.organic_initialization='resolved_soc_fraction_csv_v1';
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
