function result=stagePlantInput(organicAccepted,microbeAccepted,domAccepted, ...
        plantAccepted,plantCandidate,built)
%STAGEPLANTINPUT Stage current-step PlantCN litter into the active soil owner.
%
%   Plant cumulative litter differences [g C/N m^-2] are converted to
%   [g cm^-2]. Leaf litter enters the top element; root litter follows the
%   configured beta depth distribution. The source-MEND type-1 partition sends
%   its DOM fraction directly to the sole DOC/DON owner. The complementary
%   fraction enters the native litter AOM front end, where winter-wheat PFT
%   material quality remains authoritative. Nothing is committed here.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-08-24 Jizixing Zhu -- include pruning and split background
%           wood mortality between surface and coarse-root litter pathways.
%   Change: 2026-08-28 Codex -- route all PlantCN returns to native litter
%           AOM C/N stores and remove Daisy/direct-POC ownership branches.
%   Change: 2026-09-08 Jizixing Zhu -- apportion new leaf/fine-root inputs
%   once into CoLM components; pruning and wood mortality enter CWD.
% Change: 2026-09-13 Jizixing Zhu -- apply the source-MEND type-1 DOM
%   fraction to accepted PlantCN litter while preserving AOM material quality.
    result=struct('organic',organicAccepted,'microbes',microbeAccepted, ...
        'dom',domAccepted,'carbon_added_g_cm2',0, ...
        'nitrogen_added_g_cm2',0,'aboveground_carbon_g_m2',0, ...
        'belowground_carbon_g_m2',0,'aboveground_nitrogen_g_m2',0, ...
        'belowground_nitrogen_g_m2',0, ...
        'direct_doc_carbon_g_cm2',0,'direct_don_nitrogen_g_cm2',0, ...
        'mend_dom_input_fraction',0,'triggered',false);
    if isempty(plantCandidate), return; end
    a=plantAccepted.accum;
    mortalityC=incrementOptional(plantCandidate,'cum_C_mortality',a);
    mortalityN=incrementOptional(plantCandidate,'cum_N_mortality',a);
    [mortalityAboveC,mortalityRootC]=splitWoodMortality( ...
        mortalityC,plantCandidate,'C');
    [mortalityAboveN,mortalityRootN]=splitWoodMortality( ...
        mortalityN,plantCandidate,'N');
    leafC=increment(plantCandidate,'cum_C_leaf_litter',a.cum_C_leaf_litter) ...
        +incrementOptional(plantCandidate,'cum_C_prune_litter',a) ...
        +mortalityAboveC;
    rootC=increment(plantCandidate,'cum_C_root_litter',a.cum_C_root_litter);
    rootC=rootC+mortalityRootC;
    leafN=increment(plantCandidate,'cum_N_leaf_litter',a.cum_N_leaf_litter) ...
        +incrementOptional(plantCandidate,'cum_N_prune_litter',a) ...
        +mortalityAboveN;
    rootN=increment(plantCandidate,'cum_N_root_litter',a.cum_N_root_litter);
    rootN=rootN+mortalityRootN;
    if leafC+rootC+leafN+rootN<=0, return; end
    assert(isfield(result.organic,'aom') ...
        &&isfield(result.organic.aom,'litter'), ...
        'stagePlantInput:missingNativeAom', ...
        'Plant litter requires the native COMPASS litter AOM store.');
    dzTop=flip(built.domain.dz(:));
    cProfile=zeros(size(dzTop)); nProfile=zeros(size(dzTop));
    w=rootWeights(dzTop,built.parameters.litter);
    cProfile=cProfile+(rootC/1e4).*w./dzTop;
    nProfile=nProfile+(rootN/1e4).*w./dzTop;
    assert(isfield(built.parameters,'soil_microbes') ...
        &&isfield(built.parameters.soil_microbes,'input_fraction'), ...
        'stagePlantInput:missingMendInputFraction', ...
        'Plant litter coupling requires the source-MEND input fraction.');
    inputFraction=built.parameters.soil_microbes.input_fraction(:);
    assert(numel(inputFraction)==3&&all(isfinite(inputFraction)) ...
        &&all(inputFraction>=0)&&abs(sum(inputFraction)-1)<=1e-12, ...
        'stagePlantInput:badMendInputFraction', ...
        'The source-MEND POM1/POM2/DOM input fractions must sum to one.');
    domFraction=inputFraction(3);
    aomFraction=1-domFraction;
    assert(isfield(result.dom,'doc_total_g_cm3') ...
        &&isfield(result.dom,'don_total_g_cm3'), ...
        'stagePlantInput:missingDomOwner', ...
        'A nonzero source-MEND DOM input requires the sole DOC/DON owner.');
    directDocTop=domFraction.*cProfile;
    directDonTop=domFraction.*nProfile;
    directDocTop(1)=directDocTop(1)+domFraction.*leafC./1e4./dzTop(1);
    directDonTop(1)=directDonTop(1)+domFraction.*leafN./1e4./dzTop(1);
    result.dom.doc_total_g_cm3=result.dom.doc_total_g_cm3+flip(directDocTop);
    result.dom.don_total_g_cm3=result.dom.don_total_g_cm3+flip(directDonTop);
    if isfield(result.organic.aom.litter,'components')
        p=built.parameters.soil_organic.aom;
        pool=soil.litter.componentTotals(result.organic.aom.litter,true);
        woodC=mortalityAboveC+incrementOptional(plantCandidate,'cum_C_prune_litter',a);
        woodN=mortalityAboveN+incrementOptional(plantCandidate,'cum_N_prune_litter',a);
        names={'metabolic','cellulose','lignin','cwd'};
        aboveC=aomFraction.*[p.litter_leaf_fractions(:).*(leafC-woodC);woodC]/1e4;
        aboveN=aomFraction.*[p.litter_leaf_fractions(:).*(leafN-woodN);woodN]/1e4;
        belowC=aomFraction.*[p.litter_root_fractions(:).*(rootC-mortalityRootC);mortalityRootC]/1e4;
        belowN=aomFraction.*[p.litter_root_fractions(:).*(rootN-mortalityRootN);mortalityRootN]/1e4;
        for j=1:4
            b=pool.components.(names{j});
            b.surface_C=b.surface_C+aboveC(j); b.surface_N=b.surface_N+aboveN(j);
            b.element_C=b.element_C+flip(belowC(j).*w./dzTop);
            b.element_N=b.element_N+flip(belowN(j).*w./dzTop);
            pool.components.(names{j})=b;
        end
        result.organic.aom.litter=soil.litter.componentTotals(pool,false);
    else
    result.organic.aom.litter.surface_C = ...
        result.organic.aom.litter.surface_C+aomFraction.*leafC/1e4;
    result.organic.aom.litter.surface_N = ...
        result.organic.aom.litter.surface_N+aomFraction.*leafN/1e4;
    result.organic.aom.litter.element_C = ...
        result.organic.aom.litter.element_C+flip(aomFraction.*cProfile);
    result.organic.aom.litter.element_N = ...
        result.organic.aom.litter.element_N+flip(aomFraction.*nProfile);
    end
    result.carbon_added_g_cm2=leafC/1e4+sum(cProfile.*dzTop);
    result.nitrogen_added_g_cm2=leafN/1e4+sum(nProfile.*dzTop);
    result.organic.cumulative.input_C_g_cm2 = ...
        result.organic.cumulative.input_C_g_cm2+result.carbon_added_g_cm2;
    result.organic.cumulative.input_N_g_cm2 = ...
        result.organic.cumulative.input_N_g_cm2+result.nitrogen_added_g_cm2;
    result.organic.cumulative=addCumulative(result.organic.cumulative, ...
        'input_litter_C_g_cm2',result.carbon_added_g_cm2);
    result.organic.cumulative=addCumulative(result.organic.cumulative, ...
        'input_litter_N_g_cm2',result.nitrogen_added_g_cm2);
    result.aboveground_carbon_g_m2=leafC;
    result.belowground_carbon_g_m2=rootC;
    result.aboveground_nitrogen_g_m2=leafN;
    result.belowground_nitrogen_g_m2=rootN;
    result.direct_doc_carbon_g_cm2=domFraction*result.carbon_added_g_cm2;
    result.direct_don_nitrogen_g_cm2=domFraction*result.nitrogen_added_g_cm2;
    result.mend_dom_input_fraction=domFraction;
    aomCarbon=aomFraction*result.carbon_added_g_cm2;
    aomNitrogen=aomFraction*result.nitrogen_added_g_cm2;
    assert(abs(aomCarbon+result.direct_doc_carbon_g_cm2 ...
        -result.carbon_added_g_cm2)<=1e-14*max(1,result.carbon_added_g_cm2) ...
        &&abs(aomNitrogen+result.direct_don_nitrogen_g_cm2 ...
        -result.nitrogen_added_g_cm2)<=1e-14*max(1,result.nitrogen_added_g_cm2), ...
        'stagePlantInput:mendInputClosure', ...
        'Source-MEND PlantCN input partition does not conserve C/N.');
    result.triggered=true;
end

function cumulative=addCumulative(cumulative,name,amount)
% Add a material-specific ledger field while accepting older checkpoints.
% Editor: Codex; Date: 2026-08-28.
    if ~isfield(cumulative,name), cumulative.(name)=0; end
    cumulative.(name)=cumulative.(name)+amount;
end

function [above,below]=splitWoodMortality(total,candidate,prefix)
% Gap mortality uses one common fractional rate for stem and coarse-root pools,
% so their post-mortality pool ratio exactly recovers the flux split.
% Editor: Jizixing Zhu; Date: 2026-08-24.
    abovePool=candidateLast(candidate,[prefix '_livestem']) ...
        +candidateLast(candidate,[prefix '_deadstem']);
    belowPool=candidateLast(candidate,[prefix '_livecroot']) ...
        +candidateLast(candidate,[prefix '_deadcroot']);
    if abovePool+belowPool<=0
        above=total; below=0; return;
    end
    below=total*belowPool/(abovePool+belowPool);
    above=total-below;
end

function value=candidateLast(candidate,name)
% Return the final candidate value for one PlantCN pool.
% Editor: Jizixing Zhu; Date: 2026-08-24.
    value=0;
    if isfield(candidate,name)&&isnumeric(candidate.(name)) ...
            &&~isempty(candidate.(name))&&isfinite(candidate.(name)(end))
        value=max(candidate.(name)(end),0);
    end
end

function value=increment(candidate,name,accepted)
% Compute one required accepted-to-candidate cumulative increment.
% Editor: Jizixing Zhu; Date: 2026-08-21.
    value=max(0,candidate.(name)-accepted);
end

function value=incrementOptional(candidate,name,acceptedAccum)
% Compute an optional accepted-to-candidate cumulative increment.
% Editor: Jizixing Zhu; Date: 2026-08-24.
    value=0;
    if isfield(candidate,name) && isfield(acceptedAccum,name)
        value=max(0,candidate.(name)-acceptedAccum.(name));
    end
end

function weights=rootWeights(dz,p)
% Build normalized root-distributed litter weights.
% Editor: Jizixing Zhu; Date: 2026-08-21.
    top=[0;cumsum(dz(1:end-1))]; bottom=cumsum(dz);
    weights=p.root_distribution_beta.^top-p.root_distribution_beta.^bottom;
    weights(bottom>p.root_max_depth_cm)=0;
    if sum(weights)<=0, weights(1)=1; end
    weights=weights/sum(weights);
end
