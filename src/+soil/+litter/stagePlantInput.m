function result=stagePlantInput(organicAccepted,microbeAccepted,domAccepted, ...
        plantAccepted,plantCandidate,built)
%STAGEPLANTINPUT Stage current-step PlantCN litter into the active soil owner.
%
%   Plant cumulative litter differences [g C/N m^-2] are converted to
%   [g cm^-2]. Leaf litter enters the top element; root litter follows the
%   configured beta depth distribution. COMPASS PlantCN supplies C/N amounts.
%   Leaf/fine-root returns enter metabolic/cellulose/lignin AOM; pruning and
%   wood/coarse-root mortality enter CWD. AOM stores are bottom-up g/cm3 bulk
%   soil with surface g/cm2. No direct MEND or DOM increment is made here.
%   AOM decomposition/leaching is advanced by the caller; no state is committed.
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
% Change: 2026-10-10 Wenjie Liu -- stage full CoLM plant C/N returns once
%   into MEND POM1/POM2/total DOM; preserve supplied total N with input C:N.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- route organ-specific PlantCN returns to
%   AOM litter components, retaining exact external C/N and cumulative ledgers.
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
    leafOnlyC=increment(plantCandidate,'cum_C_leaf_litter',a.cum_C_leaf_litter);
    leafOnlyN=increment(plantCandidate,'cum_N_leaf_litter',a.cum_N_leaf_litter);
    fineRootC=increment(plantCandidate,'cum_C_root_litter',a.cum_C_root_litter);
    fineRootN=increment(plantCandidate,'cum_N_root_litter',a.cum_N_root_litter);
    woodAboveC=incrementOptional(plantCandidate,'cum_C_prune_litter',a)+mortalityAboveC;
    woodAboveN=incrementOptional(plantCandidate,'cum_N_prune_litter',a)+mortalityAboveN;
    leafC=leafOnlyC+woodAboveC;rootC=fineRootC+mortalityRootC;
    leafN=leafOnlyN+woodAboveN;rootN=fineRootN+mortalityRootN;
    if leafC+rootC+leafN+rootN<=0, return; end
    dzTop=flip(built.domain.dz(:));
    w=rootWeights(dzTop,built.parameters.litter);
    dzBottom=built.domain.dz(:);
    topWeight=zeros(size(w));topWeight(1)=1;
    leafProfileC=flip((leafOnlyC/1e4).*topWeight./dzTop);
    leafProfileN=flip((leafOnlyN/1e4).*topWeight./dzTop);
    rootProfileC=flip((fineRootC/1e4).*w./dzTop);
    rootProfileN=flip((fineRootN/1e4).*w./dzTop);
    woodProfileC=flip((woodAboveC.*topWeight+mortalityRootC.*w)./1e4./dzTop);
    woodProfileN=flip((woodAboveN.*topWeight+mortalityRootN.*w)./1e4./dzTop);
    litter=result.organic.aom.litter;before=litter;
    if isfield(litter,'components')
        p=built.parameters.soil_organic.aom;
        leafFractions=validateFractions(p.litter_leaf_fractions);
        rootFractions=validateFractions(p.litter_root_fractions);
        names={'metabolic','cellulose','lignin'};
        for j=1:3
            k=names{j};
            litter.components.(k).element_C=litter.components.(k).element_C ...
                +leafFractions(j).*leafProfileC+rootFractions(j).*rootProfileC;
            litter.components.(k).element_N=litter.components.(k).element_N ...
                +leafFractions(j).*leafProfileN+rootFractions(j).*rootProfileN;
        end
        litter.components.cwd.element_C=litter.components.cwd.element_C+woodProfileC;
        litter.components.cwd.element_N=litter.components.cwd.element_N+woodProfileN;
        litter=soil.litter.componentTotals(litter,false);
    else
        litter.element_C=litter.element_C+leafProfileC+rootProfileC+woodProfileC;
        litter.element_N=litter.element_N+leafProfileN+rootProfileN+woodProfileN;
    end
    result.organic.aom.litter=litter;
    result.carbon_added_g_cm2=(leafC+rootC)/1e4;
    result.nitrogen_added_g_cm2=(leafN+rootN)/1e4;
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
    addedAomC=litter.surface_C-before.surface_C ...
        +sum((litter.element_C-before.element_C).*dzBottom);
    addedAomN=litter.surface_N-before.surface_N ...
        +sum((litter.element_N-before.element_N).*dzBottom);
    assert(abs(addedAomC ...
        -result.carbon_added_g_cm2)<=1e-14*max(1,result.carbon_added_g_cm2) ...
        &&abs(addedAomN ...
        -result.nitrogen_added_g_cm2)<=1e-14*max(1,result.nitrogen_added_g_cm2), ...
        'stagePlantInput:aomInputClosure', ...
        'PlantCN litter AOM increments do not conserve supplied C/N.');
    result.triggered=true;
end

function fractions=validateFractions(fractions)
% Validate organ-specific metabolic/cellulose/lignin input shares.
% Editor: Wenjie Liu; Date: 2026-10-17.
    assert(isnumeric(fractions)&&isreal(fractions)&&numel(fractions)==3 ...
        &&all(isfinite(fractions(:)))&&all(fractions(:)>=0) ...
        &&abs(sum(fractions)-1)<=1e-12,'stagePlantInput:aomFractions', ...
        'Litter chemical fractions must be nonnegative and sum to one.');
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
