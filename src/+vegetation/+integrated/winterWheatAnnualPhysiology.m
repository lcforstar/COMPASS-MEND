function [C, N, memory, diag] = winterWheatAnnualPhysiology( ...
    C, N, acceptedC, acceptedN, phen, env, physiology, memory, growthRespFraction, grainCN, dt)
%WINTERWHEATANNUALPHYSIOLOGY Conservative standing senescence and stem WSC.
% Pools [g C or N m^-2] are caller-owned state; phen/env are forcing;
% physiology, growthRespFraction and grainCN are compiled parameters.
% memory is accepted canopy-area event memory, returned as pending state.
% acceptedC/N are start-of-step donors after any explicit grain-fill N
% remobilization. New tissue is ineligible for same-round senescence.
% diag contains applied internal transfers and explicit growth CO2 [g C m^-2].
% Reference mechanisms and adaptation limits:
% docs/plant/winter_wheat_annual_physiology.md (APSIM; Takahashi et al. 2001).
% Editor: Wenjie Liu; Date: 2026-09-30.
    validateattributes(dt, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(growthRespFraction, {'numeric'}, ...
        {'scalar','real','finite','>=',0,'<',1});
    validateattributes(grainCN, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(physiology.senescent_leaf_cn, {'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(physiology.stem_wsc_retrans_daily_fraction, {'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',1});
    fieldsC = {'C_leaf','C_leaf_senescent','C_stem_wsc','C_fruit'};
    fieldsN = {'N_leaf','N_leaf_senescent','N_labile','N_fruit'};
    for i=1:numel(fieldsC)
        validateattributes(C.(fieldsC{i}), {'numeric'}, ...
            {'scalar','real','finite','nonnegative'});
    end
    for i=1:numel(fieldsN)
        validateattributes(N.(fieldsN{i}), {'numeric'}, ...
            {'scalar','real','finite','nonnegative'});
    end
    validateattributes(acceptedC.C_leaf, {'numeric'}, ...
        {'scalar','real','finite','nonnegative'});
    validateattributes(acceptedC.C_stem_wsc, {'numeric'}, ...
        {'scalar','real','finite','nonnegative'});
    validateattributes(acceptedN.N_leaf, {'numeric'}, ...
        {'scalar','real','finite','nonnegative'});
    validateattributes(memory.crop_green_lai_previous, {'numeric'}, ...
        {'scalar','real','finite','nonnegative'});
    assert(isscalar(memory.crop_green_lai_initialized) ...
        && ismember(memory.crop_green_lai_initialized,[false true]), ...
        'winterWheatAnnualPhysiology:invalidMemory', ...
        'Accepted crop area memory must have a logical initialization flag.');
    diag = struct('C_leaf_senesced',0,'N_leaf_senesced',0, ...
        'N_senescence_to_labile',0,'C_wsc_to_grain',0, ...
        'N_wsc_to_grain',0,'C_wsc_growth_resp',0);
    lai = env.LAI;
    validateattributes(lai, {'numeric'}, {'scalar','real','finite','nonnegative'});
    harvest = isfield(phen,'crop_harvest_event') && phen.crop_harvest_event;
    planting = isfield(phen,'crop_planting_event') && phen.crop_planting_event;
    if planting
        memory.crop_green_lai_initialized = false;
        memory.crop_green_lai_previous = 0;
    end
    if ~harvest && phen.crop_live && memory.crop_green_lai_initialized ...
            && memory.crop_green_lai_previous > 0
        fraction = min(1,max(0,1-lai/memory.crop_green_lai_previous));
        % Only start-of-step leaf C is eligible: newly synthesized tissue is
        % not created and senesced within the same allocation round.
        dC = min(C.C_leaf, fraction*acceptedC.C_leaf);
        dN = min(N.N_leaf, fraction*acceptedN.N_leaf);
        retainedN = min(dN,dC/physiology.senescent_leaf_cn);
        C.C_leaf = C.C_leaf-dC;
        C.C_leaf_senescent = C.C_leaf_senescent+dC;
        N.N_leaf = N.N_leaf-dN;
        N.N_leaf_senescent = N.N_leaf_senescent+retainedN;
        N.N_labile = N.N_labile+dN-retainedN;
        diag.C_leaf_senesced = dC;
        diag.N_leaf_senesced = dN;
        diag.N_senescence_to_labile = dN-retainedN;
    end
    mature = isfield(phen,'crop_mature_flag') && phen.crop_mature_flag;
    if ~harvest && phen.crop_live && phen.crop_phase == 3 && ~mature
        fraction = -expm1(log1p(-physiology.stem_wsc_retrans_daily_fraction)*dt/86400);
        supply = min(C.C_stem_wsc, acceptedC.C_stem_wsc)*fraction;
        yield = 1-growthRespFraction;
        existingCapacity = max(0, grainCN*N.N_fruit-C.C_fruit);
        sink = existingCapacity + grainCN*N.N_labile;
        debit = min(supply,sink/yield);
        grainC = debit*yield;
        grainN = min(N.N_labile,max(0,grainC-existingCapacity)/grainCN);
        respiration = debit-grainC;
        C.C_stem_wsc = C.C_stem_wsc-debit;
        C.C_fruit = C.C_fruit+grainC;
        N.N_labile = N.N_labile-grainN;
        N.N_fruit = N.N_fruit+grainN;
        diag.C_wsc_to_grain = grainC;
        diag.N_wsc_to_grain = grainN;
        diag.C_wsc_growth_resp = respiration;
    end
    % A zero before first positive canopy must leave the seed leaf intact.
    % Harvest resets the memory, preventing across-season area comparisons.
    if harvest
        memory.crop_green_lai_initialized = false;
        memory.crop_green_lai_previous = 0;
    elseif phen.crop_live && lai > 0
        memory.crop_green_lai_initialized = true;
        memory.crop_green_lai_previous = lai;
    elseif memory.crop_green_lai_initialized
        memory.crop_green_lai_previous = lai;
    end
end
