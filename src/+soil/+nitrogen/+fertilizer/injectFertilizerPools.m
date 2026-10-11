function [SoilVariables, fert_diag] = injectFertilizerPools(SoilVariables, pools, DEPTHSTART, DEPTHEND, SUMDELTZ, NL, NN, DeltZ, NitrogenVariables)
%INJECTFERTILIZERPOOLS  Inject pre-resolved per-pool fertilizer N into soil concentration arrays.
%
%   [SoilVariables, fert_diag] = soil.nitrogen.fertilizer.injectFertilizerPools( ...
%       SoilVariables, pools, DEPTHSTART, DEPTHEND, SUMDELTZ, NL, NN, DeltZ, NitrogenVariables)
%   adds each pool's concentration over the [DEPTHSTART, DEPTHEND] injection band into
%   the matching state array. `pools` is the STRICT per-pool split from
%   soil.nitrogen.fertilizer.matchFertilizerEvent, so multiple coincident events with
%   different pools/ratios are allocated exactly (the legacy F_conc+fR aggregate could
%   not express that split). Solid NH4 (nh4_granule) enters the slow-release granule
%   pool; nh4_liquid enters dissolved NH4; nitrate and urea enter their fixed pools.
%
%   This is the new-chain counterpart of soil.nitrogen.applyFertilizer (which stays for
%   the legacy COMPASS.m path): same column-integral mass convention, but driven by an
%   explicit per-pool concentration struct instead of a controller + form ratios.
%
%   Inputs:
%     SoilVariables    [inout]  struct - C_NH4_new / C_NO3_new / C_Urea_new /
%                               C_NH4_granule [NN x1, ug/cm3] updated in place;
%                               .Theta_LL [NL x2] used for the mass integral
%     pools            [in]     struct(.nh4_granule,.nh4_liquid,.no3_liquid,.urea)
%                               [ug/cm3] concentration to inject per pool
%     DEPTHSTART/END   [in]     [cm]     injection band top/bottom (from surface)
%     SUMDELTZ         [in]     [cm, NL x1]  cumulative layer thickness from surface
%     NL, NN           [in]     [-]      layer / node counts
%     DeltZ            [in]     [cm, NL x1]  layer thickness
%     NitrogenVariables[in]     struct - .RHOKG [NL x1] bulk-density sorption factor
%
%   Returns:
%     SoilVariables    updated concentrations
%     fert_diag        struct(.applied,.nh4_added,.no3_added,.urea_added,.nh4_to_granule)
%                      [ug/cm2] column-integral masses added this step (same units as
%                      calculateStorage, so the N conservation ledger can treat them as
%                      external inputs; granule->dissolved release is internal).
%   Projection: integrate linear nodal basis functions over the exact fertilizer
%   band, then divide each allocated mass by its lumped storage capacity.
%   Urea/granules use bulk-volume capacity; dissolved pools use theta (+RHOKG).
%   Invalid bands or positive input into zero liquid capacity fail explicitly.
%
%   Reference: mass convention mirrors soil.nitrogen.applyFertilizer.m:51-114.
%   Editor: Jizixing Zhu; Date: 2026-07-02.
%   Change: 2026-09-07 Jizixing Zhu -- conservative band-to-node mass projection
%     and storage-weighted urea diagnostics on nonuniform grids.
%   Change: 2026-09-21 Wenjie Liu -- expose the exact nodal increment of
%     each management pool for passive source tracing; total pools are
%     unchanged and remain the only drivers of nitrogen dynamics.
    fert_diag = struct('applied', true, 'nh4_added', 0, 'no3_added', 0, ...
                       'urea_added', 0, 'nh4_to_granule', 0);

    dz = DeltZ(:);
    assert(NN==NL+1 && numel(dz)==NL && all(isfinite(dz)&dz>0), ...
        'fertilizer:badGrid','Invalid fertilizer grid.');
    z = [0;cumsum(flip(dz))];
    assert(numel(SUMDELTZ)==NL && max(abs(SUMDELTZ(:)-z(2:end)))<1e-8, ...
        'fertilizer:badGrid','Cumulative depths disagree with layer thicknesses.');
    assert(isfinite(DEPTHSTART) && isfinite(DEPTHEND) && ...
        DEPTHSTART>=0 && DEPTHEND>DEPTHSTART && DEPTHEND<=z(end), ...
        'fertilizer:badBand','Fertilizer band must lie within the soil column.');
    loadTop = zeros(NN,1);
    for j=1:NL
        a=max(z(j),DEPTHSTART); b=min(z(j+1),DEPTHEND);
        if b<=a, continue; end
        h=z(j+1)-z(j);
        deep=((b-z(j))^2-(a-z(j))^2)/(2*h);
        loadTop(j)=loadTop(j)+(b-a)-deep;
        loadTop(j+1)=loadTop(j+1)+deep;
    end
    loadNode=flip(loadTop); % [cm], sum = injection thickness
    bulkCapacity=[dz;0]/2+[0;dz]/2;

    Theta_LL = SoilVariables.Theta_LL;
    RHOKG = NitrogenVariables.RHOKG;
    [cv_nh4_pre, cv_no3_pre, cv_gran_pre] = columnStorage(SoilVariables, Theta_LL, RHOKG, NL, DeltZ);

    liquidCapacity=[dz.*Theta_LL(:,1);0]/2+[0;dz.*Theta_LL(:,2)]/2;
    nh4Capacity=liquidCapacity+[dz.*RHOKG(:);0]/2+[0;dz.*RHOKG(:)]/2;
    nh4Before=SoilVariables.C_NH4_new;
    no3Before=SoilVariables.C_NO3_new;
    ureaBefore=SoilVariables.C_Urea_new;
    granuleBefore=SoilVariables.C_NH4_granule;
    SoilVariables.C_NH4_granule = SoilVariables.C_NH4_granule + ...
        projectMass(pools.nh4_granule,loadNode,bulkCapacity);
    SoilVariables.C_NH4_new = SoilVariables.C_NH4_new + ...
        projectMass(pools.nh4_liquid,loadNode,nh4Capacity);
    SoilVariables.C_NO3_new = SoilVariables.C_NO3_new + ...
        projectMass(pools.no3_liquid,loadNode,liquidCapacity);
    SoilVariables.C_Urea_new = SoilVariables.C_Urea_new + ...
        projectMass(pools.urea,loadNode,bulkCapacity);

    [cv_nh4_post, cv_no3_post, cv_gran_post] = columnStorage(SoilVariables, Theta_LL, RHOKG, NL, DeltZ);
    % nh4_added = TOTAL NH4 injected (dissolved + granule); granule->dissolved is later
    % internal redistribution (Sdis) handled by nh4FertDissolution in the N solve.
    fert_diag.nh4_added      = (cv_nh4_post - cv_nh4_pre) + (cv_gran_post - cv_gran_pre);
    fert_diag.no3_added      = cv_no3_post - cv_no3_pre;
    fert_diag.nh4_to_granule = cv_gran_post - cv_gran_pre;
    fert_diag.urea_added = bulkCapacity'*(SoilVariables.C_Urea_new-ureaBefore);
    fert_diag.tracer_increment=struct( ...
        'nh4_liquid',SoilVariables.C_NH4_new-nh4Before, ...
        'no3_liquid',SoilVariables.C_NO3_new-no3Before, ...
        'urea',SoilVariables.C_Urea_new-ureaBefore, ...
        'nh4_granule',SoilVariables.C_NH4_granule-granuleBefore, ...
        'units','ug N cm-3 corresponding pool', ...
        'grid_order','bottom_up_nodes');
end

function delta = projectMass(concentration,loadNode,capacity)
% Convert allocated nodal mass [ug/cm2] to state concentration [ug/cm3].
% Editor: Jizixing Zhu; Date: 2026-09-07.
    assert(isscalar(concentration)&&isfinite(concentration)&&concentration>=0, ...
        'fertilizer:badInput','Fertilizer concentration must be finite and nonnegative.');
    delta=zeros(size(loadNode));
    selected=loadNode>0 & concentration>0;
    assert(all(isfinite(capacity(selected)) & capacity(selected)>0), ...
        'fertilizer:zeroCapacity','Cannot inject dissolved fertilizer into zero capacity.');
    delta(selected)=concentration.*loadNode(selected)./capacity(selected);
end

% ---------------------------------------------------------------------------
function [cv_nh4, cv_no3, cv_gran] = columnStorage(sv, Theta_LL, RHOKG, NL, DeltZ)
% Column-integral N storage [ug/cm2], matching calculateStorage / applyFertilizer:
% dissolved NH4 weighted by (Theta_LL + RHOKG) (sorption), NO3 by Theta_LL, granule by
% thickness only (surface-pool concentration, no water scaling).
% Editor: Jizixing Zhu; Date: 2026-07-02.
    cv_nh4 = 0; cv_no3 = 0; cv_gran = 0;
    for ML = 1:NL
        cv_nh4 = cv_nh4 + DeltZ(ML) * ( ...
            (Theta_LL(ML,1) + RHOKG(ML)) * sv.C_NH4_new(ML) + ...
            (Theta_LL(ML,2) + RHOKG(ML)) * sv.C_NH4_new(ML+1)) / 2;
        cv_no3 = cv_no3 + DeltZ(ML) * ( ...
            Theta_LL(ML,1) * sv.C_NO3_new(ML) + Theta_LL(ML,2) * sv.C_NO3_new(ML+1)) / 2;
        cv_gran = cv_gran + DeltZ(ML) * ( ...
            sv.C_NH4_granule(ML) + sv.C_NH4_granule(ML+1)) / 2;
    end
end
