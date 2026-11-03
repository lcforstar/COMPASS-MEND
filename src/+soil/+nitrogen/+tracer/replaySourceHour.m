function [tagAfter,exchange,routes]=replaySourceHour( ...
        tagBefore,totalBefore,sourceLedger,rootDin,totalAfter,bulkDensityBottom)
%REPLAYSOURCEHOUR Replay accepted source-MEND gross N routes for one label.
%   [TAGAFTER,EXCHANGE,ROUTES] = REPLAYSOURCEHOUR(TAGBEFORE,TOTALBEFORE,
%   SOURCELEDGER,ROOTDIN,TOTALAFTER,BULKDENSITYBOTTOM) advances a passive
%   fertilizer-N label through the already accepted production
%   source_hourly MEND reaction. It never evaluates a rate and cannot feed
%   back to total nitrogen dynamics.
%
%   TAGBEFORE, TOTALBEFORE and TOTALAFTER contain:
%     .organic_n             COMPASS organic-N owner, bottom-up [g N cm-3]
%     .don_total_g_cm3       sole total DON owner, bottom-up
%     .don_sorbed_g_cm3      dynamic sorbed subpool, when mend_dynamic
%   The organic owner must contain PON1/PON2/MON, MB/BD, ENZP1/ENZP2/EM,
%   ENZN, and the mobile/attached/retained phases of MB, BD and ENZN.
%   QON, when present, must be identically zero because equilibrium sorbed
%   DON is algebraic rather than a second prognostic owner.
%
%   SOURCELEDGER is the accepted top-down mend.source_hour_ledger. Nitrogen
%   routes use mg N/g dry soil; enzyme_production_C/enzyme_turnover_C and
%   their enzyme_production_by_phase_C/enzyme_turnover_by_phase_C arrays
%   use mg C/g dry soil. The latter are required NL-by-3 applied gross
%   amounts ordered [mobile, attached, retained], converted using enzyme_cn
%   in the accepted ledger (missing in archived ledgers defaults to source 3).
%   Endpoint differences are never a fallback for missing phase ledgers.
%   ROOTDIN.total_mg_g and ROOTDIN.tag_mg_g are
%   NL-by-2 top-down post-root inventories ordered [NH4, NO3].
%   BULKDENSITYBOTTOM is NL-by-1 [g dry soil cm-3 bulk soil]. TAGAFTER has
%   the same owner units/order as TAGBEFORE. EXCHANGE and ROUTES use
%   top-down [mg N g-1 dry soil].
%
%   Labeled transfers use accepted gross amounts and donor label fractions.
%   Dynamic sorption is the explicit exception to start-donor attribution:
%   it uses bounded implicit same-hour mixing after non-sorption DON routes.
%   First DIN
%   uptake precedes MBA/MBD stoichiometric release; second uptake uses the
%   resulting species-specific DIN fractions. Intracellular-enzyme
%   production/turnover is replayed by cell phase. Net pool changes and
%   final-state label fractions are never used to infer route attribution.
%   Inputs outside the source_hourly + supported DOM sorption contract,
%   nonfinite values, overdraws, and label bounds fail without clipping.
%
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-09-21 Wenjie Liu -- preserve accepted MBA release when
% oxygen scaling produces a ULP-level negative reconstructed surplus.
% Change: 2026-09-22 Wenjie Liu -- bounded complementary attribution of
% aggregate reaction labels to immutable resized microbial phase totals.
% Change: 2026-10-12 Wenjie Liu -- consume the producer's applied gross
% enzyme phase ledger; retain endpoint/aggregate audits without subtraction
% of almost equal enzyme inventories to infer a small positive synthesis.
% Change: 2026-10-28 Wenjie Liu -- replay dynamic free/sorbed DON gross
% exchanges using bounded simultaneous mixing; sorbed N remains within DON.

    [tag0,total0,total1,bdTop,n]=normalizeOwners( ...
        tagBefore,totalBefore,totalAfter,bulkDensityBottom);
    [dinTotal,dinTag]=validateDin(rootDin,n);
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- replay configured enzyme C:N metadata.
    ledger=validateLedger(sourceLedger,n);
    sourceLedger=ledger;
    nf=ledger.nitrogen_flux;

    dynamic=strcmp(ownerMode(totalBefore),'mend_dynamic');
    if ~dynamic
        assertZeroFlux(nf,'DOM2QOM');
        assertZeroFlux(nf,'QOM2DOM');
        assertZeroFlux(nf,'DOM2QOMnet');
    end

    % Resolve the source's signed MBA compensation into the same two gross
    % phase directions used by carryReactionPhases. The reverse transfer is
    % evaluated after active-to-dormant transfer, matching production order.
    mba2mbd=getAmount(nf,'MBA2MBD',n,true);
    mbd2mba=getAmount(nf,'MBD2MBA',n,false);
    activeToDormant=max(mba2mbd,0);
    dormantToActive=mbd2mba+max(-mba2mbd,0);
    [phaseWork,phaseRoute]=switchPhases(total0.phase,tag0.phase, ...
        activeToDormant,dormantToActive);

    % Non-sorption Euler routes use reaction-start donor label fractions.
    % Dynamic sorption instead uses accepted same-hour implicit mixing below.
    % Separate destination updates retain simultaneous opposing transfers.
    tf=struct();
    tf.POMdec1=tagRoute(getAmount(nf,'POMdec1',n,false), ...
        tag0.POM1,total0.POM1,'POM1 decomposition');
    tf.POMdec2=tagRoute(getAmount(nf,'POMdec2',n,false), ...
        tag0.POM2,total0.POM2,'POM2 decomposition');
    tf.POMdec2DOM1=tagRoute(getAmount(nf,'POMdec2DOM1',n,false), ...
        tag0.POM1,total0.POM1,'POM1 to DON');
    tf.POMdec2DOM2=tagRoute(getAmount(nf,'POMdec2DOM2',n,false), ...
        tag0.POM2,total0.POM2,'POM2 to DON');
    tf.POMdec2MOM1=tagRoute(getAmount(nf,'POMdec2MOM1',n,false), ...
        tag0.POM1,total0.POM1,'POM1 to MOM');
    tf.POMdec2MOM2=tagRoute(getAmount(nf,'POMdec2MOM2',n,false), ...
        tag0.POM2,total0.POM2,'POM2 to MOM');
    tf.MOM2DOM=tagRoute(getAmount(nf,'MOM2DOM',n,false), ...
        tag0.MOM,total0.MOM,'MOM to DON');
    tf.DOM2MBA=tagRoute(getAmount(nf,'DOM2MBA',n,false), ...
        tag0.DOM,total0.DOM,'DON to MBA');
    tf.DOM2QOM=zeros(n,1);tf.QOM2DOM=zeros(n,1);
    tf.DOM2QOMnet=tf.DOM2QOM-tf.QOM2DOM;
    tf.MBA2DOM=tagRoute(getAmount(nf,'MBA2DOM',n,false), ...
        tag0.MBA,total0.MBA,'MBA to DON');
    tf.MBA2POM1=tagRoute(getAmount(nf,'MBA2POM1',n,false), ...
        tag0.MBA,total0.MBA,'MBA to POM1');
    tf.MBA2POM2=tagRoute(getAmount(nf,'MBA2POM2',n,false), ...
        tag0.MBA,total0.MBA,'MBA to POM2');
    tf.MBA2EP1=tagRoute(getAmount(nf,'MBA2EP1',n,false), ...
        tag0.MBA,total0.MBA,'MBA to ENZP1');
    tf.MBA2EP2=tagRoute(getAmount(nf,'MBA2EP2',n,false), ...
        tag0.MBA,total0.MBA,'MBA to ENZP2');
    tf.MBA2EM=tagRoute(getAmount(nf,'MBA2EM',n,false), ...
        tag0.MBA,total0.MBA,'MBA to ENZM');
    tf.EP2DOM1=tagRoute(getAmount(nf,'EP2DOM1',n,false), ...
        tag0.ENZP1,total0.ENZP1,'ENZP1 to DON');
    tf.EP2DOM2=tagRoute(getAmount(nf,'EP2DOM2',n,false), ...
        tag0.ENZP2,total0.ENZP2,'ENZP2 to DON');
    tf.EM2DOM=tagRoute(getAmount(nf,'EM2DOM',n,false), ...
        tag0.ENZM,total0.ENZM,'ENZM to DON');
    tf.MBA_to_MBD=phaseRoute.tag_active_to_dormant;
    tf.MBD_to_MBA=phaseRoute.tag_dormant_to_active;

    validateSourceIdentities(nf,n);

    % Kinetic first uptake sees the post-root inventory at reaction start.
    firstTotal=[getAmount(sourceLedger,'first_NH4',n,false), ...
        getAmount(sourceLedger,'first_NO3',n,false)];
    firstTag=tagDinWithdrawal(firstTotal,dinTotal,dinTag,'first DIN uptake');

    phi=getFraction(sourceLedger,'phi_n',n);
    domToMbaTotal=getAmount(nf,'DOM2MBA',n,false);
    directReleaseEstimate=(1-phi).*domToMbaTotal;
    mbaReleaseTotal=getAmount(sourceLedger,'Nmn_MBA',n,false);
    [directReleaseTotal,surplusReleaseTotal]=splitAcceptedMbaRelease( ...
        mbaReleaseTotal,directReleaseEstimate);
    directReleaseTag=tagRoute(directReleaseTotal,tf.DOM2MBA, ...
        domToMbaTotal,'MBA direct release');

    % Replay pool changes before the stoichiometric overflow/release stage.
    t=struct(); q=struct();
    [t.POM1,q.POM1]=pairedUpdate(total0.POM1,tag0.POM1, ...
        -getAmount(nf,'POMdec1',n,false)+getAmount(nf,'MBA2POM1',n,false), ...
        -tf.POMdec1+tf.MBA2POM1,'POM1');
    [t.POM2,q.POM2]=pairedUpdate(total0.POM2,tag0.POM2, ...
        -getAmount(nf,'POMdec2',n,false)+getAmount(nf,'MBA2POM2',n,false), ...
        -tf.POMdec2+tf.MBA2POM2,'POM2');
    [t.MOM,q.MOM]=pairedUpdate(total0.MOM,tag0.MOM, ...
        getAmount(nf,'POMdec2MOM1',n,false)+getAmount(nf,'POMdec2MOM2',n,false) ...
        -getAmount(nf,'MOMdec',n,false), ...
        tf.POMdec2MOM1+tf.POMdec2MOM2-tf.MOM2DOM,'MOM');
    [t.DOM,q.DOM]=pairedUpdate(total0.DOM,tag0.DOM, ...
        getAmount(nf,'POMdec2DOM1',n,false)+getAmount(nf,'POMdec2DOM2',n,false) ...
        +getAmount(nf,'MOM2DOM',n,false)+getAmount(nf,'MBA2DOM',n,false) ...
        +getAmount(nf,'EP2DOM1',n,false)+getAmount(nf,'EP2DOM2',n,false) ...
        +getAmount(nf,'EM2DOM',n,false)-getAmount(nf,'DOM2MBA',n,false), ...
        tf.POMdec2DOM1+tf.POMdec2DOM2+tf.MOM2DOM+tf.MBA2DOM ...
        +tf.EP2DOM1+tf.EP2DOM2+tf.EM2DOM-tf.DOM2MBA,'DON');
    t.QOM=total0.QOM;q.QOM=tag0.QOM;sorptionRoundoff=zeros(n,1);
    if dynamic
        % The accepted source map permits sorption borrowing from all DON
        % inputs within this hour, and checks only its final phase stocks.
        % Include those gross inputs in the implicit label mixing inventory;
        % rejecting the pre-input intermediate would reject valid totals.
        [tf.DOM2QOM,tf.QOM2DOM,freeLabel,sorbedLabel,sorptionRoundoff]=sorptionLabels( ...
            t.DOM,q.DOM,total0.QOM,tag0.QOM, ...
            getAmount(nf,'DOM2QOM',n,true),getAmount(nf,'QOM2DOM',n,true), ...
            total1.DOM,total1.QOM);
        tf.DOM2QOMnet=tf.DOM2QOM-tf.QOM2DOM;
        assertClose(getAmount(nf,'DOM2QOMnet',n,true), ...
            getAmount(nf,'DOM2QOM',n,true)-getAmount(nf,'QOM2DOM',n,true), ...
            'dynamic sorption gross net');
        t.DOM=total1.DOM;t.QOM=total1.QOM;
        q.DOM=freeLabel;q.QOM=sorbedLabel;
    end
    [t.ENZP1,q.ENZP1]=pairedUpdate(total0.ENZP1,tag0.ENZP1, ...
        getAmount(nf,'MBA2EP1',n,false)-getAmount(nf,'EP2DOM1',n,false), ...
        tf.MBA2EP1-tf.EP2DOM1,'ENZP1');
    [t.ENZP2,q.ENZP2]=pairedUpdate(total0.ENZP2,tag0.ENZP2, ...
        getAmount(nf,'MBA2EP2',n,false)-getAmount(nf,'EP2DOM2',n,false), ...
        tf.MBA2EP2-tf.EP2DOM2,'ENZP2');
    [t.ENZM,q.ENZM]=pairedUpdate(total0.ENZM,tag0.ENZM, ...
        getAmount(nf,'MBA2EM',n,false)-getAmount(nf,'EM2DOM',n,false), ...
        tf.MBA2EM-tf.EM2DOM,'ENZM');

    mortalityTotal=getAmount(nf,'MBA2DOM',n,false) ...
        +getAmount(nf,'MBA2POM1',n,false)+getAmount(nf,'MBA2POM2',n,false);
    mortalityTag=tf.MBA2DOM+tf.MBA2POM1+tf.MBA2POM2;
    mbOtherLossTotal=mortalityTotal+getAmount(nf,'MBA2EP1',n,false) ...
        +getAmount(nf,'MBA2EP2',n,false)+getAmount(nf,'MBA2EM',n,false);
    mbOtherLossTag=mortalityTag+tf.MBA2EP1+tf.MBA2EP2+tf.MBA2EM;

    mbBeforeSurplusTotal=total0.MBA+getAmount(nf,'DOM2MBA',n,false) ...
        -mbOtherLossTotal-activeToDormant+dormantToActive ...
        +sum(firstTotal,2)-directReleaseTotal;
    mbBeforeSurplusTag=tag0.MBA+tf.DOM2MBA-mbOtherLossTag ...
        -tf.MBA_to_MBD+tf.MBD_to_MBA+sum(firstTag,2)-directReleaseTag;
    validatePair(mbBeforeSurplusTotal,mbBeforeSurplusTag,'MBA before surplus release');

    mbdBeforeReleaseTotal=total0.MBD+activeToDormant-dormantToActive;
    mbdBeforeReleaseTag=tag0.MBD+tf.MBA_to_MBD-tf.MBD_to_MBA;
    validatePair(mbdBeforeReleaseTotal,mbdBeforeReleaseTag,'MBD before release');

    mbdReleaseTotal=getAmount(sourceLedger,'Nmn_MBD',n,false);
    mbdReleaseTag=tagRoute(mbdReleaseTotal,mbdBeforeReleaseTag, ...
        mbdBeforeReleaseTotal,'MBD stoichiometric release');
    t.MBD=mbdBeforeReleaseTotal-mbdReleaseTotal;
    q.MBD=mbdBeforeReleaseTag-mbdReleaseTag;
    validatePair(t.MBD,q.MBD,'MBD after release');

    surplusReleaseTag=tagRoute(surplusReleaseTotal,mbBeforeSurplusTag, ...
        mbBeforeSurplusTotal,'MBA surplus release');
    mbAfterReleaseTotal=mbBeforeSurplusTotal-surplusReleaseTotal;
    mbAfterReleaseTag=mbBeforeSurplusTag-surplusReleaseTag;
    validatePair(mbAfterReleaseTotal,mbAfterReleaseTag,'MBA after release');

    % Release occurs before source stoichiometric second uptake. Thus this
    % stage deliberately uses the updated DIN label fraction, not the start
    % or final fraction.
    releaseTotal=directReleaseTotal+surplusReleaseTotal+mbdReleaseTotal;
    releaseTag=directReleaseTag+surplusReleaseTag+mbdReleaseTag;
    dinAvailableTotal=[dinTotal(:,1)-firstTotal(:,1)+releaseTotal, ...
        dinTotal(:,2)-firstTotal(:,2)];
    dinAvailableTag=[dinTag(:,1)-firstTag(:,1)+releaseTag, ...
        dinTag(:,2)-firstTag(:,2)];
    validatePair(dinAvailableTotal,dinAvailableTag,'DIN before second uptake');
    secondTotal=[getAmount(sourceLedger,'second_NH4',n,false), ...
        getAmount(sourceLedger,'second_NO3',n,false)];
    secondTag=tagDinWithdrawal(secondTotal,dinAvailableTotal,dinAvailableTag, ...
        'second DIN uptake');
    dinEndTotal=dinAvailableTotal-secondTotal;
    dinEndTag=dinAvailableTag-secondTag;
    validatePair(dinEndTotal,dinEndTag,'DIN after source reaction');

    mbBeforeEnzymeTotal=mbAfterReleaseTotal+sum(secondTotal,2);
    mbBeforeEnzymeTag=mbAfterReleaseTag+sum(secondTag,2);
    validatePair(mbBeforeEnzymeTotal,mbBeforeEnzymeTag,'MBA before ENZN exchange');

    % Mirror carryReactionPhases for N: explicit switch identity first, then
    % total-group resize to the gross-derived aggregate. Attribute the new
    % aggregate label with bounded complementary phase capacities; independent
    % rescaling of a heterogeneous label can exceed a phase's total inventory.
    [mbPhaseTotal,mbPhaseTag,mbAllocation]=resizePhaseGroup( ...
        phaseWork.total.MBA,phaseWork.tag.MBA, ...
        mbBeforeEnzymeTotal,mbBeforeEnzymeTag,true,'MBA');
    [mbdPhaseTotal,mbdPhaseTag,mbdAllocation]=resizePhaseGroup( ...
        phaseWork.total.MBD,phaseWork.tag.MBD,t.MBD,q.MBD,false,'MBD');

    productionN=getAmount(sourceLedger,'enzyme_production_C',n,false)./ledger.enzyme_cn;
    turnoverN=getAmount(sourceLedger,'enzyme_turnover_C',n,false)./ledger.enzyme_cn;
    [productionPhase,turnoverPhase]=enzymePhaseAmounts( ...
        sourceLedger,n,total0.phase.ENZN,total1.phase.ENZN,productionN,turnoverN);
    expectedMbBeforePhase=total1.phase.MBA+productionPhase-turnoverPhase;
    assertClose(mbPhaseTotal,expectedMbBeforePhase,'MBA phase carry');
    assertClose(mbdPhaseTotal,total1.phase.MBD,'MBD phase carry');

    productionTagPhase=tagRoute(productionPhase,mbPhaseTag,mbPhaseTotal, ...
        'MBA to intracellular enzyme by phase');
    turnoverTagPhase=tagRoute(turnoverPhase,tag0.phase.ENZN, ...
        total0.phase.ENZN,'intracellular enzyme to MBA by phase');
    mbPhaseTag=mbPhaseTag-productionTagPhase+turnoverTagPhase;
    enznPhaseTag=tag0.phase.ENZN+productionTagPhase-turnoverTagPhase;
    validatePair(total1.phase.MBA,mbPhaseTag,'MBA final phases');
    validatePair(total1.phase.MBD,mbdPhaseTag,'MBD final phases');
    validatePair(total1.phase.ENZN,enznPhaseTag,'ENZN final phases');

    t.MBA=sum(total1.phase.MBA,2);
    q.MBA=sum(mbPhaseTag,2);
    t.MBD=sum(total1.phase.MBD,2);
    q.MBD=sum(mbdPhaseTag,2);
    t.ENZN=sum(total1.phase.ENZN,2);
    q.ENZN=sum(enznPhaseTag,2);
    assertClose(t.MBA,total1.MBA,'MBA final aggregate');
    assertClose(t.MBD,total1.MBD,'MBD final aggregate');
    assertClose(t.ENZN,total1.ENZN,'ENZN final aggregate');
    validatePair(t.MBA,q.MBA,'MBA final aggregate tag');
    validatePair(t.MBD,q.MBD,'MBD final aggregate tag');
    validatePair(t.ENZN,q.ENZN,'ENZN final aggregate tag');

    % Every other total endpoint is independently reconstructed from gross
    % accepted routes. This is a contract check, not a source of tag values.
    for name={'POM1','POM2','MOM','DOM','QOM','ENZP1','ENZP2','ENZM'}
        key=name{1};assertClose(t.(key),total1.(key),[key ' gross replay']);
        validatePair(total1.(key),q.(key),[key ' final tag']);
    end

    tagAfter=writeOwners(tagBefore,q,mbPhaseTag,mbdPhaseTag,enznPhaseTag,bdTop);

    tf.ENZN_production=sum(productionTagPhase,2);
    tf.ENZN_turnover=sum(turnoverTagPhase,2);
    tf.ENZN_production_by_phase=productionTagPhase;
    tf.ENZN_turnover_by_phase=turnoverTagPhase;
    tf.MBA_direct_release=directReleaseTag;
    tf.MBA_surplus_release=surplusReleaseTag;
    tf.MBD_release=mbdReleaseTag;
    tf.first_NH4=firstTag(:,1);tf.first_NO3=firstTag(:,2);
    tf.second_NH4=secondTag(:,1);tf.second_NO3=secondTag(:,2);

    beforeInventory=poolInventory(tag0)+sum(dinTag,2);
    afterInventory=poolInventory(q)+sum(dinEndTag,2);
    residual=afterInventory-beforeInventory;
    assert(all(abs(residual)<=5e-11.*max(1,abs(beforeInventory))), ...
        'soil:nitrogen:tracer:sourceClosure', ...
        'Fertilizer-N source-hour gross replay does not close by layer.');

    exchange=struct( ...
        'first_uptake_mg_g',firstTag, ...
        'second_uptake_mg_g',secondTag, ...
        'total_uptake_mg_g',firstTag+secondTag, ...
        'mba_direct_release_mg_g',directReleaseTag, ...
        'mba_surplus_release_mg_g',surplusReleaseTag, ...
        'mbd_release_mg_g',mbdReleaseTag, ...
        'total_nh4_release_mg_g',releaseTag, ...
        'net_nh4_sink_mg_g',firstTag(:,1)+secondTag(:,1)-releaseTag, ...
        'net_no3_sink_mg_g',firstTag(:,2)+secondTag(:,2), ...
        'din_after_tag_mg_g',dinEndTag, ...
        'din_after_total_mg_g',dinEndTotal, ...
        'units','mg_N_g_dry_soil','grid_order','top_down');
    routes=struct('nitrogen_flux',tf, ...
        'sorption_attribution_roundoff_mg_g',sorptionRoundoff, ...
        'phase_switch',phaseRoute, ...
        'phase_attribution',struct('MBA',mbAllocation,'MBD',mbdAllocation), ...
        'closure_residual_mg_g',residual, ...
        'units','mg_N_g_dry_soil','grid_order','top_down', ...
        'reaction_method','source_hourly', ...
        'sorption_model',ownerMode(totalBefore), ...
        'attribution','accepted_gross_flux_times_stage_donor_tag_fraction');
end

function [tag0,total0,total1,bdTop,n]=normalizeOwners(tag,before,after,bdBottom)
% Convert sole COMPASS owners into source top-down mg/g working pools.
    assert(isstruct(tag)&&isfield(tag,'organic_n')&&isfield(tag,'don_total_g_cm3') ...
        &&isstruct(before)&&isfield(before,'organic_n')&&isfield(before,'don_total_g_cm3') ...
        &&isstruct(after)&&isfield(after,'organic_n')&&isfield(after,'don_total_g_cm3'), ...
        'soil:nitrogen:tracer:sourceOwner', ...
        'Tracer and total source states require organic_n and don_total_g_cm3 owners.');
    bd=bdBottom(:);n=numel(bd);
    assert(n>0&&all(isfinite(bd))&&all(bd>0), ...
        'soil:nitrogen:tracer:bulkDensity','Bulk density must be finite and positive.');
    bdTop=flip(bd);
    tag0=readOwner(tag,bdTop,n,'tag');
    total0=readOwner(before,bdTop,n,'total before');
    total1=readOwner(after,bdTop,n,'total after');
    validatePoolBounds(total0,tag0,'accepted fertilizer tag');
    mode=ownerMode(before);
    assert(ismember(mode,{'linear_equilibrium','mend_dynamic'}) ...
        &&strcmp(mode,ownerMode(tag))&&strcmp(mode,ownerMode(after)), ...
        'soil:nitrogen:tracer:sorptionModel', ...
        'Source-hour total and label sorption modes must agree.');
end

function out=readOwner(bundle,bdTop,n,label)
    o=bundle.organic_n;
    map={ ...
        'POM1','PON1';'POM2','PON2';'MOM','MON'; ...
        'MBA','MB_N';'MBD','BD_N';'ENZP1','ENZP1_N'; ...
        'ENZP2','ENZP2_N';'ENZM','EM_N';'ENZN','ENZN_N'};
    out=struct();
    for j=1:size(map,1)
        field=map{j,2};
        out.(map{j,1})=toMgG(requiredColumn(o,field,n,label),bdTop);
    end
    out.DOM=toMgG(requiredColumn(bundle,'don_total_g_cm3',n,label),bdTop);
    out.QOM=zeros(n,1);
    if strcmp(ownerMode(bundle),'mend_dynamic')
        out.QOM=toMgG(requiredColumn(bundle,'don_sorbed_g_cm3',n,label),bdTop);
        assert(all(out.QOM<=out.DOM),'soil:nitrogen:tracer:dynamicPhase', ...
            'Sorbed DON must not exceed total DON.');
        out.DOM=out.DOM-out.QOM;
    end
    if isfield(o,'QON')
        qon=requiredColumn(o,'QON',n,label);
        assert(all(qon==0),'soil:nitrogen:tracer:qonOwner', ...
            'QON must remain zero in linear-equilibrium source replay.');
    end
    names={'mobile','attached','retained'};
    groups={'MB','BD','ENZN'};
    for g=1:numel(groups)
        matrix=zeros(n,3);
        for p=1:3
            field=[groups{g} '_' names{p} '_N'];
            matrix(:,p)=toMgG(requiredColumn(o,field,n,label),bdTop);
        end
        out.phase.(sourceGroup(groups{g}))=matrix;
        assertClose(sum(matrix,2),out.(sourceGroup(groups{g})), ...
            [label ' ' groups{g} ' phase sum']);
    end
    values=poolInventory(out);
    assert(all(isfinite(values))&&all(values>=0), ...
        'soil:nitrogen:tracer:ownerValue','%s owner contains invalid N.',label);
end

function name=sourceGroup(name)
    if strcmp(name,'MB'),name='MBA';elseif strcmp(name,'BD'),name='MBD';end
end

function v=requiredColumn(s,name,n,label)
    assert(isstruct(s)&&isfield(s,name),'soil:nitrogen:tracer:ownerField', ...
        '%s is missing %s.',label,name);
    v=s.(name);
    assert(isnumeric(v)&&isreal(v)&&isequal(size(v),[n,1]) ...
        &&all(isfinite(v))&&all(v>=0), ...
        'soil:nitrogen:tracer:ownerValue','Invalid %s.%s owner.',label,name);
    v=double(v);
end

function value=toMgG(bottomUp,bdTop)
    value=flip(bottomUp)./bdTop.*1000;
end

function value=fromMgG(topDown,bdTop)
    value=flip(topDown.*bdTop./1000);
end

function [total,tag]=validateDin(rootDin,n)
    assert(isstruct(rootDin)&&isfield(rootDin,'total_mg_g') ...
        &&isfield(rootDin,'tag_mg_g'),'soil:nitrogen:tracer:rootDin', ...
        'rootDin requires total_mg_g and tag_mg_g.');
    total=double(rootDin.total_mg_g);tag=double(rootDin.tag_mg_g);
    assert(isequal(size(total),[n,2])&&isequal(size(tag),[n,2]) ...
        &&all(isfinite([total(:);tag(:)])) ...
        &&all(total(:)>=0)&&all(tag(:)>=0)&&all(tag(:)<=total(:)), ...
        'soil:nitrogen:tracer:rootDin','Invalid post-root total/tag DIN inventories.');
end

function ledger=validateLedger(ledger,n)
    if ~isfield(ledger,'enzyme_cn'),ledger.enzyme_cn=3;end
    assert(isnumeric(ledger.enzyme_cn)&&isreal(ledger.enzyme_cn) ...
        &&isscalar(ledger.enzyme_cn)&&isfinite(ledger.enzyme_cn)&&ledger.enzyme_cn>0, ...
        'soil:nitrogen:tracer:enzymeCN','Accepted enzyme C:N must be positive and finite.');
    assert(isstruct(ledger)&&isfield(ledger,'nitrogen_flux') ...
        &&isstruct(ledger.nitrogen_flux), ...
        'soil:nitrogen:tracer:sourceLedger','Accepted source_hour_ledger is required.');
    required={'Nim_NH4','Nim_NO3','Nmn_MBA','Nmn_MBD', ...
        'first_NH4','first_NO3','second_NH4','second_NO3', ...
        'phi_n','enzyme_production_C','enzyme_turnover_C'};
    for j=1:numel(required)
        getAmount(ledger,required{j},n,strcmp(required{j},'phi_n'));
    end
    assertClose(getAmount(ledger,'Nim_NH4',n,false), ...
        getAmount(ledger,'first_NH4',n,false)+getAmount(ledger,'second_NH4',n,false), ...
        'NH4 uptake stages');
    assertClose(getAmount(ledger,'Nim_NO3',n,false), ...
        getAmount(ledger,'first_NO3',n,false)+getAmount(ledger,'second_NO3',n,false), ...
        'NO3 uptake stages');
end

function amount=getAmount(s,name,n,allowSigned)
    assert(isstruct(s)&&isfield(s,name),'soil:nitrogen:tracer:ledgerField', ...
        'Accepted source ledger is missing %s.',name);
    amount=double(s.(name));
    assert(isreal(amount)&&isequal(size(amount),[n,1])&&all(isfinite(amount)), ...
        'soil:nitrogen:tracer:ledgerValue','Invalid accepted route %s.',name);
    if ~allowSigned
        requireNonnegative(amount,name);
    end
end

function value=getFraction(s,name,n)
    value=getAmount(s,name,n,true);
    assert(all(value>=0&value<=1),'soil:nitrogen:tracer:ledgerFraction', ...
        'Accepted source fraction %s must lie in [0,1].',name);
end

function assertZeroFlux(s,name)
    if ~isfield(s,name),return;end
    value=double(s.(name));
    assert(all(isfinite(value(:)))&&all(abs(value(:))<=1e-13), ...
        'soil:nitrogen:tracer:sorptionGross', ...
        'Linear-equilibrium replay cannot consume nonzero %s.',name);
end

function validateSourceIdentities(nf,n)
    p1=getAmount(nf,'POMdec1',n,false);
    p2=getAmount(nf,'POMdec2',n,false);
    assertClose(p1,getAmount(nf,'POMdec2DOM1',n,false) ...
        +getAmount(nf,'POMdec2MOM1',n,false),'POM1 gross split');
    assertClose(p2,getAmount(nf,'POMdec2DOM2',n,false) ...
        +getAmount(nf,'POMdec2MOM2',n,false),'POM2 gross split');
    assertClose(getAmount(nf,'MOMdec',n,false), ...
        getAmount(nf,'MOM2DOM',n,false),'MOM gross transfer');
    mortality=getAmount(nf,'MBA2DOM',n,false) ...
        +getAmount(nf,'MBA2POM1',n,false)+getAmount(nf,'MBA2POM2',n,false);
    assertClose(getAmount(nf,'MBA_mortality',n,false),mortality, ...
        'MBA mortality gross split');
    pm=mortality+getAmount(nf,'MBA2EP1',n,false) ...
        +getAmount(nf,'MBA2EP2',n,false)+getAmount(nf,'MBA2EM',n,false);
    assertClose(getAmount(nf,'MBA_PM',n,false),pm,'MBA gross outflow');
end

function route=tagRoute(amount,tagDonor,totalDonor,label)
    assert(isequal(size(amount),size(totalDonor)) ...
        &&isequal(size(tagDonor),size(totalDonor)), ...
        'soil:nitrogen:tracer:routeShape','Route shape mismatch for %s.',label);
    requireNonnegative(amount,label);
    validatePair(totalDonor,tagDonor,[label ' donor']);
    tolerance=2e-11.*max(1,abs(totalDonor));
    assert(all(amount<=totalDonor+tolerance,'all'), ...
        'soil:nitrogen:tracer:routeOverdraw','Gross route overdraws %s.',label);
    fraction=zeros(size(totalDonor));positive=totalDonor>0;
    fraction(positive)=tagDonor(positive)./totalDonor(positive);
    assert(~any(amount(~positive)>0), ...
        'soil:nitrogen:tracer:emptyDonor','Positive %s route has an empty donor.',label);
    route=amount.*fraction;
    assert(all(isfinite(route),'all')&&all(route>=0,'all') ...
        &&all(route<=amount,'all'), ...
        'soil:nitrogen:tracer:routeValue','Invalid tagged route for %s.',label);
end

function [direct,surplus]=splitAcceptedMbaRelease(accepted,directEstimate)
% Preserve the accepted total release when two equivalent multiplication
% orders differ by a few floating-point units.  Production scales Nmn_MBA
% after source evaluation, whereas replay reconstructs the direct component
% from the separately scaled DOM2MBA route.  A roundoff-negative remainder
% is therefore physical zero, but a materially negative split is invalid.
    assert(isequal(size(accepted),size(directEstimate)) ...
            &&all(isfinite([accepted(:);directEstimate(:)])) ...
            &&all(accepted(:)>=0)&&all(directEstimate(:)>=0), ...
        'soil:nitrogen:tracer:negativeRoute', ...
        'Invalid accepted/direct MBA release values.');
    direct=directEstimate;
    surplus=accepted-direct;
    scale=max(abs(accepted),abs(direct));
    roundoffTolerance=8.*eps(scale);
    roundoffNegative=surplus<0&surplus>=-roundoffTolerance;
    direct(roundoffNegative)=accepted(roundoffNegative);
    surplus(roundoffNegative)=0;
    assert(all(surplus>=0,'all'), ...
        'soil:nitrogen:tracer:negativeRoute', ...
        'MBA surplus release is negative beyond floating-point roundoff.');
end

function tagged=tagDinWithdrawal(amount,total,tag,label)
    tagged=tagRoute(amount,tag,total,label);
end

function [work,route]=switchPhases(total,tag,a2d,d2a)
    validatePair(total.MBA,tag.MBA,'MBA initial phases');
    validatePair(total.MBD,tag.MBD,'MBD initial phases');
    [aTotal,aTag]=withdrawByPhase(total.MBA,tag.MBA,a2d,'MBA to MBD');
    total.MBA=total.MBA-aTotal;tag.MBA=tag.MBA-aTag;
    total.MBD=total.MBD+aTotal;tag.MBD=tag.MBD+aTag;
    [dTotal,dTag]=withdrawByPhase(total.MBD,tag.MBD,d2a,'MBD to MBA');
    total.MBD=total.MBD-dTotal;tag.MBD=tag.MBD-dTag;
    total.MBA=total.MBA+dTotal;tag.MBA=tag.MBA+dTag;
    validatePair(total.MBA,tag.MBA,'MBA switched phases');
    validatePair(total.MBD,tag.MBD,'MBD switched phases');
    work=struct('total',total,'tag',tag);
    route=struct('total_active_to_dormant',a2d, ...
        'tag_active_to_dormant',sum(aTag,2), ...
        'total_dormant_to_active',d2a, ...
        'tag_dormant_to_active',sum(dTag,2), ...
        'tag_active_to_dormant_by_phase',aTag, ...
        'tag_dormant_to_active_by_phase',dTag);
end

function [totalRoute,tagRouteValue]=withdrawByPhase(total,tag,amount,label)
    donor=sum(total,2);requireNonnegative(amount,label);
    tolerance=2e-11.*max(1,abs(donor));
    assert(all(amount<=donor+tolerance), ...
        'soil:nitrogen:tracer:phaseOverdraw','%s overdraws its phase donor.',label);
    totalRoute=zeros(size(total));positive=donor>0;
    if any(positive)
        totalRoute(positive,:)=total(positive,:).*amount(positive)./donor(positive);
    end
    assert(~any(amount(~positive)>0), ...
        'soil:nitrogen:tracer:emptyPhaseDonor','%s has an empty phase donor.',label);
    tagRouteValue=tagRoute(totalRoute,tag,total,label);
end

function [totalOut,tagOut,allocation]=resizePhaseGroup(totalIn,tagIn,totalTarget,tagTarget,allowNew,label)
    totalTarget=totalTarget(:);tagTarget=tagTarget(:);
    validatePair(totalIn,tagIn,[label ' pre-resize phases']);
    validatePair(totalTarget,tagTarget,[label ' resize target']);
    current=sum(totalIn,2);
    totalOut=zeros(size(totalIn));
    positive=current>0;
    if any(positive)
        totalOut(positive,:)=bsxfun(@times,totalIn(positive,:), ...
            totalTarget(positive)./current(positive));
    end
    empty=~positive&totalTarget>0;
    assert(allowNew||~any(empty), ...
        'soil:nitrogen:tracer:phaseAppearance', ...
        '%s appeared without a donor phase.',label);
    totalOut(empty,2)=totalTarget(empty);
    % Change: 2026-09-22 Wenjie Liu -- preserve aggregate gross attribution
    % while respecting each immutable total phase's receiving capacity.
    [tagOut,allocation]=soil.nitrogen.tracer.allocatePhaseLabels(totalIn,tagIn,totalOut,tagTarget);
    validatePair(totalOut,tagOut,[label ' resized phases']);
end

function [production,turnover]=enzymePhaseAmounts(ledger,n,before,after,productionTotal,turnoverTotal)
% Applied C routes come from the same cell-phase loop as the total update.
% Never substitute N biomass shares: production uses the producer's C shares.
    production=getPhaseCarbonAmount(ledger,'enzyme_production_by_phase_C',n)./ledger.enzyme_cn;
    turnover=getPhaseCarbonAmount(ledger,'enzyme_turnover_by_phase_C',n)./ledger.enzyme_cn;
    assertClose(sum(production,2),productionTotal,'intracellular-enzyme production');
    assertClose(sum(turnover,2),turnoverTotal,'intracellular-enzyme turnover');
    assertClose(before+production-turnover,after,'intracellular-enzyme phase endpoint');
end

function amount=getPhaseCarbonAmount(ledger,name,n)
    assert(isfield(ledger,name),'soil:nitrogen:tracer:ledgerField', ...
        'Accepted source ledger is missing applied phase route %s.',name);
    amount=ledger.(name);
    assert(isnumeric(amount)&&isreal(amount)&&isequal(size(amount),[n,3]) ...
        &&all(isfinite(amount),'all')&&all(amount>=0,'all'), ...
        'soil:nitrogen:tracer:ledgerValue', ...
        'Applied phase route %s must be finite nonnegative NL-by-3 mg C/g.',name);
    amount=double(amount);
end

function [totalOut,tagOut]=pairedUpdate(totalIn,tagIn,totalDelta,tagDelta,label)
    totalOut=totalIn+totalDelta;tagOut=tagIn+tagDelta;
    validatePair(totalOut,tagOut,[label ' gross update']);
end

function validatePoolBounds(total,tag,label)
    pools={'POM1','POM2','MOM','DOM','QOM','MBA','MBD','ENZP1','ENZP2','ENZM','ENZN'};
    for j=1:numel(pools)
        key=pools{j};validatePair(total.(key),tag.(key),[label ' ' key]);
    end
    for key={'MBA','MBD','ENZN'}
        k=key{1};validatePair(total.phase.(k),tag.phase.(k),[label ' ' k ' phases']);
    end
end

function validatePair(total,tag,label)
    assert(isequal(size(total),size(tag))&&all(isfinite([total(:);tag(:)])) ...
        &&all(total(:)>=0)&&all(tag(:)>=0), ...
        'soil:nitrogen:tracer:tagValue','%s contains invalid total/tag values.',label);
    tolerance=5e-11.*max(1,abs(total));
    assert(all(tag<=total+tolerance,'all'), ...
        'soil:nitrogen:tracer:tagBound','%s fertilizer tag exceeds total N.',label);
end

function requireNonnegative(value,label)
    tolerance=5e-12.*max(1,abs(value));
    assert(all(isfinite(value),'all')&&all(value>=-tolerance,'all'), ...
        'soil:nitrogen:tracer:negativeRoute','%s is negative.',label);
end

function assertClose(actual,expected,label)
    assert(isequal(size(actual),size(expected))&&all(isfinite([actual(:);expected(:)])), ...
        'soil:nitrogen:tracer:comparison','Invalid values in %s comparison.',label);
    tolerance=5e-10.*max(1,max(abs(actual),abs(expected)));
    assert(all(abs(actual-expected)<=tolerance,'all'), ...
        'soil:nitrogen:tracer:totalParity', ...
        'Accepted total endpoint is inconsistent with %s.',label);
end

function value=poolInventory(s)
    value=s.POM1+s.POM2+s.MOM+s.DOM+s.QOM+s.MBA+s.MBD ...
        +s.ENZP1+s.ENZP2+s.ENZM+s.ENZN;
end

function out=writeOwners(template,q,mb,mbd,enzn,bdTop)
    out=template;o=out.organic_n;
    map={ ...
        'PON1','POM1';'PON2','POM2';'MON','MOM'; ...
        'MB_N','MBA';'BD_N','MBD';'ENZP1_N','ENZP1'; ...
        'ENZP2_N','ENZP2';'EM_N','ENZM';'ENZN_N','ENZN'};
    for j=1:size(map,1),o.(map{j,1})=fromMgG(q.(map{j,2}),bdTop);end
    o.PON=o.PON1+o.PON2;o.EP_N=o.ENZP1_N+o.ENZP2_N;
    phases={'mobile','attached','retained'};
    matrices={mb,mbd,enzn};groups={'MB','BD','ENZN'};
    for g=1:3
        for p=1:3
            o.([groups{g} '_' phases{p} '_N'])=fromMgG(matrices{g}(:,p),bdTop);
        end
    end
    if isfield(o,'QON'),o.QON(:)=0;end
    out.organic_n=o;
    out.don_total_g_cm3=fromMgG(q.DOM,bdTop);
    if strcmp(ownerMode(template),'mend_dynamic')
        out.don_sorbed_g_cm3=fromMgG(q.QOM,bdTop);
        out.don_total_g_cm3=fromMgG(q.DOM+q.QOM,bdTop);
    end
end

function mode=ownerMode(bundle)
    mode='linear_equilibrium';
    if isfield(bundle,'sorption_model'),mode=char(string(bundle.sorption_model));end
end

function [adsTag,desTag,freeTagNext,sorbedTagNext,roundoff]=sorptionLabels( ...
        free,freeTag,sorbed,sorbedTag,ads,des,freeNext,sorbedNext)
% Attribute two accepted opposing transfers without overdraw or net inference.
% Fractions solve the same two-box implicit mixing contract as frozen DON
% transport: (old+incoming)*fraction - incoming*other_fraction = old label.
% This also permits accepted gross recycling larger than a starting donor.
    validatePair(free,freeTag,'free DON after accepted non-sorption gross routes');
    validatePair(sorbed,sorbedTag,'sorbed DON before exchange');
    validatePair(freeNext,zeros(size(free)),'accepted free DON after exchange');
    validatePair(sorbedNext,zeros(size(free)),'accepted sorbed DON after exchange');
    assertClose(free-ads+des,freeNext,'accepted free DON sorption endpoint');
    assertClose(sorbed+ads-des,sorbedNext,'accepted sorbed DON sorption endpoint');
    % The source map can emit signed adsorption above Qmax. Normalize each
    % signed accepted route to its physical direction, preserving its net
    % and its original signed public route label; never clamp total fluxes.
    forward=max(ads,0)+max(-des,0);backward=max(des,0)+max(-ads,0);
    adsTag=zeros(size(free));desTag=adsTag;
    freeTagNext=adsTag;sorbedTagNext=adsTag;
    for i=1:numel(free)
        matrix=[free(i)+backward(i),-backward(i);-forward(i),sorbed(i)+forward(i)];
        rhs=[freeTag(i);sorbedTag(i)];
        if free(i)+sorbed(i)==0
            assert(ads(i)==0&&des(i)==0,'soil:nitrogen:tracer:emptySorption', ...
                'An empty DOM pair cannot carry accepted sorption exchange.');
            continue
        end
        empty=all(matrix==0,2);
        matrix(empty,empty)=eye(nnz(empty));
        fraction=matrix\rhs;
        assert(all(isfinite(fraction))&&all(fraction>=-512*eps&fraction<=1+512*eps), ...
            'soil:nitrogen:tracer:sorptionFraction','Invalid bounded sorption label fraction.');
        % Only numerical roundoff may cross the exact mixing bounds.
        fraction=min(max(fraction,0),1);
        adsTag(i)=max(ads(i),0)*fraction(1)-max(-ads(i),0)*fraction(2);
        desTag(i)=max(des(i),0)*fraction(2)-max(-des(i),0)*fraction(1);
        % Stable at exact depletion: subtracting opposing gross labels can
        % otherwise create a negative ULP despite an accepted zero stock.
        freeTagNext(i)=freeNext(i)*fraction(1);
        sorbedTagNext(i)=sorbedNext(i)*fraction(2);
    end
    freeGross=freeTag-adsTag+desTag;sorbedGross=sorbedTag+adsTag-desTag;
    tol=2048*eps(max([free;sorbed;freeNext;sorbedNext;realmin]));
    assert(all(abs(freeTagNext-freeGross)<=tol) ...
        &&all(abs(sorbedTagNext-sorbedGross)<=tol), ...
        'soil:nitrogen:tracer:sorptionRoundoff', ...
        'Stable sorption endpoints disagree with accepted gross label transfers.');
    roundoff=(freeTagNext-freeGross)+(sorbedTagNext-sorbedGross);
end
