function [organic,dom,diagnostic]=mendOrganicNitrogenCandidate( ...
        organic,dom,carbonBefore,carbonAfter,substepLedger,theta,bd, ...
        parameters,domParameters)
%MENDORGANICNITROGENCANDIDATE Advance organic N on accepted MEND RK stages.
%   Organic and DOM states are bottom-up; MEND carbon pools/fluxes are
%   top-down [mg C g-1]. SUBSTEPLEDGER contains only accepted carbon RK4
%   substeps. THETA and BD are bottom-up layer profiles; DOMPARAMETERS owns
%   the bottom-up DON Kd profile. The reaction repartitions total DON at
%   each RK stage and consumes only its aqueous fraction. Enzyme production
%   and turnover use fixed C:N parameters. Total
%   organic N remains conserved without a second accepted DOM owner; DIN
%   requests are handled separately by the soil-N solver.
%   Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-12 Codex -- derive net mineral DIN immobilization or
%   mineralization from the accepted MEND microbial C:N envelope.
% Change: 2026-09-12 Codex -- replace beginning-of-step flux replay with
%   synchronous organic-N integration on every accepted carbon RK4 stage.
% Change: 2026-09-13 Codex -- fix enzyme C:N and aqueous DON uptake on the
%   same frozen sorption interval as the accepted carbon RK trajectory.
% Change: 2026-09-13 Codex -- release (1-YgN) of DON uptake to NH4 and
%   calculate the MBA DIN deficit separately from MBD C:N correction.
% Change: 2026-09-13 Codex -- repartition stage-total DON under Kd equilibrium
%   in lockstep with the carbon stage's aqueous DOC donor.
% Change: 2026-09-14 Codex -- carry dormant/active N through matching physical
%   phases on each accepted carbon substep before C:N release adjustments.
    bd=bd(:);n=numel(bd);bdTop=flip(bd);
    required={'PON1','PON2','MON','MB_N','BD_N','ENZP1_N','ENZP2_N','EM_N'};
    for i=1:numel(required)
        k=required{i}; assert(isfield(organic,k), ...
            'soil:microbe:mendNitrogenState','Missing MEND N pool %s.',k);
    end
    assert(isfield(dom,'don_total_g_cm3')&&numel(dom.don_total_g_cm3)==n, ...
        'soil:microbe:mendNitrogenState','MEND requires total DON ownership.');
    assert(numel(theta)==n&&all(isfinite(theta(:)))&&all(theta(:)>=0) ...
        &&isfield(domParameters,'kd_don_cm3_g'), ...
        'soil:microbe:mendNitrogenParameters', ...
        'MEND requires the bottom-up water and DON sorption profiles.');
    assert(isfield(parameters,'CN_ENZP')&&numel(parameters.CN_ENZP)==2 ...
        &&all(isfinite(parameters.CN_ENZP(:))) ...
        &&all(parameters.CN_ENZP(:)>0) ...
        &&isfield(parameters,'CN_ENZM') ...
        &&isscalar(parameters.CN_ENZM) ...
        &&isfinite(parameters.CN_ENZM)&&parameters.CN_ENZM>0, ...
        'soil:microbe:mendNitrogenParameters', ...
        'MEND enzyme C:N parameters must be finite and positive.');
    assert(isfield(parameters,'CN_MB_min') ...
        &&isfield(parameters,'CN_MB_max') ...
        &&isfield(parameters,'wdorm'), ...
        'soil:microbe:mendNitrogenParameters', ...
        'MEND microbial C:N bounds and NUE exponent are required.');
    cnMin=expand(parameters.CN_MB_min,n);
    cnMax=expand(parameters.CN_MB_max,n);
    assert(all(cnMin>0)&all(cnMax>cnMin), ...
        'soil:microbe:mendNitrogenParameters', ...
        'MEND microbial C:N bounds must satisfy 0<CN_min<CN_max.');
    donPhase=soil.dom.equilibrate(ones(n,1),theta(:), ...
        bd,domParameters.kd_don_cm3_g);
    aqueousFraction=flip(theta(:).*donPhase.aqueous_g_cm3_water(:));
    doc=flip(dom.doc_total_g_cm3(:))./bdTop*1000;
    N=struct(); for i=1:numel(required), k=required{i}; N.(k)=flip(organic.(k)(:))./bdTop*1000; end
    N.DON=flip(dom.don_total_g_cm3(:))./bdTop*1000;
    initialOrganicN=N.DON;
    for i=1:numel(required), initialOrganicN=initialOrganicN+N.(required{i}); end
    assert(isequal(size(carbonBefore.POM1),size(carbonAfter.POM1)), ...
        'soil:microbe:mendNitrogenState','MEND carbon endpoints must align.');
    assertEnzymeStoichiometry(carbonBefore,N,parameters);
    [N,donUptake,nueMineralization,phaseLedger]=advanceOnCarbonLedger( ...
        N,substepLedger,doc,aqueousFraction,parameters);
    organic=soil.microbe.carryReactionPhases(organic,phaseLedger,'N',bd);
    don=N.DON;
    dinRequest=max(carbonAfter.MB(:)./cnMax-N.MB_N,0);
    mbaRelease=max(N.MB_N-carbonAfter.MB(:)./cnMin,0);
    mbdRelease=max(N.BD_N-carbonAfter.BD(:)./cnMin,0);
    N.MB_N=N.MB_N-mbaRelease;
    N.BD_N=N.BD_N-mbdRelease;
    dinRelease=nueMineralization+mbaRelease+mbdRelease;
    donAssimilated=donUptake-nueMineralization;
    assert(all(donAssimilated>=-1e-12), ...
        'soil:microbe:mendNitrogenClosure', ...
        'MEND intrinsic NUE mineralized more N than DON uptake.');
    vals=[]; for i=1:numel(required), vals=[vals; N.(required{i})(:)]; end %#ok<AGROW>
    assert(all(isfinite(vals))&&all(vals>=-1e-10), ...
        'soil:microbe:mendNitrogenOverdraw','MEND donor-N transfer overdraw.');
    assert(all(don>=-1e-12), ...
        'soil:microbe:mendNitrogenOverdraw','MEND DON uptake exceeded the donor pool.');
    assertNitrogenState(N,required);
    assertEnzymeStoichiometry(carbonAfter,N,parameters);
    finalOrganicN=don+dinRelease;
    for i=1:numel(required), finalOrganicN=finalOrganicN+N.(required{i}); end
    nitrogenResidual=finalOrganicN-initialOrganicN;
    assert(max(abs(nitrogenResidual))<=2e-10, ...
        'soil:microbe:mendNitrogenClosure', ...
        'MEND donor transfer and mineralization failed N closure.');
    for i=1:numel(required), k=required{i}; organic.(k)=flip(N.(k).*bdTop/1000); end
    organic=soil.microbe.scaleReactionPhases(organic,'N');
    organic.PON=organic.PON1+organic.PON2; organic.EP_N=organic.ENZP1_N+organic.ENZP2_N;
    dom.don_total_g_cm3=flip(don.*bdTop/1000);
    organic.qon_ownership='dom_total_sole_owner';
    diagnostic=struct('don_total_g_cm3',dom.don_total_g_cm3, ...
        'nitrogen_residual_mg_g',nitrogenResidual, ...
        'din_request_mg_g',dinRequest,'din_release_mg_g',dinRelease, ...
        'din_release_mba_mg_g',mbaRelease, ...
        'din_release_mbd_mg_g',mbdRelease, ...
        'nue_mineralization_mg_g',nueMineralization, ...
        'don_uptake_mg_g',donUptake, ...
        'don_assimilated_mg_g',donAssimilated, ...
        'source_model','MEND','same_stage',true, ...
        'stage_count',4*numel(substepLedger), ...
        'substep_count',numel(substepLedger));
end

function [state,donUptake,nueMineralization,phaseLedger]=advanceOnCarbonLedger( ...
        state,ledger,totalDoc,aqueousFraction,parameters)
% Advance N with the same RK4 abscissae and weights as accepted carbon.
% Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-13 Codex -- repartition stage-total DON before donor uptake.
% Change: 2026-09-13 Codex -- integrate source NUE mineralization at each
%   accepted carbon stage and return it as a mineral-N ledger amount.
    assert(iscell(ledger)&&~isempty(ledger), ...
        'soil:microbe:mendNitrogenLedger', ...
        'MEND organic N requires a nonempty accepted carbon substep ledger.');
    first=ledger{1}.stages(1).state;
    assert(all(isfinite(totalDoc-first.DOC(:))) ...
        &&all(totalDoc-first.DOC(:)>=-1e-12), ...
        'soil:microbe:mendNitrogenLedger', ...
        'Accepted total DOC must contain the first-stage aqueous DOC.');
    fields=fieldnames(state);
    assertNitrogenState(state,fields);
    donUptake=zeros(size(totalDoc));
    nueMineralization=zeros(size(totalDoc));
    phaseLedger=cell(size(ledger));
    for j=1:numel(ledger)
        item=ledger{j};
        assert(isfield(item,'stages')&&numel(item.stages)==4 ...
            &&isfield(item,'dt_s')&&isscalar(item.dt_s)&&item.dt_s>0, ...
            'soil:microbe:mendNitrogenLedger', ...
            'Each accepted MEND carbon substep must contain four RK stages.');
        s=item.stages; h=item.dt_s/3600;
        beforePhase=struct('MB',state.MB_N,'BD',state.BD_N);
        [k1,u1,m1,p1]=nitrogenDerivative( ...
            s(1).state,state,s(1).flux,aqueousFraction,parameters);
        n2=addScaledChecked(state,k1,h/2,fields);
        [k2,u2,m2,p2]=nitrogenDerivative( ...
            s(2).state,n2,s(2).flux,aqueousFraction,parameters);
        n3=addScaledChecked(state,k2,h/2,fields);
        [k3,u3,m3,p3]=nitrogenDerivative( ...
            s(3).state,n3,s(3).flux,aqueousFraction,parameters);
        n4=addScaledChecked(state,k3,h,fields);
        [k4,u4,m4,p4]=nitrogenDerivative( ...
            s(4).state,n4,s(4).flux,aqueousFraction,parameters);
        for i=1:numel(fields)
            key=fields{i};
            state.(key)=state.(key)+h/6.*(k1.(key)+2.*k2.(key) ...
                +2.*k3.(key)+k4.(key));
        end
        assertNitrogenState(state,fields);
        donUptake=donUptake+h/6.*(u1+2.*u2+2.*u3+u4);
        nueMineralization=nueMineralization+h/6.*(m1+2.*m2+2.*m3+m4);
        phaseLedger{j}=struct('before',beforePhase, ...
            'after',struct('MB',state.MB_N,'BD',state.BD_N), ...
            'a_to_d',h/6.*(p1.a_to_d+2*p2.a_to_d+2*p3.a_to_d+p4.a_to_d), ...
            'd_to_a',h/6.*(p1.d_to_a+2*p2.d_to_a+2*p3.d_to_a+p4.d_to_a), ...
            'dormant_loss',zeros(size(state.BD_N)));
    end
end

function [d,donUptake,nueMineralization,phaseSwitch]=nitrogenDerivative( ...
        c,n,f,aqueousFraction,parameters)
% Conservative MEND organic-N tendency [mg N g-1 h-1].
% Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-13 Codex -- apply aqueous DON:N and fixed enzyme C:N.
% Change: 2026-09-13 Codex -- partition DON uptake by intrinsic NUE.
% Change: 2026-09-13 Codex -- use the current total DON's aqueous fraction.
    pom1=transfer(f.POM1_decomp,n.PON1,c.POM1);
    pom2=transfer(f.POM2_decomp,n.PON2,c.POM2);
    mom=transfer(f.MOM_to_DOM,n.MON,c.MOC);
    mortality=transfer(f.MBA_mortality,n.MB_N,c.MB);
    mbaDom=transfer(f.MBA_to_DOM,n.MB_N,c.MB);
    mbaPom1=transfer(f.MBA_to_POM1,n.MB_N,c.MB);
    mbaPom2=transfer(f.MBA_to_POM2,n.MB_N,c.MB);
    aqueousDon=aqueousFraction.*max(n.DON,0);
    donUptake=transfer(f.DOM_to_MBA,aqueousDon,c.DOC);
    yNg=soil.microbe.mendNitrogenUseEfficiency( ...
        c.MB,n.MB_N,parameters.CN_MB_min, ...
        parameters.CN_MB_max,parameters.wdorm);
    nueMineralization=(1-yNg).*donUptake;
    enzP1Prod=f.MBA_to_ENZP1./parameters.CN_ENZP(1);
    enzP2Prod=f.MBA_to_ENZP2./parameters.CN_ENZP(2);
    enzMProd=f.MBA_to_ENZM./parameters.CN_ENZM;
    enzP1Loss=f.ENZP1_to_DOM./parameters.CN_ENZP(1);
    enzP2Loss=f.ENZP2_to_DOM./parameters.CN_ENZP(2);
    enzMLoss=f.ENZM_to_DOM./parameters.CN_ENZM;
    toDormant=transfer(f.MBA_to_MBD,n.MB_N,c.MB);
    toActive=transfer(f.MBD_to_MBA,n.BD_N,c.BD);
    phaseSwitch=struct('a_to_d',toDormant,'d_to_a',toActive);
    pom1Dom=fraction(f.POM1_to_DOM,f.POM1_decomp).*pom1;
    pom2Dom=fraction(f.POM2_to_DOM,f.POM2_decomp).*pom2;
    d.PON1=-pom1+mbaPom1;
    d.PON2=-pom2+mbaPom2;
    d.MON=pom1-pom1Dom+pom2-pom2Dom-mom;
    d.MB_N=donUptake-nueMineralization ...
        -mortality-enzP1Prod-enzP2Prod-enzMProd ...
        -toDormant+toActive;
    d.BD_N=toDormant-toActive;
    d.ENZP1_N=enzP1Prod-enzP1Loss;
    d.ENZP2_N=enzP2Prod-enzP2Loss;
    d.EM_N=enzMProd-enzMLoss;
    d.DON=pom1Dom+pom2Dom+mom+mbaDom+enzP1Loss+enzP2Loss ...
        +enzMLoss-donUptake;
    residual=zeros(size(n.DON)); names=fieldnames(d);
    for i=1:numel(names), residual=residual+d.(names{i}); end
    assert(max(abs(residual+nueMineralization))<=2e-13, ...
        'soil:microbe:mendNitrogenClosure', ...
        'MEND RK-stage organic plus mineral-N tendency failed closure.');
end

function q=transfer(carbonFlux,donorN,donorC)
    q=zeros(size(carbonFlux)); ok=donorC>0;
    q(ok)=carbonFlux(ok).*donorN(ok)./donorC(ok);
end

function value=fraction(numerator,denominator)
    value=zeros(size(numerator)); ok=denominator>0;
    value(ok)=numerator(ok)./denominator(ok);
end

function out=addScaledChecked(state,k,scale,fields)
    out=state;
    for i=1:numel(fields)
        key=fields{i}; out.(key)=state.(key)+scale.*k.(key);
    end
    assertNitrogenState(out,fields);
end

function assertNitrogenState(state,fields)
    for i=1:numel(fields)
        value=state.(fields{i});
        assert(all(isfinite(value))&&all(value>=-1e-12), ...
            'soil:microbe:mendNitrogenStageNegative', ...
            'MEND RK-stage nitrogen pool %s became invalid.',fields{i});
    end
end

function assertEnzymeStoichiometry(c,n,p)
% Existing enzyme pools must already follow the fixed MEND C:N contract.
% Editor: Jizixing Zhu; Date: 2026-09-13.
    ratio={p.CN_ENZP(1),p.CN_ENZP(2),p.CN_ENZM};
    carbon={'ENZP1','ENZP2','EM'};
    nitrogen={'ENZP1_N','ENZP2_N','EM_N'};
    for i=1:3
        difference=n.(nitrogen{i})(:)-c.(carbon{i})(:)./ratio{i};
        assert(all(abs(difference)<=1e-10.*max(1,c.(carbon{i})(:))), ...
            'soil:microbe:mendEnzymeStoichiometry', ...
            'MEND enzyme %s violates fixed C:N; rebuild the checkpoint.', ...
            carbon{i});
    end
end

function value=expand(value,n)
    value=value(:);
    if isscalar(value), value=repmat(value,n,1); end
    assert(numel(value)==n&&all(isfinite(value)), ...
        'soil:microbe:mendNitrogenParameters', ...
        'MEND microbial C:N parameter must be scalar or NL-by-1.');
end
