function [out, nn, newdin, ledger] = mendSourceHour(c,n,din,p,addC,addN,ihr,dt_s,coupling)
%MENDSOURCEHOUR Source MEND discrete hourly organic C/N and DIN exchange.
% Inputs/outputs are structs of aligned column pools [mg element/g dry soil]
% with source archive field names. Each entry is an independent box; no grid
% reversal. Parameters are scalar or matching columns, already evaluated at
% the hourly environment. Additions are hourly amounts. CO2 is cumulative.
% Source: MOD_MEND.F90 (92323c7), subMEND MM=0, fDOM=1, iHR=0/1,
% organic balances and aggregate intracellular enzyme ledger. Enzyme CN=3,
% microbial CN bounds 2--14, source plant uptake=0. DIN output precedes
% inorganic chemistry. No transport, oxygen/root allocation or state commit.
% Undefined inputs, non-hour durations and negative proposals throw errors;
% no silent clipping or redistribution. See MEND_MATLAB_HOURLY_KERNEL_CONTRACT.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-09-14 Jizixing Zhu -- optional root-reserved total DIN inventory
%   caps first uptake independently of free-substrate kinetics. Default remains
%   exact source. Coupled mode omits source plant competition denominator terms.
% Change: 2026-10-09 Wenjie Liu -- bounded substrate C transfers carry a shared
%   donor N budget, with exact depletion and remaining-first balances. DON
%   must be nonnegative; dynamic sorption gross-recycling formulas are retained.
    assert(isscalar(dt_s)&&isfinite(dt_s)&&dt_s==3600, ...
        'soil:microbe:sourceHourClock','Source discrete map requires exactly 3600 s.');
    assert(isscalar(ihr)&&ismember(ihr,[0 1]), ...
        'soil:microbe:sourceHourMode','Only source iHR=0/1 is supported.');
    pools={'POMO','POMH','MOM','QOM','DOMT','MBA','MBD', ...
        'ENZP1','ENZP2','ENZM','ENZNmt'};
    assert(isfield(c,'MBA_C')&&iscolumn(c.MBA_C), ...
        'soil:microbe:sourceHourShape','Source states must be nonempty columns.');
    count=numel(c.MBA_C);
    coupled=nargin>=9&&~isempty(coupling);
    inventory=din;
    if coupled
        validateField(coupling,'NH4',count);validateField(coupling,'NO3',count);
        assert(all(coupling.NH4>=din.NH4)&all(coupling.NO3>=din.NO3), ...
            'soil:microbe:sourceHourInventory','Free DIN exceeds post-root total inventory.');
        inventory=coupling;
    end
    for j=1:numel(pools)
        validateField(c,[pools{j} '_C'],count);
        validateField(n,[pools{j} '_N'],count);
    end
    validateField(c,'CO2',count);
    for key={'NH4','NO3','NO2'},validateField(din,key{1},count);end
    for key={'POMinp1','POMinp2','DOMinp'}
        validateField(addC,key{1},count);validateField(addN,key{1},count);
    end
    for key={'NH4_dep','NO3_dep'},validateField(addN,key{1},count);end
    positive={'KP1','KP2','KM','KD','Qmax','Yg','KsNH4_MB','KsNO3_MB'};
    if ~coupled,positive=[positive {'KsNH4_VG','KsNO3_VG'}];end
    nonnegative={'VP1','VP2','VM','Vg','Vm','Vm_dorm','Kads','Kdes', ...
        'rMORT','pEP','pEM','rEP1','rEP2','rEM','VmA2D','VmD2A', ...
        'wdorm','VNup_MB','VNup_VG','fpENZN','fD','gD'};
    for key=[positive nonnegative]
        name=key{1};validateField(p,name,count,true);
        if isscalar(p.(name)),p.(name)=repmat(double(p.(name)),count,1);end
    end
    for key=positive
        assert(all(p.(key{1})>0),'soil:microbe:sourceHourParameter', ...
            'Source parameter %s must be positive.',key{1});
    end
    assert(all(p.Yg<=1)&&all(p.fD<=1)&&all(p.gD<=1)&&all(p.VNup_VG==0), ...
        'soil:microbe:sourceHourParameter','Invalid fraction or enabled source plant uptake.');
    assert(all(abs(n.ENZNmt_N-c.ENZNmt_C./3)<=1e-12), ...
        'soil:microbe:sourceHourEnzymeCN','Intracellular enzyme C:N must be 3.');
    assert(all(c.POMO_C+c.POMH_C>0), ...
        'soil:microbe:sourceHourPOM','Source enzyme allocation requires positive total POM.');
    ratio=struct();
    for j=1:numel(pools)-1
        name=pools{j};a=c.([name '_C']);b=n.([name '_N']);
        assert(all(a>0|b==0),'soil:microbe:sourceHourDonorCN', ...
            'N-only organic donor %s has undefined source C:N.',name);
        r=zeros(size(a));hasCarbon=a>0;r(hasCarbon)=b(hasCarbon)./a(hasCarbon);
        ratio.(name)=r;
    end

    f=struct();
    for j=1:2
        if j==1,pool='POMO_C';else,pool='POMH_C';end
        suffix=num2str(j);
        dec=mm(p.(['VP' suffix]),c.(pool),c.(['ENZP' suffix '_C']),p.(['KP' suffix]));
        f.(['POMdec' suffix])=dec;
        f.(['POMdec2DOM' suffix])=p.fD.*dec;
        f.(['POMdec2MOM' suffix])=dec-f.(['POMdec2DOM' suffix]);
    end
    f.MOMdec=mm(p.VM,c.MOM_C,c.ENZM_C,p.KM);f.MOM2DOM=f.MOMdec;
    d=c.DOMT_C;b=c.MBA_C;q=c.QOM_C;
    f.DOM2MBA=mm((p.Vg+p.Vm)./p.Yg,d,b,p.KD);
    if ihr==0
        f.CO2_growth=mm(p.Vg.*(1./p.Yg-1),d,b,p.KD);
        f.CO2_maintn_MBA=mm(p.Vm.*(1./p.Yg-1),d,b,p.KD);
    else
        rate=p.Vg+p.Vm;share=zeros(size(rate));active=rate>0;
        share(active)=p.Vm(active)./rate(active);
        f.CO2_growth=f.DOM2MBA.*(1-p.Yg).*(1-share);
        f.CO2_maintn_MBA=f.DOM2MBA.*(1-p.Yg).*share;
    end
    f.CO2_maintn_MBD=p.Vm_dorm.*c.MBD_C;
    available=d-f.DOM2MBA;
    ads=p.Kads.*available.*(1-q./p.Qmax);des=p.Kdes.*q./p.Qmax;
    limitDes=des>q+ads;des(limitDes)=q(limitDes)+ads(limitDes);
    limitAds=~limitDes&(ads>available+des);
    ads(limitAds)=available(limitAds)+des(limitAds);
    f.DOM2QOM=ads;f.QOM2DOM=des;f.DOM2QOMnet=ads-des;
    f.MBA_mortality=p.rMORT.*b;f.MBA2DOM=p.gD.*f.MBA_mortality;
    f.MBA2POM1=.1.*(1-p.gD).*f.MBA_mortality;
    f.MBA2POM2=.9.*(1-p.gD).*f.MBA_mortality;
    share=c.POMO_C./(c.POMO_C+c.POMH_C);
    f.MBA2EP1=share.*p.pEP.*p.Vm.*b;
    f.MBA2EP2=(1-share).*p.pEP.*p.Vm.*b;
    f.MBA2EM=p.pEM.*p.Vm.*b;
    f.EP2DOM1=p.rEP1.*c.ENZP1_C;f.EP2DOM2=p.rEP2.*c.ENZP2_C;
    f.EM2DOM=p.rEM.*c.ENZM_C;
    phi=d./(d+p.KD);
    f.MBA2MBD=(1-phi).*p.VmA2D.*b;f.MBD2MBA=phi.*p.VmD2A.*c.MBD_C;
    f.MBA_PM=f.MBA_mortality+f.MBA2EP1+f.MBA2EP2+f.MBA2EM;
    mb=b+f.DOM2MBA-f.CO2_growth-f.CO2_maintn_MBA-f.MBA_PM-f.MBA2MBD+f.MBD2MBA;
    f.MBA2MBD=f.MBA2MBD+min(mb,0); % Source signed compensation, not a mass clamp.
    out=c;
    out.POMO_C=(c.POMO_C-f.POMdec1)+addC.POMinp1+f.MBA2POM1;
    out.POMH_C=(c.POMH_C-f.POMdec2)+addC.POMinp2+f.MBA2POM2;
    out.MOM_C=(c.MOM_C-f.MOMdec)+f.POMdec2MOM1+f.POMdec2MOM2;
    out.QOM_C=q+f.DOM2QOMnet;
    out.DOMT_C=(available-f.DOM2QOMnet)+addC.DOMinp ...
        +f.POMdec2DOM1+f.POMdec2DOM2+f.MOMdec ...
        +f.MBA2DOM+f.EP2DOM1+f.EP2DOM2+f.EM2DOM;
    out.MBA_C=max(mb,0);
    out.MBD_C=c.MBD_C-f.CO2_maintn_MBD+f.MBA2MBD-f.MBD2MBA;
    out.ENZP1_C=c.ENZP1_C+f.MBA2EP1-f.EP2DOM1;
    out.ENZP2_C=c.ENZP2_C+f.MBA2EP2-f.EP2DOM2;
    out.ENZM_C=c.ENZM_C+f.MBA2EM-f.EM2DOM;
    out.CO2=c.CO2+f.CO2_growth+f.CO2_maintn_MBA+f.CO2_maintn_MBD;

    nf=struct();
    for j=1:2
        if j==1,pool='POMO';else,pool='POMH';end
        suffix=num2str(j);key=['POMdec' suffix];
        nf.(key)=boundedDonorN(f.(key),c.([pool '_C']),n.([pool '_N']));
        nf.(['POMdec2DOM' suffix])=p.fD.*nf.(key);
        nf.(['POMdec2MOM' suffix])=nf.(key)-nf.(['POMdec2DOM' suffix]);
    end
    nf.MOMdec=boundedDonorN(f.MOMdec,c.MOM_C,n.MOM_N);nf.MOM2DOM=nf.MOMdec;
    nf.DOM2MBA=boundedDonorN(f.DOM2MBA,d,n.DOMT_N);
    donor={'QOM2DOM','QOM';'DOM2QOM','DOMT';'MBA2DOM','MBA'; ...
        'MBA2POM1','MBA';'MBA2POM2','MBA';'MBA_mortality','MBA'; ...
        'MBA2MBD','MBA';'MBD2MBA','MBD'};
    for j=1:size(donor,1),nf.(donor{j,1})=f.(donor{j,1}).*ratio.(donor{j,2});end
    % Uptake and adsorption share the same pre-step DOM donor. For adsorption
    % within the residual C stock, allocate its N from that SAME remainder.
    % Gross recycling can borrow C from contemporaneous desorption; retain
    % the source donor-ratio formula there rather than impose a new N limiter.
    domRemainingN=n.DOMT_N-nf.DOM2MBA;
    boundedAds=f.DOM2QOM>=0&f.DOM2QOM<=available;
    nf.DOM2QOM(boundedAds)=boundedDonorN(f.DOM2QOM(boundedAds), ...
        available(boundedAds),domRemainingN(boundedAds));
    boundedDes=f.QOM2DOM>=0&f.QOM2DOM<=q;
    nf.QOM2DOM(boundedDes)=boundedDonorN(f.QOM2DOM(boundedDes), ...
        q(boundedDes),n.QOM_N(boundedDes));
    domRemainingN=domRemainingN-nf.DOM2QOM+nf.QOM2DOM;
    for key={'MBA2EP1','MBA2EP2','MBA2EM','EP2DOM1','EP2DOM2','EM2DOM'}
        nf.(key{1})=f.(key{1})./3;
    end
    nf.DOM2QOMnet=nf.DOM2QOM-nf.QOM2DOM;
    nf.MBA_PM=nf.MBA_mortality+nf.MBA2EP1+nf.MBA2EP2+nf.MBA2EM;
    cn=2.*ones(size(b));hasN=n.MBA_N>0;cn(hasN)=b(hasN)./n.MBA_N(hasN);
    cn(~hasN&b>0)=Inf;
    phiN=((min(14,max(2,cn))-2)./12).^p.wdorm;
    releaseMB=(1-phiN).*nf.DOM2MBA;
    nh4=din.NH4;no3=din.NO3;
    denominator=1+nh4./p.KsNH4_MB+no3./p.KsNO3_MB;
    if ~coupled,denominator=denominator+nh4./p.KsNH4_VG+no3./p.KsNO3_VG;end
    split=zeros(size(nh4));hasDIN=nh4+no3>0;split(hasDIN)=nh4(hasDIN)./(nh4(hasDIN)+no3(hasDIN));
    uptake4=p.VNup_MB.*split.*phiN.*nh4.*n.MBA_N./(p.KsNH4_MB.*(denominator+n.MBA_N./p.KsNH4_MB));
    uptake3=p.VNup_MB.*(1-split).*phiN.*no3.*n.MBA_N./(p.KsNO3_MB.*(denominator+n.MBA_N./p.KsNO3_MB));
    raw4=uptake4;raw3=uptake3;
    if coupled
        uptake4=min(uptake4,inventory.NH4);uptake3=min(uptake3,inventory.NO3);
    end
    first4=uptake4;first3=uptake3;
    nn=n;
    nn.POMO_N=(n.POMO_N-nf.POMdec1)+addN.POMinp1+nf.MBA2POM1;
    nn.POMH_N=(n.POMH_N-nf.POMdec2)+addN.POMinp2+nf.MBA2POM2;
    nn.MOM_N=(n.MOM_N-nf.MOMdec)+nf.POMdec2MOM1+nf.POMdec2MOM2;
    nn.QOM_N=n.QOM_N+nf.DOM2QOMnet;
    nn.DOMT_N=domRemainingN+addN.DOMinp+nf.POMdec2DOM1+nf.POMdec2DOM2+nf.MOMdec ...
        +nf.MBA2DOM+nf.EP2DOM1+nf.EP2DOM2+nf.EM2DOM;
    nn.MBA_N=n.MBA_N+nf.DOM2MBA-nf.MBA_PM-nf.MBA2MBD+nf.MBD2MBA+uptake4+uptake3-releaseMB;
    nn.MBD_N=n.MBD_N+nf.MBA2MBD-nf.MBD2MBA;
    nn.ENZP1_N=n.ENZP1_N+nf.MBA2EP1-nf.EP2DOM1;
    nn.ENZP2_N=n.ENZP2_N+nf.MBA2EP2-nf.EP2DOM2;
    nn.ENZM_N=n.ENZM_N+nf.MBA2EM-nf.EM2DOM;
    releaseBD=max(nn.MBD_N-out.MBD_C./2,0);nn.MBD_N=nn.MBD_N-releaseBD;
    overflowBD=max(out.MBD_C-nn.MBD_N.*14,0);out.MBD_C=out.MBD_C-overflowBD;
    surplus=max(nn.MBA_N-out.MBA_C./2,0);nn.MBA_N=nn.MBA_N-surplus;releaseMB=releaseMB+surplus;
    available4=inventory.NH4+addN.NH4_dep-uptake4+releaseMB+releaseBD;
    available3=inventory.NO3+addN.NO3_dep-uptake3;
    assert(all(available4>=0)&all(available3>=0), ...
        'soil:microbe:sourceHourDIN','Source first uptake overdraws DIN.');
    deficit=max(out.MBA_C./14-nn.MBA_N,0);available=available4+available3;
    second=min(deficit,available);s4=zeros(size(second));s3=s4;take=second>0;
    s4(take)=second(take).*available4(take)./available(take);
    s3(take)=second(take).*available3(take)./available(take);
    uptake4=uptake4+s4;uptake3=uptake3+s3;nn.MBA_N=nn.MBA_N+second;
    overflowMB=max(out.MBA_C-nn.MBA_N.*14,0);out.MBA_C=out.MBA_C-overflowMB;
    out.CO2=out.CO2+overflowMB+overflowBD;
    somN=n.POMO_N+n.POMH_N+n.MOM_N+n.QOM_N+n.DOMT_N;
    assert(all(somN>0),'soil:microbe:sourceHourSOM','Source N-enzyme allocation requires positive SOM N.');
    factor=p.fpENZN.*(nh4+no3+din.NO2)./somN;
    enzProduction=factor.*(f.MBA2EP1+f.MBA2EP2+f.MBA2EM);
    enzTurnover=p.rEM.*c.ENZNmt_C;
    out.ENZNmt_C=c.ENZNmt_C+enzProduction-enzTurnover;
    out.MBA_C=out.MBA_C+enzTurnover-enzProduction;
    nn.ENZNmt_N=out.ENZNmt_C./3;
    nn.MBA_N=nn.MBA_N+(enzTurnover-enzProduction)./3;
    newdin=struct('NH4',inventory.NH4+addN.NH4_dep+releaseMB+releaseBD-uptake4, ...
        'NO3',inventory.NO3+addN.NO3_dep-uptake3);
    rc=out.CO2-c.CO2-addC.POMinp1-addC.POMinp2-addC.DOMinp;
    rn=newdin.NH4-inventory.NH4+newdin.NO3-inventory.NO3-addN.POMinp1-addN.POMinp2 ...
        -addN.DOMinp-addN.NH4_dep-addN.NO3_dep;
    for j=1:numel(pools)
        ck=[pools{j} '_C'];nk=[pools{j} '_N'];
        rc=rc+out.(ck)-c.(ck);rn=rn+nn.(nk)-n.(nk);
        validateResult(out.(ck),ck);validateResult(nn.(nk),nk);
    end
    validateResult(out.CO2,'CO2');validateResult(newdin.NH4,'NH4');validateResult(newdin.NO3,'NO3');
    assert(all(abs(rc)<1e-10)&all(abs(rn)<1e-10), ...
        'soil:microbe:sourceHourClosure','Source-hour C/N balance failed.');
    ledger=struct('Nim_NH4',uptake4,'Nim_NO3',uptake3,'Nmn_MBA',releaseMB, ...
        'Nmn_MBD',releaseBD,'CO2_ovflow_MBA',overflowMB,'CO2_ovflow_MBD',overflowBD, ...
        'raw_first_NH4',raw4,'raw_first_NO3',raw3, ...
        'first_NH4',first4,'first_NO3',first3,'second_NH4',s4,'second_NO3',s3, ...
        'second_uptake',second,'second_deficit',deficit,'phi_n',phiN, ...
        'enzyme_production_C',enzProduction,'enzyme_turnover_C',enzTurnover, ...
        'C_residual',rc,'N_residual',rn,'carbon_flux',f,'nitrogen_flux',nf);
end

function value=mm(v,s,e,k)
% Source per-hour MM cap; individual calls retain independent donor caps.
% Editor: Jizixing Zhu; Date: 2026-09-14.
    value=min(v.*s.*e./(k+s),s);
end

function value=boundedDonorN(carbonFlux,carbonStock,nitrogenStock)
% Paired N [mg N/g soil] for a nonnegative, donor-bounded C transfer. Columns
% may be empty subsets. Fraction-first arithmetic cannot exceed donor N;
% explicit depletion also leaves zero N when its paired C stock is exhausted.
% This is not a limiter for signed microbial or gross sorption transfers.
% Editor: Wenjie Liu; Date: 2026-10-09.
    assert(all(carbonFlux>=0&carbonFlux<=carbonStock), ...
        'soil:microbe:sourceHourDonorBudget','Unbounded substrate carbon transfer.');
    value=zeros(size(carbonFlux));active=carbonStock>0;
    value(active)=nitrogenStock(active).*(carbonFlux(active)./carbonStock(active));
    depleted=carbonFlux==carbonStock;value(depleted)=nitrogenStock(depleted);
end

function validateField(s,name,count,allowScalar)
% Validate explicit source columns; parameters may also be scalar.
% Editor: Jizixing Zhu; Date: 2026-09-14.
    if nargin<4,allowScalar=false;end
    assert(isstruct(s)&&isfield(s,name),'soil:microbe:sourceHourField', ...
        'Missing source-hour field %s.',name);
    v=s.(name);
    assert(isnumeric(v)&&isreal(v)&&iscolumn(v)&&all(isfinite(v))&&all(v>=0) ...
        &&(numel(v)==count||(allowScalar&&isscalar(v))), ...
        'soil:microbe:sourceHourValue','Invalid source-hour field %s.',name);
end

function validateResult(v,name)
% Never manufacture mass by clamping a result. DON must pass the same strict
% nonnegativity contract as the next source-hour input.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-10-09 Wenjie Liu -- reject every negative DON proposal.
    lower=-1e-12;if strcmp(name,'DOMT_N'),lower=0;end
    assert(all(isfinite(v))&&all(v>=lower), ...
        'soil:microbe:sourceHourNegative','Invalid source-hour result %s.',name);
end
