function report=biogeochemistryBudget(built,startY,endY,options)
%BIOGEOCHEMISTRYBUDGET Independent coupled accepted-interval C/N balance.
% Inputs: immutable built, two accepted full Y endpoints, and options with
% carbon_tolerance_g_cm2/nitrogen_tolerance_g_cm2 (explicit absolute mass
% tolerances); start_is_fresh=true permits missing first-step N cumulative
% fields ONLY at the zero-flux initial endpoint. Outputs use g element/cm2.
% 物理存量变化等于真实外部输入减输出。全耦合时根吸氮和凋落为内部转移；
% 土壤控制体中根吸氮为外部输出。矿化与尿素水解始终为内部转移。
% C outputs include accepted plant respiration and MEND/AOM respiration;
% no prognostic CO2-C pool is owned by the current soil-gas model. N gaseous
% losses are NH3, nitrification N2O and total denitrification N. All inorganic,
% dissolved-organic and mobile-particle face routes are accounted explicitly.
% Negative xsmr and its writeoff are reported separately: a closed accounting
% balance cannot certify physical C closure. No closure residual supplies a
% missing input/output. Unknown/missing coverage is an error, never a zero.
% 适用范围：当前逐小时 MEND 土柱，BNF 关闭，DOM/颗粒边界拥有累计通量账本。
% PlantCN 与凋落耦合同时开启时包含植株；否则控制体仅包含土柱与地表 AOM。
% This diagnostic neither mutates state nor changes an acceptance tolerance.
% Editor: Wenjie Liu; Date: 2026-09-30.
% Change: 2026-10-17 Wenjie Liu -- read disjoint manure/straw inputs;
%   audit total-minus-litter closure without a cancellation-based route.
% Change: 2026-10-23 Wenjie Liu -- 新增土壤控制体的独立收支审计，
%   将实际根吸氮列为外部输出；保留全部真实边界、AOM 呼吸与原守恒容差。
    required={'carbon_tolerance_g_cm2','nitrogen_tolerance_g_cm2'};
    for j=1:numel(required)
        k=required{j};assert(isfield(options,k)&&isscalar(options.(k)) ...
            &&isfinite(options.(k))&&options.(k)>0, ...
            'biogeochemistryBudget:tolerance','Explicit positive %s is required.',k);
    end
    fresh=isfield(options,'start_is_fresh')&&isequal(options.start_is_fresh,true);
    sw=built.control.switches;
    assert(isfield(sw,'lPlantCN')&&isfield(sw,'lLitterInput'), ...
        'biogeochemistryBudget:coverage','PlantCN and litter switches are required.');
    includePlant=sw.lPlantCN==1&&sw.lLitterInput==1;
    assert(~isfield(sw,'lBNF')||sw.lBNF==0, ...
        'biogeochemistryBudget:coverage','BNF requires a separately verified external ledger.');
    first=diagnostics.biogeochemistryInventory(built,startY);
    last=diagnostics.biogeochemistryInventory(built,endY);
    if includePlant,a0=startY.plant_cn.accum;a=endY.plant_cn.accum;end
    o0=startY.soil_organic.cumulative;o=endY.soil_organic.cumulative;
    f0=startY.soil_nitrogen.flux;f=endY.soil_nitrogen.flux;
    if fresh
        if includePlant
            names=fieldnames(a0);
            for j=1:numel(names),assert(a0.(names{j})==0, ...
                    'biogeochemistryBudget:notFresh','Fresh PlantCN cumulative state is nonzero.');end
        end
        names=fieldnames(o0);
        for j=1:numel(names),assert(o0.(names{j})==0, ...
                'biogeochemistryBudget:notFresh','Fresh organic cumulative state is nonzero.');end
        keys={'cvTop','cvBot','cvTop1','cvBot1','cvSnit','cvSmin','cvSden','cvSvol', ...
            'cvSinkSnh4','cvSinkSno3','Cvnit_N2o','Cvdenit_N2o','Cvdenit_N2'};
        for j=1:numel(keys),assert(isfield(f0,keys{j})&&all(f0.(keys{j})(:)==0), ...
                'biogeochemistryBudget:notFresh','Fresh N flux memory is missing/nonzero.');end
    end
    ci=struct();co=struct();ni=struct();no=struct();
    if includePlant
        ci.gpp=delta(a0,a,'cum_GPP')/1e4;
        co.plant_respiration=delta(a0,a,'cum_R_a')/1e4;
        co.harvest=delta(a0,a,'cum_C_harvest')/1e4;
        no.harvest=delta(a0,a,'cum_N_harvest')/1e4;
    else
        % 只读取实际施加的根吸氮累计账本，不以植株需求或残差推算。
        no.root_uptake=(ndelta('cum_cvSinkSnh4')+ndelta('cum_cvSinkSno3'))/1e6;
    end
    % Change: 2026-10-30 Wenjie Liu -- paired prescribed root residue is an
    % external soil input, separate from the internal PlantCN litter ledger.
    rootC=0;rootN=0;
    rootKeys={'input_root_residue_C_g_cm2','input_root_residue_N_g_cm2'};
    hasRootLedger=any(isfield(o0,rootKeys))||any(isfield(o,rootKeys));
    rootScheduled=isfield(built,'sources') ...
        &&isfield(built.sources,'organic_fertilization_schedule') ...
        &&any(strcmp(built.sources.organic_fertilization_schedule.material_type,'root_residue'));
    if hasRootLedger||rootScheduled
        rootC=delta(o0,o,rootKeys{1});rootN=delta(o0,o,rootKeys{2});
    end
    ci.organic_management=delta(o0,o,'input_manure_C_g_cm2')+delta(o0,o,'input_straw_C_g_cm2')+rootC;
    ni.organic_management=delta(o0,o,'input_manure_N_g_cm2')+delta(o0,o,'input_straw_N_g_cm2')+rootN;
    co.aom_respiration=delta(o0,o,'co2_C_g_cm2');
    co.mend_respiration=last.mend_respired_carbon_g_cm2-first.mend_respired_carbon_g_cm2;
    ni.mineral_fertilizer=0;
    for key={'cum_Fert_NH4','cum_Fert_NO3','cum_Fert_Urea'}
        ni.mineral_fertilizer=ni.mineral_fertilizer+ndelta(key{1})/1e6;
    end
    no.ammonia_volatilization=ndelta('cum_cvSvol')/1e6;
    no.nitrification_n2o=ndelta('cum_Cvnit_N2o')/1e6;
    no.denitrification_total=ndelta('cum_cvSden')/1e6;
    for carrier={'NH4','NO3'}
        key=carrier{1};small=lower(key);
        ni.([small '_top'])=ndelta(['cum_top_input_' key])/1e6;
        ni.([small '_bottom'])=ndelta(['cum_bottom_input_' key])/1e6;
        no.([small '_top'])=ndelta(['cum_top_loss_' key])/1e6;
        no.([small '_bottom'])=ndelta(['cum_bottom_loss_' key])/1e6;
    end
    for carrier={'doc','don'}
        key=carrier{1};[input,output]=organicBoundary(startY.soil_dom,endY.soil_dom, ...
            key,[key '_g_cm2'],built.parameters.soil_dom.boundary);
        if strcmp(key,'doc')
            ci.doc_boundary=input;co.doc_boundary=output;
        else
            ni.don_boundary=input;no.don_boundary=output;
        end
    end
    if built.parameters.soil_microbes.transport_enabled
        p=built.parameters.soil_microbes.bioparticle.boundary;
        for group={'MB','BD','ENZN'}
            g=group{1};
            for element={'carbon','nitrogen'}
                e=element{1};[input,output]=particleBoundary(startY.soil_microbes, ...
                    endY.soil_microbes,g,e,p);
                key=['particle_' lower(g)];
                if strcmp(e,'carbon'),ci.(key)=input;co.(key)=output;
                else,ni.(key)=input;no.(key)=output;end
            end
        end
    end
    % 内部路径核对独立于总收支，不能将检查残差用于补齐外源通量。
    plantLitterC=0;plantLitterN=0;
    if includePlant
        for route={'leaf_litter','root_litter','prune_litter','mortality'}
            plantLitterC=plantLitterC+delta(a0,a,['cum_C_' route{1}])/1e4;
            plantLitterN=plantLitterN+delta(a0,a,['cum_N_' route{1}])/1e4;
        end
    end
    transfer=struct('litter_carbon_g_cm2',delta(o0,o,'input_litter_C_g_cm2')-plantLitterC, ...
        'litter_nitrogen_g_cm2',delta(o0,o,'input_litter_N_g_cm2')-plantLitterN, ...
        'organic_source_carbon_g_cm2',delta(o0,o,'input_C_g_cm2') ...
            -delta(o0,o,'input_litter_C_g_cm2')-ci.organic_management, ...
        'organic_source_nitrogen_g_cm2',delta(o0,o,'input_N_g_cm2') ...
            -delta(o0,o,'input_litter_N_g_cm2')-ni.organic_management, ...
        'nitrification_product_g_cm2',no.nitrification_n2o-0.02*ndelta('cum_cvSnit')/1e6);
    if includePlant
        transfer.root_nitrogen_g_cm2=delta(a0,a,'cum_U_N')/1e4 ...
            -(ndelta('cum_cvSinkSnh4')+ndelta('cum_cvSinkSno3'))/1e6;
    end
    for terms={ci,co,ni,no}
        assert(all(cell2mat(struct2cell(terms{1}))>=0), ...
            'biogeochemistryBudget:nonmonotone','An external route decreased or became negative.');
    end
    netC=sum(cell2mat(struct2cell(ci)))-sum(cell2mat(struct2cell(co)));
    netN=sum(cell2mat(struct2cell(ni)))-sum(cell2mat(struct2cell(no)));
    debtChange=last.plant_carbon_debt_g_cm2-first.plant_carbon_debt_g_cm2;
    writeoff=0;
    if includePlant,writeoff=delta(a0,a,'cum_C_xsmr_writeoff')/1e4;end
    residualC=last.physical_carbon_g_cm2-first.physical_carbon_g_cm2-netC;
    residualN=last.nitrogen_g_cm2-first.nitrogen_g_cm2-netN;
    accountingC=last.accounting_carbon_g_cm2-first.accounting_carbon_g_cm2-netC-writeoff;
    tolC=options.carbon_tolerance_g_cm2;tolN=options.nitrogen_tolerance_g_cm2;
    transferPass=all(abs([transfer.litter_carbon_g_cm2, ...
        transfer.organic_source_carbon_g_cm2])<=tolC) ...
        &&all(abs([transfer.litter_nitrogen_g_cm2,transfer.organic_source_nitrogen_g_cm2, ...
            transfer.nitrification_product_g_cm2])<=tolN);
    if includePlant,transferPass=transferPass&&abs(transfer.root_nitrogen_g_cm2)<=tolN;end
    report=struct('schema_version',1,'units','g element/cm2 ground', ...
        'control_volume',first.control_volume, ...
        'coverage_complete',true,'initial',first,'final',last, ...
        'carbon_inputs_g_cm2',ci,'carbon_outputs_g_cm2',co, ...
        'nitrogen_inputs_g_cm2',ni,'nitrogen_outputs_g_cm2',no, ...
        'carbon_physical_residual_g_cm2',residualC,'nitrogen_residual_g_cm2',residualN, ...
        'carbon_accounting_residual_g_cm2',accountingC, ...
        'carbon_debt_change_g_cm2',debtChange,'carbon_debt_writeoff_g_cm2',writeoff, ...
        'carbon_tolerance_g_cm2',tolC,'nitrogen_tolerance_g_cm2',tolN, ...
        'internal_transfer_residuals',transfer,'internal_transfers_passed',transferPass, ...
        'carbon_physical_passed',abs(residualC)<=tolC, ...
        'carbon_accounting_passed',abs(accountingC)<=tolC, ...
        'nitrogen_passed',abs(residualN)<=tolN, ...
        'writeoff_is_physical_input',false,'component_residuals_used',false);
    report.plant_included=includePlant;
    report.root_uptake_is_external_output=~includePlant;
    report.passed=report.carbon_physical_passed&&report.carbon_accounting_passed ...
        &&report.nitrogen_passed&&transferPass;
    report.status='FAIL';if report.passed,report.status='PASS';end
    function value=ndelta(key)
        if ~isfield(f0,key)&&fresh,start=0;else,start=scalar(f0,key);end
        value=scalar(f,key)-start;
        assert(value>=0,'biogeochemistryBudget:nonmonotone', ...
            'Accepted N cumulative route %s decreased.',key);
    end
end

function [input,output]=organicBoundary(a,b,key,exportKey,p)
    assert(isfield(p,'schema_version')&&p.schema_version==2 ...
        &&~strcmp(p.bottom_kind,'closed'), ...
        'biogeochemistryBudget:coverage','Numbered open organic boundary ledger required.');
    before=a.boundary_exchange_cumulative.(key);after=b.boundary_exchange_cumulative.(key);
    input=delta(before,after,'top_inflow_g_cm2')+delta(before,after,'bottom_inflow_g_cm2');
    output=delta(before,after,'top_outflow_g_cm2') ...
        +delta(a.bottom_export_cumulative,b.bottom_export_cumulative,exportKey);
end
function [input,output]=particleBoundary(a,b,group,element,p)
    aa=a;bb=b;
    aa.bottom_export_cumulative=a.bottom_export_cumulative.(group);
    bb.bottom_export_cumulative=b.bottom_export_cumulative.(group);
    [input,output]=organicBoundary(aa,bb,[group '_' element],[element '_g_cm2'],p);
end
function value=delta(a,b,key)
    value=scalar(b,key)-scalar(a,key);
    assert(value>=0,'biogeochemistryBudget:nonmonotone','Accepted cumulative %s decreased.',key);
end
function value=scalar(s,key)
    assert(isfield(s,key)&&isnumeric(s.(key))&&isscalar(s.(key)) ...
        &&isreal(s.(key))&&isfinite(s.(key)), ...
        'biogeochemistryBudget:missingLedger','Missing/nonfinite accepted ledger %s.',key);
    value=double(s.(key));
end
