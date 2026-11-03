function [summary,detail]=maxRelativeFastPoolChange(built,startY,endY,absoluteFloor)
%MAXRELATIVEFASTPOOLCHANGE 比较同一强迫相位的微生物、酶与 DOM 快库漂移。
% built 仅提供网格、容重和 DOM 相参数；两个接受态提供各自的液态水容量。
% 门禁逐层包括 MB/BD、胞外及胞内酶的碳氮与三相，以及 DOC/DON 总量、
% 游离量、吸附量。慢库、DIN、累计通量不进入快库收敛门禁。
% AOM 逐组分漂移独立报告，不把快库达标误称为 AOM 或全系统周期平衡。
% 数值保留 owner 原生单位；绝对地板只识别零库存，不放宽相对收敛容差。
% Editor: Wenjie Liu; Date: 2026-10-23.
% Change: 2026-10-23 Wenjie Liu -- 显式快库门禁，液相/吸附相与 AOM 补充审计。
    if nargin<4||isempty(absoluteFloor),absoluteFloor=1e-18;end
    assert(isscalar(absoluteFloor)&&isfinite(absoluteFloor)&&absoluteFloor>0, ...
        'maxRelativeFastPoolChange:floor','库存零值地板必须为正有限标量。');
    n=built.domain.NL;rho=built.parameters.soil.bulk_density(:);
    statePath=strings(0,1);elementIndex=zeros(0,1);startC=zeros(0,1);endC=zeros(0,1);
    function appendStore(path,a,b)
        assert(isnumeric(a)&&isnumeric(b)&&isreal(a)&&isreal(b)&&isequal(size(a),size(b)) ...
            &&all(isfinite(a(:)))&&all(isfinite(b(:)))&&all(a(:)>=0)&&all(b(:)>=0), ...
            'maxRelativeFastPoolChange:state','库存 %s 的布局或数值无效。',path);
        a=a(:);b=b(:);count=numel(a);
        statePath=[statePath;repmat(string(path),count,1)];
        elementIndex=[elementIndex;(1:count)'];startC=[startC;a];endC=[endC;b];
    end
    for group={'MB','BD','ENZP1','ENZP2','EM','ENZN'}
        k=group{1};appendStore("soil_microbes."+k,startY.soil_microbes.(k),endY.soil_microbes.(k));
        nk=[k '_N'];appendStore("soil_organic."+nk,startY.soil_organic.(nk),endY.soil_organic.(nk));
        if ismember(k,{'MB','BD','ENZN'})
            for phase={'mobile','attached','retained'}
                ck=[k '_' phase{1}];nk=[ck '_N'];
                appendStore("soil_microbes."+ck,startY.soil_microbes.(ck),endY.soil_microbes.(ck));
                appendStore("soil_organic."+nk,startY.soil_organic.(nk),endY.soil_organic.(nk));
            end
        end
    end
    aTheta=phaseWater(startY,n);bTheta=phaseWater(endY,n);
    for kind={'doc','don'}
        k=kind{1};p=built.parameters.soil_dom;
        a=soil.dom.phaseState(startY.soil_dom,k,aTheta,rho,p);
        b=soil.dom.phaseState(endY.soil_dom,k,bTheta,rho,p);
        for field={'total_g_cm3_bulk','free_g_cm3_bulk','sorbed_g_cm3_bulk'}
            f=field{1};appendStore("soil_dom."+k+"."+f,a.(f),b.(f));
        end
    end
    gateCount=numel(startC);
    for material={'litter','straw','manure'}
        k=material{1};a=startY.soil_organic.aom.(k);b=endY.soil_organic.aom.(k);
        assert(isfield(a,'components')==isfield(b,'components'), ...
            'maxRelativeFastPoolChange:AOM','AOM 组分布局发生变化。');
        if isfield(a,'components')
            a=soil.litter.componentTotals(a,true);b=soil.litter.componentTotals(b,true);
            for component={'metabolic','cellulose','lignin','cwd'}
                c=component{1};addAom("soil_organic.aom."+k+"."+c,a.components.(c),b.components.(c));
            end
        else
            addAom("soil_organic.aom."+k,a,b);
        end
    end
    function addAom(path,a,b)
        for field={'surface_C','surface_N','element_C','element_N'}
            f=field{1};appendStore(path+"."+f,a.(f),b.(f));
        end
    end
    change=abs(endC-startC);relative=zeros(size(change));regular=abs(startC)>absoluteFloor;
    relative(regular)=change(regular)./abs(startC(regular));
    relative(~regular&change>absoluteFloor)=Inf;
    inGate=(1:numel(startC))'<=gateCount;
    detail=table(statePath,elementIndex,startC,endC,relative,inGate, ...
        'VariableNames',{'state_path','element_index','C_start','C_end','relative_change','in_convergence_gate'});
    [maximum,index]=max(relative(1:gateCount));
    [aomMaximum,aomIndex]=max(relative(gateCount+1:end));aomIndex=aomIndex+gateCount;
    summary=struct('max_relative_change',maximum,'state_path',char(statePath(index)), ...
        'element_index',elementIndex(index),'C_start',startC(index),'C_end',endC(index), ...
        'absolute_floor_native_units',absoluteFloor,'n_compared_states',gateCount, ...
        'scope','microbial_and_enzyme_CN_phases_and_DOM_total_free_sorbed', ...
        'slow_pools_and_DIN_in_gate',false,'aom_in_gate',false, ...
        'aom_max_relative_change',aomMaximum,'aom_state_path',char(statePath(aomIndex)), ...
        'aom_element_index',elementIndex(aomIndex), ...
        'interpretation','快库门禁不证明 AOM、慢库、无机氮或完整生态系统达到周期平衡。');
end

function theta=phaseWater(Y,n)
    assert(isequal(size(Y.Theta_L),[n,2])&&all(isfinite(Y.Theta_L(:)))&&all(Y.Theta_L(:)>=0), ...
        'maxRelativeFastPoolChange:water','DOM 相核对要求各自接受态液态水。');
    theta=mean(Y.Theta_L,2);
end
