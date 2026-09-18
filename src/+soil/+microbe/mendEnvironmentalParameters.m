function p = mendEnvironmentalParameters(base, temperatureC, swpMPa, pH)
%MENDENVIRONMENTALPARAMETERS Evaluate source MEND environmental responses.
%
% Inputs are top-down column vectors. Returned kinetic fields use the MEND
% hourly clock. This is the Arrhenius (`iTmp_Func=0`) and MM
% (`iKinetics=0`) path from MOD_MEND.F90 commit 92323c7.
%
% Editor: Jizixing Zhu; Date: 2026-09-12.
    temperatureC=temperatureC(:); swpMPa=swpMPa(:); pH=pH(:);
    assert(isequal(size(temperatureC),size(swpMPa),size(pH)) ...
        &&all(isfinite([temperatureC;swpMPa;pH]))&&all(pH>0), ...
        'soil:microbe:mendEnvironment','Invalid MEND environmental drivers.');
    n=numel(temperatureC);
    p=base;
    p.VdPOM_h=[base.VdPOM_ref_h(1).*arrhenius(temperatureC,53.0) ...
        .*swpOptimal(swpMPa).*phScalar(pH,4.2,1.4), ...
        base.VdPOM_ref_h(2).*arrhenius(temperatureC,36.3) ...
        .*swpIncreasing(swpMPa,base.biome,base.som) ...
        .*phScalar(pH,5.3,1.7)];
    kmT=arrhenius(temperatureC,30.0);
    p.KsPOM_mg_g=kmT.*repmat(base.KsPOM_ref_mg_g,n,1);
    p.VdMOM_h=base.VdMOM_ref_h.*arrhenius(temperatureC,47.0) ...
        .*swpIncreasing(swpMPa,base.biome,base.som) ...
        .*phScalar(pH,4.8,1.6);
    p.KsMOM_mg_g=base.KsMOM_ref_mg_g.*kmT;
    p.Vg_h=base.Vg_ref_h.*arrhenius(temperatureC,47.0);
    p.Vm_h=base.Vm_ref_h.*arrhenius(temperatureC,20.0);
    p.KsDOM_mg_g=base.KsDOM_ref_mg_g.*kmT;
    p.Yg=min(0.9,max(0.01,base.Yg_ref ...
        -base.Yg_slope_per_C.*(temperatureC-20)));
    dry=swpA2D(swpMPa,base.SWP_A2D_MPa,base.wdorm);
    wet=swpD2A(swpMPa,base.tau.*base.SWP_A2D_MPa,base.wdorm);
    p.rMORT_h=base.rMORT_ref_h.*dry;
    p.VmD_h=repmat(base.VmD_ref_h,n,1);
    p.VmA2D_h=base.Vm_ref_h.*arrhenius(temperatureC,20.0).*dry;
    p.VmD2A_h=base.Vm_ref_h.*arrhenius(temperatureC,20.0).*wet;
    p.VNup_MB_h=base.VNup_MB_ref_h.*arrhenius(temperatureC,47.0);
    p.KsNH4_MB_mg_g=base.KsNH4_MB_ref_mg_g.*kmT;
    p.KsNO3_MB_mg_g=base.KsNO3_MB_ref_mg_g.*kmT;
    p.temperature_C=temperatureC; p.swp_MPa=swpMPa; p.pH=pH;
end

function value=arrhenius(t,eakJmol)
% Source fTArh0 at the fixed MEND 20 degC reference.
    r=8.31446261815324;
    value=exp(eakJmol.*1000./r.*(1./293.15-1./(t+273.15)));
end

function value=phScalar(pH,optimum,sensitivity)
% Source fpH0.
    value=exp(-((pH-optimum)./sensitivity).^2);
end

function value=swpIncreasing(swp,biome,som)
% Source fSWP/fSWP0 for the selected biome and organic-matter class.
    [minimum,w]=swpShape(upper(string(biome)),upper(string(som)));
    fieldCapacity=-0.033;
    value=ones(size(swp));
    value(swp<minimum)=0;
    middle=swp>=minimum&swp<fieldCapacity;
    value(middle)=1-(log(swp(middle)./fieldCapacity) ...
        ./log(minimum/fieldCapacity)).^w;
    value=min(1,max(0,value));
end

function [minimum,w]=swpShape(biome,som)
    switch som
        case "SOD", minimum=-1.71; w=1.43;
        case "SOL", minimum=-13.86; w=1.20;
        case "LIT", minimum=-36.49; w=1.04;
        otherwise
            error('soil:microbe:mendSom','Unknown MEND SOM class %s.',som);
    end
    switch biome
        case "ASM"
            if som=="SOL", minimum=-10.95; w=1.26; end
        case "MGC"
            if som=="SOD", minimum=-1.71; w=1.43;
            elseif som=="SOL", minimum=-22.61; w=1.11;
            elseif som=="LIT", minimum=-39.73; w=0.89; end
        case "MDF"
            if som=="SOL", minimum=-4.97; w=1.07;
            elseif som=="LIT", minimum=-29.00; w=1.27; end
        case "MCF"
            if som=="SOL", minimum=-8.24; w=1.40;
            elseif som=="LIT", minimum=-39.85; w=1.06; end
        otherwise
            error('soil:microbe:mendBiome','Unknown MEND biome class %s.',biome);
    end
end

function value=swpOptimal(swp)
% Source fSWP_OPT active implementation.
    minimum=-10^2.5; low=-10^-1.5; high=-10^-2.5; maximum=-10^-4;
    value=zeros(size(swp));
    i=swp>minimum&swp<=low;
    value(i)=0.625-0.25.*log10(abs(swp(i)));
    i=swp>low&swp<=high; value(i)=1;
    i=swp>high&swp<=maximum;
    value(i)=(2.5+0.4.*log10(abs(swp(i))))./1.5;
    value(swp>maximum)=0.6;
    value=min(1,max(0,value));
end

function value=swpA2D(swp,threshold,w)
% Source fSWP_A2D with an explicit zero/zero boundary.
    a=abs(swp).^w; b=abs(threshold).^w;
    value=a./(a+b);
    value(a+b==0)=0;
end

function value=swpD2A(swp,threshold,w)
% Source fSWP_D2A with an explicit zero/zero boundary.
    a=abs(swp).^w; b=abs(threshold).^w;
    value=b./(a+b);
    value(a+b==0)=1;
end
