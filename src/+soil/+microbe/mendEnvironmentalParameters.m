function p = mendEnvironmentalParameters(base, temperatureC, swpMPa, pH)
%MENDENVIRONMENTALPARAMETERS Evaluate source MEND environmental responses.
%
% Inputs are top-down column vectors. Returned kinetic fields use the MEND
% hourly clock. This is the Arrhenius (`iTmp_Func=0`) and MM
% (`iKinetics=0`) path from MOD_MEND.F90 commit 92323c7.
%
% Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-21 Wenjie Liu -- make the configured dormancy switch gate
%   dormant maintenance and active/dormant exchange in both MEND solvers.
% Change: 2026-09-29 Wenjie Liu -- accept top-down NL-by-1 kinetic profiles
%   (NL-by-2 for POM pairs), with scalar legacy expansion at this adapter.
% Change: 2026-10-12 Wenjie Liu -- validate and expand layer-capable CUE and
% dormancy thresholds with the same explicit top-down profile contract.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- bypass only the three pH scalars when disabled.
% Change: 2026-10-27 Wenjie Liu -- map dormancy and DIN kinetics per layer.
    temperatureC=temperatureC(:); swpMPa=swpMPa(:); pH=pH(:);
    assert(isequal(size(temperatureC),size(swpMPa),size(pH)) ...
        &&all(isfinite([temperatureC;swpMPa;pH]))&&all(pH>0), ...
        'soil:microbe:mendEnvironment','Invalid MEND environmental drivers.');
    n=numel(temperatureC);
    dormancyEnabled=true;
    if isfield(base,'dormancy_enabled')
        flag=base.dormancy_enabled;
        assert(isscalar(flag)&&(islogical(flag) ...
            ||(isnumeric(flag)&&isreal(flag)&&isfinite(flag) ...
            &&ismember(flag,[0 1]))), ...
            'soil:microbe:mendDormancySwitch', ...
            'dormancy_enabled must be a scalar logical or zero/one flag.');
        dormancyEnabled=logical(flag);
    end
    p=base;
    % Change: 2026-10-12 Wenjie Liu -- identification-only shape controls for
    % the source fpH0 response. Defaults (scale=1, adapt=0) reproduce the
    % frozen source optimum/sensitivity pair bitwise.
    phScale=1.0; phAdapt=0.0;
    if isfield(base,'ph_sensitivity_scale'), phScale=base.ph_sensitivity_scale; end
    if isfield(base,'ph_optimum_adapt_per_ph'), phAdapt=base.ph_optimum_adapt_per_ph; end
    assert(isscalar(phScale)&&isfinite(phScale)&&phScale>=0.5&&phScale<=4 ...
        &&isscalar(phAdapt)&&isfinite(phAdapt)&&phAdapt>=0&&phAdapt<=0.5, ...
        'soil:microbe:mendPhShape','Invalid pH-response identification controls.');
    optimum=@(center)center+phAdapt.*(pH-5);
    phEnabled=true;
    if isfield(base,'ph_response_enabled')
        flag=base.ph_response_enabled;
        assert(isscalar(flag)&&(islogical(flag) ...
            ||(isnumeric(flag)&&isreal(flag)&&isfinite(flag)&&ismember(flag,[0 1]))), ...
            'soil:microbe:mendPhSwitch','Invalid scalar pH-response switch.');
        phEnabled=logical(flag);
    end
    phPOM1=ones(size(pH));phPOM2=phPOM1;phMOM=phPOM1;
    if phEnabled
        phPOM1=phScalar(pH,optimum(4.2),1.4.*phScale);
        phPOM2=phScalar(pH,optimum(5.3),1.7.*phScale);
        phMOM=phScalar(pH,optimum(4.8),1.6.*phScale);
    end
    vd=pairProfile(base.VdPOM_ref_h,n);ks=pairProfile(base.KsPOM_ref_mg_g,n);
    p.VdPOM_h=[vd(:,1).*arrhenius(temperatureC,53.0) ...
        .*swpOptimal(swpMPa).*phPOM1, ...
        vd(:,2).*arrhenius(temperatureC,36.3) ...
        .*swpIncreasing(swpMPa,base.biome,base.som) ...
        .*phPOM2];
    kmT=arrhenius(temperatureC,30.0);
    p.KsPOM_mg_g=kmT.*ks;
    p.VdMOM_h=layerProfile(base.VdMOM_ref_h,n).*arrhenius(temperatureC,47.0) ...
        .*swpIncreasing(swpMPa,base.biome,base.som) ...
        .*phMOM;
    p.KsMOM_mg_g=layerProfile(base.KsMOM_ref_mg_g,n).*kmT;
    vm=layerProfile(base.Vm_ref_h,n);
    p.Vg_h=layerProfile(base.Vg_ref_h,n).*arrhenius(temperatureC,47.0);
    p.Vm_h=vm.*arrhenius(temperatureC,20.0);
    p.KsDOM_mg_g=layerProfile(base.KsDOM_ref_mg_g,n).*kmT;
    p.Yg=min(0.9,max(0.01,layerProfile(base.Yg_ref,n) ...
        -base.Yg_slope_per_C.*(temperatureC-20)));
    swpA2DMPa=layerProfile(base.SWP_A2D_MPa,n);
    tau=layerProfile(base.tau,n);
    wdorm=layerProfile(base.wdorm,n);
    p.wdorm=wdorm;
    % Preserve the original scalar representation for legacy scalar callers.
    if isscalar(base.wdorm)&&(~isfield(base,'reaction_method') ...
            ||strcmp(string(base.reaction_method),"rk4_equilibrium"))
        p.wdorm=base.wdorm;
    end
    dry=swpA2D(swpMPa,swpA2DMPa,wdorm);
    wet=swpD2A(swpMPa,tau.*swpA2DMPa,wdorm);
    p.rMORT_h=layerProfile(base.rMORT_ref_h,n).*dry;
    p.VmD_h=layerProfile(base.VmD_ref_h,n);
    p.VmA2D_h=vm.*arrhenius(temperatureC,20.0).*dry;
    p.VmD2A_h=vm.*arrhenius(temperatureC,20.0).*wet;
    if ~dormancyEnabled
        p.VmD_h(:)=0;
        p.VmA2D_h(:)=0;
        p.VmD2A_h(:)=0;
    end
    p.VNup_MB_h=layerProfile(base.VNup_MB_ref_h,n).*arrhenius(temperatureC,47.0);
    p.KsNH4_MB_mg_g=layerProfile(base.KsNH4_MB_ref_mg_g,n).*kmT;
    p.KsNO3_MB_mg_g=layerProfile(base.KsNO3_MB_ref_mg_g,n).*kmT;
    p.temperature_C=temperatureC; p.swp_MPa=swpMPa; p.pH=pH;
    p.ph_scalar_pom=[phPOM1,phPOM2];p.ph_scalar_mom=phMOM;
end

function value=layerProfile(value,n)
% Only scalar legacy values or explicit top-down columns are accepted.
    assert(isnumeric(value)&&isreal(value)&&all(isfinite(value),'all') ...
        &&(isscalar(value)||isequal(size(value),[n,1])), ...
        'soil:microbe:mendKineticShape','Expected scalar or top-down NL-by-1 kinetic parameter.');
    if isscalar(value),value=repmat(value,n,1);end
end

function value=pairProfile(value,n)
% Columns are POM1/POM2, rows are top-down layers; never linear indexing.
    assert(isnumeric(value)&&isreal(value)&&all(isfinite(value),'all') ...
        &&(isequal(size(value),[1,2])||isequal(size(value),[n,2])), ...
        'soil:microbe:mendKineticShape','Expected 1-by-2 or top-down NL-by-2 POM parameter.');
    if size(value,1)==1,value=repmat(value,n,1);end
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
            if som=="SOD"
                minimum=-1.71; w=1.43;
            elseif som=="SOL"
                minimum=-22.61; w=1.11;
            elseif som=="LIT"
                minimum=-39.73; w=0.89;
            end
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
