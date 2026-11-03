function [work,input_g_cm2,externalNH4_g_cm3,diagnostic]=stageManureRoutes( ...
        accepted,context,dz)
%STAGEMANUREROUTES Stage one passive manure-source cohort on applied routes.
% accepted is a manure_organic/manure_tan tracer. context.management carries
% actual accepted candidate event increments: bottom-up NL-by-1 g N cm^-3
% manure_fast_don_N_g_cm3, manure_aom_organic_N_added_g_cm3 and
% manure_tan_N_added_g_cm3, plus manure_surface_aom_N_added_g_cm2.
% context.aom_before/after are the total manure AOM before/after this step's
% decomposition, after management staging. context.manure_routes contains
% producer-owned gross element decomposition/leaching [g N cm^-3] and
% surface decomposition/leaching/PON/TAN [g N cm^-2].
% context.partition_pon_N_g_cm3 is the actual donor of context.aom_partition
% (all-material PON in legacy mixed mode, manure PON in material-specific
% mode). These applied N_PON1/2 routes are consumed without inferring gross
% reactions from endpoint differences. Newly tagged DON precedes transport.
% externalNH4 is an input to the later frozen ADE, never same-step root/MEND
% inventory. All changes are passive; no total pool or event flag is edited.
% Editor: Wenjie Liu; Date: 2026-10-30.
    dz=dz(:);n=numel(dz);work=accepted;
    required={'source_kind','management','aom_before','aom_after', ...
        'manure_routes','partition_pon_N_g_cm3','aom_partition'};
    assert(isstruct(context)&&all(isfield(context,required)) ...
        &&isfield(accepted,'source_kind')&&strcmp(accepted.source_kind,context.source_kind) ...
        &&ismember(context.source_kind,{'manure_organic','manure_tan'}), ...
        'soil:nitrogen:tracer:manureContext','Applied manure cohort context is incomplete.');
    assert(all(isfinite(dz)&dz>0)&&isfield(accepted,'aom'), ...
        'soil:nitrogen:tracer:manureOwner','Manure AOM tracer and positive layer thickness are required.');
    management=context.management;
    fast=vector(management,'manure_fast_don_N_g_cm3',n);
    added=vector(management,'manure_aom_organic_N_added_g_cm3',n);
    addedTAN=vector(management,'manure_tan_N_added_g_cm3',n);
    surface=scalar(management,'manure_surface_aom_N_added_g_cm2');
    input_g_cm2=0;externalNH4_g_cm3=zeros(n,1);
    tagPon=zeros(n,1);tagDon=zeros(n,1);tagTan=zeros(n,1);tagFast=zeros(n,1);
    if strcmp(context.source_kind,'manure_tan')
        externalNH4_g_cm3=addedTAN;
        input_g_cm2=sum(addedTAN.*dz);
        % TAN labels can recycle into MEND pools, but never into AOM.
        assert(accepted.aom.manure.surface_N==0 ...
            &&all(accepted.aom.manure.element_N==0), ...
            'soil:nitrogen:tracer:TANAom','The native total system has no TAN-to-AOM route.');
        residual=0;
    else
        before=context.aom_before;after=context.aom_after;r=context.manure_routes;
        totalElement=vector(before,'element_N',n);totalSurface=scalar(before,'surface_N');
        stagedElement=accepted.aom.manure.element_N(:)+added;
        stagedSurface=accepted.aom.manure.surface_N+surface;
        fraction=donorFraction(stagedElement,totalElement,'incorporated manure');
        surfaceFraction=donorFraction(stagedSurface,totalSurface,'surface manure');
        decomp=vector(r,'element_decomp_N_g_cm3',n);
        leach=vector(r,'element_leach_N_g_cm3',n);
        surfaceDecomp=scalar(r,'surface_decomp_N_g_cm2');
        surfaceLeach=scalar(r,'surface_leach_N_g_cm2');
        surfacePon=scalar(r,'surface_pon_N_g_cm2');
        surfaceTan=scalar(r,'surface_tan_N_g_cm2');
        assert(abs(surfacePon+surfaceTan-surfaceDecomp)<= ...
            512*eps(max([surfaceDecomp;realmin])), ...
            'soil:nitrogen:tracer:surfacePartition','Surface AOM gross N outputs do not close.');
        remaining=totalElement-decomp-leach;
        surfaceRemaining=totalSurface-surfaceDecomp-surfaceLeach;
        compareRemaining(remaining,vector(after,'element_N',n),'element manure AOM');
        compareRemaining(surfaceRemaining,scalar(after,'surface_N'),'surface manure AOM');
        work.aom.manure.element_N=remaining.*fraction;
        work.aom.manure.surface_N=surfaceRemaining.*surfaceFraction;
        tagPon=decomp.*fraction;
        tagDon=leach.*fraction;
        tagPon(end)=tagPon(end)+surfacePon.*surfaceFraction./dz(end);
        tagDon(end)=tagDon(end)+surfaceLeach.*surfaceFraction./dz(end);
        tagTan(end)=surfaceTan.*surfaceFraction./dz(end);
        partitionDonor=vector(context,'partition_pon_N_g_cm3',n);
        ponFraction=donorFraction(tagPon,partitionDonor,'applied AOM PON partition');
        one=vector(context.aom_partition,'N_PON1',n);
        two=vector(context.aom_partition,'N_PON2',n);
        assert(all(abs(one+two-partitionDonor)<=512*eps(max([partitionDonor;realmin]))), ...
            'soil:nitrogen:tracer:particulatePartition','Applied PON partition differs from its actual donor.');
        work.organic_n.PON1=work.organic_n.PON1+one.*ponFraction;
        work.organic_n.PON2=work.organic_n.PON2+two.*ponFraction;
        work.organic_n.PON=work.organic_n.PON1+work.organic_n.PON2;
        tagFast=fast;
        work.dom.don_total_g_cm3=work.dom.don_total_g_cm3+tagDon+tagFast;
        externalNH4_g_cm3=tagTan;
        input_g_cm2=sum((added+fast).*dz)+surface;
        priorAom=accepted.aom.manure.surface_N+sum(accepted.aom.manure.element_N(:).*dz);
        nextAom=work.aom.manure.surface_N+sum(work.aom.manure.element_N.*dz);
        residual=priorAom+sum(added.*dz)+surface-nextAom-sum((tagPon+tagDon+tagTan).*dz);
        assert(abs(residual)<=2e-12*max(1,input_g_cm2), ...
            'soil:nitrogen:tracer:manureAomClosure','Manure-source AOM gross-route balance failed.');
    end
    diagnostic=struct('source_kind',context.source_kind, ...
        'aom_to_pon_g_cm3',tagPon,'aom_to_don_g_cm3',tagDon, ...
        'aom_to_tan_g_cm3',tagTan,'direct_fast_don_g_cm3',tagFast, ...
        'external_nh4_source_g_cm3',externalNH4_g_cm3, ...
        'input_g_cm2',input_g_cm2,'aom_closure_residual_g_cm2',residual, ...
        'aom_to_pon_gross_g_cm2',sum(tagPon.*dz), ...
        'aom_to_don_gross_g_cm2',sum(tagDon.*dz), ...
        'aom_to_tan_gross_g_cm2',sum(tagTan.*dz), ...
        'direct_fast_don_gross_g_cm2',sum(tagFast.*dz));
end

function value=vector(s,name,n)
    assert(isstruct(s)&&isfield(s,name),'soil:nitrogen:tracer:manureGrossField', ...
        'Required applied manure gross field %s is absent.',name);
    value=s.(name)(:);
    assert(isnumeric(value)&&isreal(value)&&numel(value)==n ...
        &&all(isfinite(value)&value>=0),'soil:nitrogen:tracer:manureGrossDomain', ...
        'Applied manure field %s must be nonnegative NL-by-1.',name);
end

function value=scalar(s,name)
    value=vector(s,name,1);
end

function fraction=donorFraction(tag,total,name)
    tolerance=512*eps(max([total(:);realmin]));
    assert(all(isfinite(tag))&&all(tag>=0&tag<=total+tolerance), ...
        'soil:nitrogen:tracer:manureDonor','Manure label exceeds %s donor.',name);
    fraction=zeros(size(total));active=total>0;
    assert(all(tag(~active)==0),'soil:nitrogen:tracer:manureZeroDonor', ...
        'Zero %s donor contains labelled N.',name);
    fraction(active)=tag(active)./total(active);
    assert(all(fraction>=0&fraction<=1+512*eps), ...
        'soil:nitrogen:tracer:manureFraction','Invalid %s label fraction.',name);
    % Preserve the calculated label mass, including bounded roundoff near
    % unity. The accepted-state bounds and unchanged column gate audit it;
    % do not silently discard a fraction of the passive source.
end

function compareRemaining(actual,total,name)
    tolerance=512*eps(max([abs(actual(:));abs(total(:));realmin]));
    assert(all(isfinite(actual)&actual>=0)&&all(abs(actual-total)<=tolerance), ...
        'soil:nitrogen:tracer:manureGrossEndpoint', ...
        'Explicit gross routes do not reproduce %s endpoint.',name);
end
