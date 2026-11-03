function [candidate,report]=accelerateSlowPools(y0,y1,y2,maxFactor)
%ACCELERATESLOWPOOLS Bounded periodic Aitken update for POC/MOC initialization.
%
%   Only POC and MOC are replaced. Their paired PON/MON stores receive the
%   identical multiplicative update, preserving each element's current C:N.
%   An update is admissible only where two annual increments have the same
%   sign and contract. This is a COMPASS adaptation of CoLM-SASU's spin-up-
%   controller ownership principle; it is not CoLM's flux-matrix inversion.
%
%   Editor: Jizixing Zhu; Date: 2026-09-02.
% Change: 2026-09-14 Codex -- keep authoritative POM1/POM2 and PON1/PON2
%   subpools synchronized when the optional POC acceleration is enabled.
    assert(isscalar(maxFactor)&&isfinite(maxFactor)&&maxFactor>=1, ...
        'accelerateSlowPools:factor','maxFactor must be finite and >= 1.');
    candidate=y2; count=0; factors=[];
    totalDelta=zeros(size(y2.soil_microbes.Tot));
    pairs={{'POC','PON'},{'MOC','MON'}};
    for i=1:numel(pairs)
        cName=pairs{i}{1}; nName=pairs{i}{2};
        x0=y0.soil_microbes.(cName)(:);
        x1=y1.soil_microbes.(cName)(:);
        x2=y2.soil_microbes.(cName)(:);
        d1=x1-x0; d2=x2-x1;
        ratio=zeros(size(x2)); active=d1.*d2>0 & abs(d2)<abs(d1) & abs(d1)>0;
        ratio(active)=d2(active)./d1(active);
        factor=ones(size(x2));
        factor(active)=min(maxFactor,1./max(1-ratio(active),eps));
        target=max(0,x2+(factor-1).*d2);
        scale=ones(size(x2)); nonzero=x2>0;
        scale(nonzero)=target(nonzero)./x2(nonzero);
        scale(~nonzero&target>0)=1; % no paired N basis: reject creation
        target(~nonzero)=x2(~nonzero);
        candidate.soil_microbes.(cName)=reshape(target, ...
            size(y2.soil_microbes.(cName)));
        totalDelta=totalDelta+reshape(target-x2,size(totalDelta));
        n=y2.soil_organic.(nName)(:);
        candidate.soil_organic.(nName)=reshape(n.*flip(scale), ...
            size(y2.soil_organic.(nName)));
        if strcmp(cName,'POC')
            for key={'POM1','POM2'}
                name=key{1};
                candidate.soil_microbes.(name)= ...
                    y2.soil_microbes.(name).*reshape(scale, ...
                    size(y2.soil_microbes.(name)));
            end
            for key={'PON1','PON2'}
                name=key{1};
                candidate.soil_organic.(name)= ...
                    y2.soil_organic.(name).*reshape(flip(scale), ...
                    size(y2.soil_organic.(name)));
            end
            candidate.soil_microbes.POC= ...
                candidate.soil_microbes.POM1+candidate.soil_microbes.POM2;
            candidate.soil_organic.PON= ...
                candidate.soil_organic.PON1+candidate.soil_organic.PON2;
        end
        count=count+nnz(abs(target-x2)>0);
        factors=[factors;factor(active)]; %#ok<AGROW>
    end
    candidate.soil_microbes.Tot=y2.soil_microbes.Tot+totalDelta;
    report=struct('applied',count>0,'accelerated_elements',count, ...
        'method','bounded_periodic_aitken_slow_pools', ...
        'maximum_factor',max([1;factors]), ...
        'paired_nitrogen_scaled',true,'production_equations_modified',false, ...
        'exact_colm_sasu',false);
end
