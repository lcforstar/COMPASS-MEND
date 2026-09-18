function [concentration, diagnostic] = colmVerticalTransport( ...
        concentration, sourceIncrement, dz_m, dt_s, advection_m_s, diffusivity_m2_s)
%COLMVERTICALTRANSPORT CoLM Patankar/tridiagonal vertical mixing port.
%
%   Inputs are top-down soil-cell concentrations and source increments in the
%   same concentration unit, dz [m], dt [s], advection [m s^-1], and diffusion
%   [m^2 s^-1]. The A(Pe)=max(0,(1-0.1|Pe|)^5), harmonic diffusivity,
%   coefficient assembly, source placement, and zero-gradient ghost equations
%   follow MOD_BGC_Soil_BiogeochemLittVertTransp.F90 verbatim in operation
%   order. Source SHA-256:
%   70E59F3BB69B40099517FED4F4B482887AA3A0A7977F973087A510122D9DC731.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
    concentration=concentration(:); sourceIncrement=sourceIncrement(:);
    dz_m=dz_m(:); n=numel(concentration);
    assert(n>=1&&numel(sourceIncrement)==n&&numel(dz_m)==n ...
        &&all(isfinite(concentration))&&all(concentration>=0) ...
        &&all(isfinite(sourceIncrement))&&all(dz_m>0)&&dt_s>0 ...
        &&isscalar(advection_m_s)&&isfinite(advection_m_s) ...
        &&isscalar(diffusivity_m2_s)&&isfinite(diffusivity_m2_s) ...
        &&diffusivity_m2_s>=0, ...
        'soil:litter:badInput','CoLM vertical-transport inputs are invalid.');
    before=sum(concentration.*dz_m);
    epsilon=1e-30;
    adv=repmat(max(abs(advection_m_s),epsilon)*signOrOne(advection_m_s),n+1,1);
    diffus=repmat(max(diffusivity_m2_s,epsilon),n+1,1);
    zi=[0;cumsum(dz_m)];
    z=zi(1:end-1)+dz_m/2;
    z=[z;z(end)+dz_m(end)/2];
    dzNode=zeros(n+1,1); dzNode(1)=z(1);
    for j=2:n+1, dzNode(j)=z(j)-z(j-1); end
    dm=zeros(n+1,1); dp=zeros(n+1,1);
    fm=zeros(n+1,1); fp=zeros(n+1,1);
    pem=zeros(n+1,1); pep=zeros(n+1,1);
    for j=1:n+1
        if j==1
            w=(z(j+1)-zi(j+1))/dzNode(j+1);
            dh=harmonic(diffus(j),diffus(j+1),w);
            dp(j)=dh/dzNode(j+1); fm(j)=adv(j); fp(j)=adv(j+1);
            pem(j)=0; pep(j)=fp(j)/dp(j);
        elseif j>=n+1
            w=(zi(j)-z(j-1))/dzNode(j);
            dh=harmonic(diffus(j),diffus(j-1),w);
            dm(j)=dh/dzNode(j); dp(j)=dm(j);
            fm(j)=adv(j); fp(j)=0;
            pem(j)=fm(j)/dm(j); pep(j)=fp(j)/dp(j);
        else
            wm=(zi(j)-z(j-1))/dzNode(j);
            wp=(z(j+1)-zi(j+1))/dzNode(j+1);
            dm(j)=harmonic(diffus(j),diffus(j-1),wm)/dzNode(j);
            dp(j)=harmonic(diffus(j),diffus(j+1),wp)/dzNode(j+1);
            fm(j)=adv(j); fp(j)=adv(j+1);
            pem(j)=fm(j)/dm(j); pep(j)=fp(j)/dp(j);
        end
    end
    a=zeros(n+2,1); b=zeros(n+2,1); c=zeros(n+2,1); r=zeros(n+2,1);
    old=[0;concentration;0];
    b(1)=1; c(1)=-1;
    for j=1:n
        idx=j+1; ap0=dz_m(j)/dt_s;
        a(idx)=-(dm(j)*A(pem(j))+max(fm(j),0));
        c(idx)=-(dp(j)*A(pep(j))+max(-fp(j),0));
        b(idx)=-a(idx)-c(idx)+ap0;
        if j==1
            r(idx)=sourceIncrement(j)*dz_m(j)/dt_s ...
                +(ap0-adv(j))*old(idx);
        else
            r(idx)=sourceIncrement(j)*dz_m(j)/dt_s+ap0*old(idx);
        end
    end
    a(end)=-1; b(end)=1;
    solved=tridiagonal(a,b,c,r);
    concentration=solved(2:end-1);
    if any(~isfinite(concentration))||any(concentration<-1e-12)
        error('soil:litter:invalidSolution', ...
            'CoLM vertical transport produced invalid concentrations.');
    end
    concentration=max(concentration,0);
    after=sum(concentration.*dz_m);
    diagnostic=struct('storage_before',before, ...
        'storage_after',after,'storage_change',after-before, ...
        'source_hash', ...
        '70E59F3BB69B40099517FED4F4B482887AA3A0A7977F973087A510122D9DC731');
end

function value=A(pe)
    value=max(0,(1-0.1*abs(pe))^5);
end

function value=harmonic(left,right,weight)
    if left>0&&right>0
        value=1/((1-weight)/left+weight/right);
    else
        value=0;
    end
end

function s=signOrOne(value)
    s=sign(value); if s==0, s=1; end
end

function x=tridiagonal(a,b,c,r)
    n=numel(b); cp=zeros(n,1); rp=zeros(n,1);
    pivot=b(1);
    assert(abs(pivot)>eps,'soil:litter:singularMatrix','Singular CoLM matrix.');
    cp(1)=c(1)/pivot; rp(1)=r(1)/pivot;
    for i=2:n
        pivot=b(i)-a(i)*cp(i-1);
        assert(abs(pivot)>eps,'soil:litter:singularMatrix','Singular CoLM matrix.');
        if i<n, cp(i)=c(i)/pivot; end
        rp(i)=(r(i)-a(i)*rp(i-1))/pivot;
    end
    x=zeros(n,1); x(n)=rp(n);
    for i=n-1:-1:1, x(i)=rp(i)-cp(i)*x(i+1); end
end
