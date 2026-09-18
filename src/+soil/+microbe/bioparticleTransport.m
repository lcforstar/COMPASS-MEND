function candidate=bioparticleTransport(microbes,organicN,thetaBottom, ...
        bulkDensityBottom,dzBottom,qBottom,parameters,dt_s)
%BIOPARTICLETRANSPORT Transport active/dormant microbial C/N between phases.
%
%   Mobile biomass follows the closed-column advection-dispersion operator;
%   attached/mobile exchange and irreversible straining follow the PFLOTRAN
%   BIOPARTICLE transport-only source structure. Growth, respiration, death,
%   and dormancy are excluded. C and N use the identical linear operator and
%   are independently conserved. State units remain MEMC [mg g^-1] for C and
%   COMPASS [g cm^-3] for N.
%   The irreversible retained store follows the first-order straining term of
%   Gargiulo et al. (2007), J. Contam. Hydrol. 92, 255--273,
%   https://doi.org/10.1016/j.jconhyd.2007.01.009. PFLOTRAN BIOPARTICLE is
%   retained as the attachment/detachment transport interface reference.
%
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-08-30 Codex -- document the separate literature ownership of
%   irreversible straining and enable its nonzero site-calibration path.
% Change: 2026-09-14 Jizixing Zhu -- carry explicitly stored intracellular
%   enzyme phases with the same cell-particle C/N transport/exchange operator.

    n=numel(dzBottom); bd=bulkDensityBottom(:); dz=dzBottom(:);
    theta=thetaBottom(:); q=qBottom(:);
    assert(all([numel(bd),numel(theta),numel(q)]==n), ...
        'bioparticleTransport:shape','BIOPARTICLE driver shapes must equal NL.');
    % Change: 2026-09-11 Jizixing Zhu -- reaction owner must already close
    % phase sums. Transport cannot repair them by inventing phase transfers.
    groups={'MB','BD'};
    if isfield(microbes,'ENZN')
        assert(isfield(organicN,'ENZN_N'),'bioparticleTransport:enzymeOwner', ...
            'Intracellular enzyme C requires its paired N owner.');
        groups{end+1}='ENZN';
    end
    for key=groups
        k=key{1};
        csum=microbes.([k '_mobile'])+microbes.([k '_attached'])+microbes.([k '_retained']);
        nsum=organicN.([k '_mobile_N'])+organicN.([k '_attached_N'])+organicN.([k '_retained_N']);
        assert(all(abs(csum-microbes.(k))<=1e-11.*max(1,abs(microbes.(k)))) ...
            &&all(abs(nsum-organicN.([k '_N']))<=1e-13.*max(1,abs(organicN.([k '_N'])))), ...
            'bioparticleTransport:phaseClosure','Reaction phase sums must match C/N aggregates.');
    end
    beforeC=columnCarbon(microbes,bd,dz);
    beforeN=columnNitrogen(organicN,dz);
    diagnostics=struct();
    for prefix=groups
        name=prefix{1};
        c=struct('mobile',flip(microbes.([name '_mobile'])).*bd./1000, ...
            'attached',flip(microbes.([name '_attached'])).*bd./1000, ...
            'retained',flip(microbes.([name '_retained'])).*bd./1000);
        nn=struct('mobile',organicN.([name '_mobile_N']), ...
            'attached',organicN.([name '_attached_N']), ...
            'retained',organicN.([name '_retained_N']));
        [c,nn,d]=advancePair(c,nn,theta,bd,dz,q,parameters,dt_s);
        microbes.([name '_mobile'])=flip(c.mobile./bd.*1000);
        microbes.([name '_attached'])=flip(c.attached./bd.*1000);
        microbes.([name '_retained'])=flip(c.retained./bd.*1000);
        organicN.([name '_mobile_N'])=nn.mobile;
        organicN.([name '_attached_N'])=nn.attached;
        organicN.([name '_retained_N'])=nn.retained;
        microbes.(name)=microbes.([name '_mobile']) ...
            +microbes.([name '_attached'])+microbes.([name '_retained']);
        organicN.([name '_N'])=organicN.([name '_mobile_N']) ...
            +organicN.([name '_attached_N'])+organicN.([name '_retained_N']);
        diagnostics.(name)=d;
    end
    afterC=columnCarbon(microbes,bd,dz); afterN=columnNitrogen(organicN,dz);
    if isfield(microbes,'ENZN')
        microbes.Tot=microbes.POC+microbes.MOC+microbes.MB+microbes.BD ...
            +microbes.EP+microbes.EM+microbes.ENZN;
    end
    residualC=afterC-beforeC; residualN=afterN-beforeN;
    assert(abs(residualC)<=1e-12*max(1,beforeC) ...
        &&abs(residualN)<=1e-12*max(1,beforeN), ...
        'bioparticleTransport:closure', ...
        'BIOPARTICLE C/N transport failed closed-column conservation.');
    candidate=struct('microbes',microbes,'organic_n',organicN, ...
        'diagnostics',diagnostics,'carbon_residual_g_cm2',residualC, ...
        'nitrogen_residual_g_cm2',residualN, ...
        'reaction_terms_disabled',true,'elapsed_s',dt_s);
end

function [c,n,d]=advancePair(c,n,theta,bd,dz,q,p,dt)
% Advance one active/dormant phase pair with separate C and N ledgers.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    [c.mobile,dc]=soil.dom.transportConservative(c.mobile,theta,bd,0,dz,q, ...
        p.diffusion_cm2_s,p.dispersion_cm,dt,p.transport_substep_s);
    [n.mobile,dn]=soil.dom.transportConservative(n.mobile,theta,bd,0,dz,q, ...
        p.diffusion_cm2_s,p.dispersion_cm,dt,p.transport_substep_s);
    elapsed=0; count=0;
    while elapsed<dt-1e-12
        step=min(p.transport_substep_s,dt-elapsed);
        attach=p.attachment_rate_s.*c.mobile.*step;
        detach=p.detachment_rate_s.*c.attached.*step;
        strain=p.straining_rate_s.*c.mobile.*step;
        mobileOut=attach+strain;
        scale=min(1,c.mobile./max(mobileOut,realmin));
        attach=attach.*scale; strain=strain.*scale;
        detach=min(detach,c.attached);
        ratioMobile=safeRatio(n.mobile,c.mobile);
        ratioAttached=safeRatio(n.attached,c.attached);
        nAttach=min(n.mobile,attach.*ratioMobile);
        nStrain=min(n.mobile-nAttach,strain.*ratioMobile);
        nDetach=min(n.attached,detach.*ratioAttached);
        c.mobile=c.mobile-attach-strain+detach;
        c.attached=c.attached+attach-detach; c.retained=c.retained+strain;
        n.mobile=n.mobile-nAttach-nStrain+nDetach;
        n.attached=n.attached+nAttach-nDetach; n.retained=n.retained+nStrain;
        elapsed=elapsed+step; count=count+1;
    end
    d=struct('carbon_transport',dc,'nitrogen_transport',dn, ...
        'exchange_substeps',count,'active_motility',false);
end


function value=columnCarbon(s,bdBottom,dz)
% Active plus dormant microbial C column [g C cm^-2].
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-09-14 Jizixing Zhu -- include intracellular enzyme cargo.
    value=sum(flip(s.MB+s.BD).*bdBottom(:)./1000.*dz(:));
    if isfield(s,'ENZN'),value=value+sum(flip(s.ENZN).*bdBottom(:)./1000.*dz(:));end
end

function value=columnNitrogen(s,dz)
% Active plus dormant microbial N column [g N cm^-2].
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-09-14 Jizixing Zhu -- include intracellular enzyme cargo N.
    value=sum((s.MB_N+s.BD_N).*dz(:));
    if isfield(s,'ENZN_N'),value=value+sum(s.ENZN_N.*dz(:));end
end

function ratio=safeRatio(n,c)
% Phase N:C with exact zero-carbon behavior.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    ratio=zeros(size(c)); active=c>0; ratio(active)=n(active)./c(active);
end
