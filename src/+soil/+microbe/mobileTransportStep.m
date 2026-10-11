function step=mobileTransportStep(theta,dz,q,p,bottom)
%MOBILETRANSPORTSTEP Common positivity clock for all open mobile C/N cargo.
% theta [-], dz [cm], q [cm/s upward], bottom.q_liquid_cm_s [into-domain].
% For kd=0, each explicit FV row loses at most dt*rate(i) of its own mass.
% rate includes outgoing upwind advection, symmetric dispersion/diffusion,
% and the bottom advective outlet. dt*max(rate)<=0.9 is a sufficient
% nonnegative-coefficient condition independent of concentration or C:N.
% This changes only the open-boundary numerical clock; closed legacy paths
% keep their existing adaptive schedule. 0.9 is numerical margin, not a fit.
% Editor: Wenjie Liu; Date: 2026-09-29.
% Change: 2026-09-29 Wenjie Liu -- include the two configured external-face
% loss coefficients, including fixed-concentration half-cell dispersion.
    theta=theta(:);dz=dz(:);q=q(:);n=numel(theta);
    assert(numel(dz)==n&&numel(q)==n&&all(isfinite([theta;dz;q])) ...
        &&all(theta>0)&&all(dz>0),'bioparticle:mobileCapacity', ...
        'Suspended particles require positive liquid capacity and layer thickness.');
    alpha=p.dispersion_cm;if isscalar(alpha),alpha=repmat(alpha,n,1);else,alpha=alpha(:);end
    assert(numel(alpha)==n&&all(isfinite(alpha))&&all(alpha>=0) ...
        &&isscalar(p.diffusion_cm2_s)&&isfinite(p.diffusion_cm2_s)&&p.diffusion_cm2_s>=0 ...
        &&isscalar(p.transport_substep_s)&&isfinite(p.transport_substep_s)&&p.transport_substep_s>0, ...
        'bioparticle:mobileParameters','Invalid mobile transport coefficients.');
    assert(ismember(bottom.kind,{'advective_outflow','configured'})&&isscalar(bottom.q_liquid_cm_s) ...
        &&isfinite(bottom.q_liquid_cm_s),'bioparticle:bottomFlux','Accepted bottom flux required.');
    [~,boundaryLoss]=soil.dom.boundaryCoefficients(bottom,theta,dz,p.diffusion_cm2_s,alpha);
    rate=zeros(n,1);
    for i=1:n-1
        qf=.5*(q(i)+q(i+1));distance=.5*(dz(i)+dz(i+1));
        mix=(.5*(theta(i)+theta(i+1))*p.diffusion_cm2_s ...
            +.5*(alpha(i)+alpha(i+1))*abs(qf))/distance;
        rate(i)=rate(i)+(max(qf,0)+mix)/(theta(i)*dz(i));
        rate(i+1)=rate(i+1)+(max(-qf,0)+mix)/(theta(i+1)*dz(i+1));
    end
    rate(1)=rate(1)+boundaryLoss(1)/(theta(1)*dz(1));
    rate(end)=rate(end)+boundaryLoss(2)/(theta(end)*dz(end));
    assert(all(isfinite(rate))&&all(rate>=0),'bioparticle:mobileRate','Invalid transport loss coefficient.');
    step=p.transport_substep_s;
    if any(rate>0),step=min(step,.9/max(rate));end
end
