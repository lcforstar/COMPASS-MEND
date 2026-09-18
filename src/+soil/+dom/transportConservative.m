function [totalNext, diagnostic] = transportConservative( ...
        total, theta, bulkDensity, kd, dz_cm, q_cm_s, ...
        diffusion_cm2_s, dispersion_cm, dt_s, maxSubstep_s)
%TRANSPORTCONSERVATIVE Closed-boundary finite-volume DOM ADE transport.
%
%   Element order is COMPASS bottom-up. q_cm_s is the accepted liquid flux at
%   each element centre, positive deeper-to-shallower. Internal-face advection
%   is upwind; dispersion/diffusion is centred. Only internal faces exchange
%   mass, so sum(total.*dz) is invariant to roundoff.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-09-02 Codex -- accept depth-resolved Kd and dispersivity;
%   internal-face dispersivity is the arithmetic mean of adjacent elements.
%   Change: 2026-09-08 Jizixing Zhu -- expose accepted per-layer DOM transport
%     change [g/cm3 bulk soil] for independent reaction/transport budget audits.
%   Change: 2026-09-14 Jizixing Zhu -- reject negative trial inventories
%     instead of clipping them; absolute clipping breaks paired C/N transport.
    total=total(:); theta=theta(:); bulkDensity=bulkDensity(:);
    dz_cm=dz_cm(:); q_cm_s=q_cm_s(:);
    n=numel(total);
    kd=expandProfile(kd,n,'kd');
    dispersion_cm=expandProfile(dispersion_cm,n,'dispersion_cm');
    assert(n>=1 && numel(theta)==n && numel(bulkDensity)==n ...
        && numel(dz_cm)==n && numel(q_cm_s)==n, ...
        'soil:dom:transportSize','DOM transport arrays must have NL values.');
    assert(all(isfinite(total))&&all(total>=0)&&all(theta>=0) ...
        && all(bulkDensity>0)&&all(dz_cm>0)&&all(isfinite(q_cm_s)) ...
        && diffusion_cm2_s>=0&&all(dispersion_cm>=0) ...
        &&dt_s>0&&maxSubstep_s>0, ...
        'soil:dom:transportValues','DOM transport inputs are invalid.');
    beforeProfile=total;
    before=sum(total.*dz_cm);
    elapsed=0; acceptedSubsteps=0; minSubstep=Inf;
    while elapsed < dt_s-1e-12
        trialDt=min(maxSubstep_s,dt_s-elapsed);
        while true
            phase=soil.dom.equilibrate(total,theta,bulkDensity,kd);
            concentration=phase.aqueous_g_cm3_water;
            flux=zeros(max(n-1,0),1);
            for face=1:n-1
                qFace=0.5*(q_cm_s(face)+q_cm_s(face+1));
                if qFace>=0, upwind=concentration(face); ...
                else, upwind=concentration(face+1); end
                thetaFace=max(0.5*(theta(face)+theta(face+1)),0);
                distance=0.5*(dz_cm(face)+dz_cm(face+1));
                dispersionFace=0.5*(dispersion_cm(face)+dispersion_cm(face+1));
                dispersive=thetaFace*diffusion_cm2_s ...
                    +dispersionFace*abs(qFace);
                flux(face)=qFace*upwind ...
                    +dispersive*(concentration(face)-concentration(face+1))/distance;
            end
            candidate=total;
            if n>1
                transfer=flux*trialDt;
                candidate(1:end-1)=candidate(1:end-1)-transfer./dz_cm(1:end-1);
                candidate(2:end)=candidate(2:end)+transfer./dz_cm(2:end);
            end
            if all(candidate>=0) && all(isfinite(candidate))
                total=candidate;
                break;
            end
            trialDt=trialDt/2;
            if trialDt<1e-6
                error('soil:dom:transportNonnegative', ...
                    'DOM transport could not find a nonnegative substep.');
            end
        end
        elapsed=elapsed+trialDt;
        acceptedSubsteps=acceptedSubsteps+1;
        minSubstep=min(minSubstep,trialDt);
    end
    after=sum(total.*dz_cm);
    totalNext=total;
    diagnostic=struct('net_layer_change_g_cm3',total-beforeProfile, ...
        'storage_before_g_cm2',before, ...
        'storage_after_g_cm2',after, ...
        'closure_residual_g_cm2',after-before, ...
        'accepted_substeps',acceptedSubsteps, ...
        'min_substep_s',minSubstep, ...
        'boundary_condition','closed_zero_mass_flux');
end

function value=expandProfile(value,n,name)
% Normalize a scalar or element profile without changing its physical value.
    assert(isnumeric(value)&&all(isfinite(value(:)))&&all(value(:)>=0) ...
        &&(isscalar(value)||numel(value)==n), ...
        'soil:dom:transportProfile', ...
        '%s must be a finite nonnegative scalar or NL-element profile.',name);
    if isscalar(value), value=repmat(double(value),n,1); ...
    else, value=double(value(:)); end
end
