function total=totalFromFree(free,theta,rho,kd)
%TOTALFROMFREE Reconstruct linear-equilibrium D+Q from an explicit D target.
% All columns use the same bottom-up element grid. free and total are
% g C/cm3 bulk soil, theta is liquid volume fraction, rho is g dry soil/cm3,
% and kd is cm3/g (scalar or column). Q=D*rho*Kd/theta; no mass is created
% by this conversion: the initializer must fund total from its carbon budget.
% Positive free mass at zero liquid water is rejected, never water-floored.
% Zero free mass in a dry cell selects zero total; an existing sorbed store
% must instead be preserved explicitly by the caller.
% Editor: Wenjie Liu; Date: 2026-10-12.
    free=free(:);theta=theta(:);rho=rho(:);
    if isscalar(kd),kd=repmat(kd,size(free));else,kd=kd(:);end
    assert(isequal(size(free),size(theta),size(rho),size(kd)) ...
        &&all(isfinite([free;theta;rho;kd]))&&all(free>=0) ...
        &&all(theta>=0)&&all(rho>0)&&all(kd>=0), ...
        'soil:dom:freeInitialValues','Free-D target and phase drivers must be finite and aligned.');
    wet=theta>0;
    assert(all(free(~wet)==0),'soil:dom:freeInitialDry', ...
        'Positive initial free DOC requires positive liquid-water volume.');
    total=zeros(size(free));
    total(wet)=free(wet).*(1+rho(wet).*kd(wet)./theta(wet));
    assert(all(isfinite(total)),'soil:dom:freeInitialOverflow','Initial D+Q is not finite.');
end
