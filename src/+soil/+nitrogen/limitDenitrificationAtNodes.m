function [limited,info]=limitDenitrificationAtNodes(rate,otherRate,cOld,thetaOld,dz,dt)
%LIMITDENITRIFICATIONATNODES Bound denitrification by the assembled nodal RHS.
% Bottom-up endpoint rates [ug N/cm3 soil/s], old nodal liquid concentration
% [ug N/cm3 water], old endpoint liquid capacity, layer dz [cm], dt [s].
% OTHER_RATE is .98*Snit-root_NO3+organic_NO3: unchanged external exchanges.
% The existing layer-inventory cap must run first. This additional cap uses
% the SAME lumped old storage and consistent linear source weights as
% transport.assembleMatrixCoefficients/assembleEquation. Both endpoints of
% an element receive the minimum adjacent nodal allowance factor, preserving
% nonnegative rates and never increasing any layer's allocated reaction.
% A negative allowance caused by OTHER_RATE is not repaired or redistributed.
% No state commit, clipping or change to transport/boundary equations occurs.
% The caller must propagate reduced actual Sden into its reaction and gas
% ledgers. This is an availability condition, not a general ADE positivity
% proof for arbitrary advection matrices or boundary conditions.
% Editor: Jizixing Zhu; Date: 2026-09-15.
    dz=dz(:);cOld=cOld(:);L=numel(dz);
    assert(isequal(size(rate),[L,2],size(otherRate),size(thetaOld)) ...
        &&numel(cOld)==L+1&&isscalar(dt)&&isfinite(dt)&&dt>0 ...
        &&all(isfinite([rate(:);otherRate(:);thetaOld(:);cOld;dz])) ...
        &&all(rate>=0,'all')&&all(thetaOld>=0,'all')&&all(cOld>=0)&&all(dz>0), ...
        'soil:nitrogen:nodalDenitrificationInputs','Invalid nodal reaction allocation inputs.');
    load=@(v)[(2*v(:,1)+v(:,2)).*dz/6;0]+[0;(v(:,1)+2*v(:,2)).*dz/6];
    stock=[thetaOld(:,1).*dz/2;0]+[0;thetaOld(:,2).*dz/2];
    available=stock.*cOld/dt+load(otherRate);
    demand=load(rate);factor=ones(L+1,1);active=demand>max(available,0);
    factor(active)=max(available(active),0)./demand(active)*(1-8*eps);
    element=min(factor(1:end-1),factor(2:end));limited=rate.*element;
    info=struct('node_available_ug_cm2_s',available,'node_requested_ug_cm2_s',demand, ...
        'node_applied_ug_cm2_s',load(limited),'element_scale',element);
end
