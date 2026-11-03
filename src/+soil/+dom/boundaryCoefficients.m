function [input,loss]=boundaryCoefficients(b,theta,dz,diffusion,alpha,rewetting)
%BOUNDARYCOEFFICIENTS Directed FV boundary exchange coefficients [cm/s].
% Outputs are [bottom;top]. input is a mass rate [g/cm2/s], loss multiplies
% the adjacent aqueous concentration [g/cm3 water]. Type 1 imposes a FACE
% concentration through half-cell dispersion/diffusion plus upwind advection;
% it does not overwrite the prognostic whole-cell total or sorbed storage.
% Type 2 top has prescribed input and no evaporative solute/cell loss;
% bottom type 2 has prescribed inflow and advective (zero-gradient) outflow.
% Type 3 bottom has zero gradient on outflow; reverse flow requires type 2.
% q at both ends is positive upward. Positive input/loss are into/out of soil.
% Editor: Wenjie Liu; Date: 2026-09-29.
    input=zeros(2,1);loss=input;
    if nargin<6,rewetting=false(size(theta));end
    if strcmp(b.kind,'closed'),return;end
    if strcmp(b.kind,'advective_outflow')
        assert(isfinite(b.q_liquid_cm_s)&&b.q_liquid_cm_s<=0, ...
            'soil:dom:bottomInflowUnsupported','Bottom inflow requires external concentration.');
        loss(1)=-b.q_liquid_cm_s;return
    end
    assert(strcmp(b.kind,'configured')&&isfinite(b.q_liquid_cm_s) ...
        &&isfinite(b.q_top_cm_s),'soil:dom:transportBoundary','Invalid organic boundary driver.');
    q=[b.q_liquid_cm_s;-b.q_top_cm_s];ends=[1 numel(theta)];
    if isscalar(alpha),alpha=repmat(alpha,numel(theta),1);end
    sides={'bottom','top'};
    for i=1:2
        f=b.(sides{i});idx=ends(i);c=f.concentration_g_cm3_water;
        valid=[1 2 3];if i==2,valid=[1 2];end
        assert(isscalar(f.type)&&ismember(f.type,valid)&&isscalar(c)&&isfinite(c)&&c>=0, ...
            'soil:dom:transportBoundary','Invalid boundary type or aqueous concentration.');
        if f.type==1
            conductance=(theta(idx)*diffusion+alpha(idx)*abs(q(i)))/(.5*dz(idx));
            input(i)=(max(q(i),0)+conductance)*c;
            loss(i)=max(-q(i),0)+conductance;
        elseif f.type==2
            input(i)=max(q(i),0)*c;
            if i==1,loss(i)=max(-q(i),0);end
        else
            assert(q(i)<=0&&c==0,'soil:dom:bottomInflowUnsupported', ...
                'Outflow-only type 3 cannot supply reverse-flow concentration; select type 2.');
            loss(i)=-q(i);
        end
        assert((input(i)==0&&loss(i)==0)||theta(idx)>0||(rewetting(idx)&&q(i)>0), ...
            'soil:dom:dryBoundary', ...
            'Active liquid boundary requires positive liquid water content.');
    end
end
