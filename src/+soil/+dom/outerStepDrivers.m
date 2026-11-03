function [theta,q,dt,faces]=outerStepDrivers(soilAdvance,thetaFinal,dtOuter,liquidDriver,parameters)
%OUTERSTEPDRIVERS Use the same outer accepted liquid QL/QMB as inorganic N.
% thetaFinal: bottom-up NL-by-1 end-of-hour liquid content; liquidDriver.QL:
% bottom-up NL-by-1 accepted-dt mean [cm/s], upward positive. dtOuter [s].
% Bottom QMB [cm/s, into column positive] comes from the water owner's
% accumulated outer mean; top uses QL(end), matching the DIN node adapter.
% Returns one immutable driver segment, not a replay of water solver tiles.
% The DOM FV owner still substeps at its configured positivity-limited dt.
% Physical open-bottom support remains liquid-only free drainage. No change
% to the water or inorganic-N solve. Missing/nonfinite drivers reject.
% Editor: Wenjie Liu; Date: 2026-10-17.
    theta=thetaFinal(:);q=liquidDriver.QL(:);dt=dtOuter;faces={};
    assert(isequal(size(theta),size(q))&&all(isfinite([theta;q])) ...
        &&all(theta>=0)&&isscalar(dt)&&isfinite(dt)&&dt>0, ...
        'soil:dom:outerDriver','Invalid end-of-step theta or outer liquid-flux driver.');
    if ~isfield(parameters,'boundary')||strcmp(parameters.boundary.bottom_kind,'closed'),return;end
    assert(isfield(soilAdvance,'water_flux_accum') ...
        &&isfield(soilAdvance.water_flux_accum,'QMB'), ...
        'soil:dom:missingOuterBottom','Open DOM requires the outer water QMB.');
    qb=soilAdvance.water_flux_accum.QMB;
    assert(isnumeric(qb)&&isreal(qb)&&isscalar(qb)&&isfinite(qb)&&qb<=0, ...
        'soil:dom:outerBottom','Physical DOM bottom currently requires liquid-only free drainage.');
    faces={struct('schema_version',1,'kind','advective_outflow', ...
        'phase_contract','free_drainage_liquid_only_zero_vapor', ...
        'position','physical_column_bottom','positive_direction','into_domain', ...
        'q_liquid_cm_s',qb,'q_vapor_cm_s',0,'q_top_cm_s',q(end),'dt_s',dt)};
end
