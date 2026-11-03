function [tiles,hourly]=waterDriverSchedule(advance,thetaFinal,dtOuter,liquidDriver,parameters)
%WATERDRIVERSCHEDULE Read the marked production hydraulic schedule for DOM.
% Older adapters supplied only outer QL/QMB and a liquid endpoint. Preserve
% their exact prescribed constant driver without inventing a starting theta.
% The returned outer summaries retain the original NL-by-1 and scalar shapes.
    sd=advance.solver_diagnostics;
    v2=isfield(sd,'dom_water_driver_version');
    hasSchedule=isfield(sd,'accepted_transport_schedule');
    if v2
        assert(isequal(sd.dom_water_driver_version,2)&&hasSchedule, ...
            'soil:dom:waterDriverVersion','Invalid or missing v2 DOM water schedule.');
        s=sd.accepted_transport_schedule;
        fields={'schema_version','tile_count','elapsed_s','dt_s','theta_liquid_initial','theta_liquid_endpoint','QL'};
        assert(isstruct(s)&&isscalar(s)&&all(isfield(s,fields))&&isequal(s.schema_version,1), ...
            'soil:dom:waterDriverSchema','Incomplete v2 DOM water schedule.');
        dt=s.dt_s(:).';theta=s.theta_liquid_endpoint;q=s.QL;initial=s.theta_liquid_initial;
        n=numel(thetaFinal); expected=[n,numel(dt)];
        assert(~isempty(dt)&&all(isfinite(dt)&dt>0)&&isequal(s.tile_count,numel(dt)) ...
            &&isfinite(s.elapsed_s)&&abs(sum(dt)-dtOuter)<=1e-9*max(1,dtOuter) ...
            &&abs(s.elapsed_s-sum(dt))<=1e-9*max(1,dtOuter) ...
            &&isequal(size(q),expected)&&isequal(size(theta),expected)&&isequal(size(initial),[n 1]) ...
            &&all(isfinite(q(:)))&&all(isfinite(theta(:))&theta(:)>=0)&&all(isfinite(initial)&initial>=0) ...
            &&isequal(theta(:,end),thetaFinal(:)), ...
            'soil:dom:waterDriverCoverage','Invalid DOM water storage, flux or elapsed time.');
        faces=soil.dom.acceptedBottomDrivers(advance,parameters,dt);
        policy='accepted_water_tiles_substep_start_theta_v2';
    else
        assert(~hasSchedule || (isstruct(sd.accepted_transport_schedule)&&isscalar(sd.accepted_transport_schedule) ...
            &&isfield(sd.accepted_transport_schedule,'schema_version')&&isequal(sd.accepted_transport_schedule.schema_version,1) ...
            &&~isfield(sd.accepted_transport_schedule,'theta_liquid_initial')), ...
            'soil:dom:waterDriverVersion','A start/end DOM water schedule lacks its producer version.');
        policy='din_outer_mean_endpoint_v1';
    end
    [endTheta,meanQ,duration,meanFaces]=soil.dom.outerStepDrivers(advance,thetaFinal,dtOuter,liquidDriver,parameters);
    hourly=struct('policy',policy,'theta',endTheta,'QL',meanQ,'dt_s',duration,'bottom_faces',{meanFaces});
    if v2
        tolerance=32*numel(dt)*eps(max([abs(q(:));abs(meanQ);realmin]));
        assert(all(abs(q*dt.'/dtOuter-meanQ)<=tolerance),'soil:dom:waterDriverMean', ...
            'The chronological fluxes must agree with the independently supplied accepted outer mean.');
        tiles=struct('policy',policy,'theta_initial',initial,'theta',theta,'QL',q,'dt_s',dt,'bottom_faces',{faces});
    else
        tiles=hourly;
        tiles.theta_initial=[];
    end
end
