function registerManureSourceTracer(reg,group)
%REGISTERMANURESOURCETRACER Register independent organic-N/TAN mass labels.
% No physical field or legacy mineral-fertilizer variable is renamed. All
% manure label diagnostics retain double precision in the NetCDF adapter.
% Editor: Wenjie Liu; Date: 2026-10-30.
    s=diagnostics.manureTracerSchema();
    for cohort=1:size(s.cohorts,1)
        key=s.cohorts{cohort,1};prefix=s.cohorts{cohort,2};label=s.cohorts{cohort,3};
        owner=['Y.soil_manure_n_tracer.' key];
        work=['workspace.manure_n_tracer.' key];
        for j=1:size(s.nodes,1)
            add(s.nodes{j,1},'ug N cm-3',{'time','node'},[owner '.' s.nodes{j,2}]);
        end
        for j=1:size(s.layers,1)
            add(s.layers{j,1},'g N cm-3',{'time','layer'},[owner '.organic_n.' s.layers{j,2}]);
        end
        add('don_profile','g N cm-3',{'time','layer'},[owner '.dom.don_total_g_cm3']);
        add('aom_manure_profile','g N cm-3',{'time','layer'},[owner '.aom.manure.element_N']);
        add('aom_manure_surface','g N cm-2',{'time'},[owner '.aom.manure.surface_N']);
        for j=1:numel(s.storage)
            name=s.storage{j};add([name '_storage'],'g N cm-2',{'time'},[work '.storage.' name '_g_cm2']);
        end
        for j=1:numel(s.routes)
            name=s.routes{j};
            add([name '_step'],'g N cm-2',{'time'},[work '.step.' name '_g_cm2']);
            add([name '_cumulative'],'g N cm-2',{'time'},[owner '.cumulative.' name '_g_cm2']);
        end
        for j=1:size(s.aom_routes,1)
            name=s.aom_routes{j,2};source=[work '.manure_routes.' s.aom_routes{j,1}];
            add([name '_profile_step'],'g N cm-3',{'time','layer'},source);
            add([name '_step'],'g N cm-2',{'time'},['column integral of ' source]);
        end
        add('aom_closure_residual_step','g N cm-2',{'time'},[work '.manure_routes.aom_closure_residual_g_cm2']);
        add('mend_closure_residual_step','g N cm-2',{'time'},[work '.step.mend_closure_residual_g_cm2']);
        add('din_closure_residual_step','g N cm-2',{'time'},[work '.step.din_closure_residual_g_cm2']);
        reg.register(struct('name',[prefix 'cumulative_closure_residual'], ...
            'long_name',[label ' cumulative column closure residual'], ...
            'units','g N cm-2','dims',{{'time'}},'group','conservation', ...
            'source',[work '.cumulative_closure_residual_g_cm2'], ...
            'source_code','src/+soil/+nitrogen/+tracer/advanceCandidate.m', ...
            'transform','scalar'));
        diagnostics.groups.registerFertilizerTracerDetails(reg,group,prefix,owner,work,label);
    end
    function add(name,units,dims,source)
        reg.register(struct('name',[prefix name],'long_name',[label ': ' name], ...
            'units',units,'dims',{dims},'group',group,'source',source, ...
            'source_code','src/+diagnostics/writeManureSourceTracer.m','transform','none', ...
            'comment','Passive exclusive source cohort; profiles surface to bottom. Source input field retains internal fertilizer_input key; initial TAN is a separate cohort.'));
    end
end
