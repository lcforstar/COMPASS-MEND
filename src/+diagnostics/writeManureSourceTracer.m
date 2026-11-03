function writeManureSourceTracer(reg,built,kt)
%WRITEMANURESOURCETRACER Write only committed passive manure cohorts/ledgers.
% Internal element stocks and AOM gross amounts are bottom-up [g N cm^-3];
% exports are top-down. Step fluxes are integrated once, not rates.
% Editor: Wenjie Liu; Date: 2026-10-30.
    assert(isfield(built.Y,'soil_manure_n_tracer') ...
        &&isfield(built.workspace,'manure_n_tracer'), ...
        'diagnostics:manureTracer:missingOwner','Committed manure labels/ledgers are missing.');
    s=diagnostics.manureTracerSchema();n=built.domain.NL;dz=built.domain.dz(:);
    for cohort=1:size(s.cohorts,1)
        key=s.cohorts{cohort,1};prefix=s.cohorts{cohort,2};
        state=built.Y.soil_manure_n_tracer.(key);
        diagnostic=built.workspace.manure_n_tracer.(key);
        for j=1:size(s.nodes,1)
            value=state.(s.nodes{j,2});
            assert(isvector(value)&&numel(value)==built.domain.NN&&all(isfinite(value)));
            reg.write([prefix s.nodes{j,1}],flip(value(:)),kt);
        end
        for j=1:size(s.layers,1)
            layer(s.layers{j,1},state.organic_n.(s.layers{j,2}));
        end
        layer('don_profile',state.dom.don_total_g_cm3);
        layer('aom_manure_profile',state.aom.manure.element_N);
        scalar('aom_manure_surface',state.aom.manure.surface_N);
        for j=1:numel(s.storage)
            name=s.storage{j};scalar([name '_storage'],diagnostic.storage.([name '_g_cm2']));
        end
        for j=1:numel(s.routes)
            name=s.routes{j};scalar([name '_step'],diagnostic.step.([name '_g_cm2']));
            scalar([name '_cumulative'],state.cumulative.([name '_g_cm2']));
        end
        for j=1:size(s.aom_routes,1)
            name=s.aom_routes{j,2};value=diagnostic.manure_routes.(s.aom_routes{j,1});
            layer([name '_profile_step'],value);scalar([name '_step'],sum(value(:).*dz));
        end
        scalar('aom_closure_residual_step',diagnostic.manure_routes.aom_closure_residual_g_cm2);
        scalar('mend_closure_residual_step',diagnostic.step.mend_closure_residual_g_cm2);
        scalar('din_closure_residual_step',diagnostic.step.din_closure_residual_g_cm2);
        scalar('cumulative_closure_residual',diagnostic.cumulative_closure_residual_g_cm2);
        diagnostics.writeFertilizerTracerDetails(reg,built,kt,prefix,state,diagnostic);
    end
    function layer(name,value)
        assert(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==n ...
            &&all(isfinite(value)),'diagnostics:manureTracer:layer','Invalid accepted manure profile.');
        reg.write([prefix name],flip(value(:)),kt);
    end
    function scalar(name,value)
        assert(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value), ...
            'diagnostics:manureTracer:scalar','Invalid accepted manure ledger.');
        reg.write([prefix name],value,kt);
    end
end
