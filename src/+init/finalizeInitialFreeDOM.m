function built=finalizeInitialFreeDOM(built)
%FINALIZEINITIALFREEDOM Fund a pending free-D target at final t0 liquid water.
% Called once by initialize after thermal alignment and measured water/DIN
% conditioning. The request is bottom-up mg C/kg dry soil; D+Q is reconstructed
% at mean(final Theta_LL,2) and funded only from native MOC, preserving total
% SOC in every element. MEND carbon stores are top-down mg C/g dry soil.
% No time step, physical flux, DON/N adjustment or water/energy change occurs.
% A deferred scaffold cannot be returned as a completed free-D initialization.
% Repeated calls after completion are exact no-ops. Frozen free-D starts remain
% unsupported, and an infeasible MOC budget is rejected without clipping.
% Editor: Wenjie Liu; Date: 2026-10-12.
    if ~isfield(built,'initialization') ...
            ||~isfield(built.initialization,'soil_dom') ...
            ||~isfield(built.initialization.soil_dom,'free_doc_pending') ...
            ||~built.initialization.soil_dom.free_doc_pending
        return;
    end
    assert(all(built.Y.soil_dom.doc_total_g_cm3==0), ...
        'finalizeInitialFreeDOM:notPending','Only an unfunded fresh DOC request may be finalized.');
    assert(isfield(built,'workspace')&&isfield(built.workspace,'derived') ...
        &&isfield(built.workspace.derived,'Theta_LL'), ...
        'finalizeInitialFreeDOM:missingWater','Final liquid-water endpoints are required.');
    if isfield(built.control.switches,'lFreezing')&&built.control.switches.lFreezing==1
        assert(all(isfinite(built.Y.T(:)))&&all(built.Y.T(:)>0) ...
            &&all(built.workspace.derived.Theta_II(:)==0), ...
            'buildSoilDOM:frozenFreeInitial', ...
            'doc_free_mgkg requires an unfrozen (>0 C) final initial state.');
    end
    rho=built.parameters.soil.bulk_density(:);
    theta=mean(built.workspace.derived.Theta_LL,2);
    target=built.initialization.soil_dom.free_doc_target_mgkg_bottom_up;
    free=target.*rho*1e-6;
    total=soil.dom.totalFromFree(free,theta,rho,built.parameters.soil_dom.kd_doc_cm3_g);
    built.Y.soil_dom.doc_total_g_cm3=total;
    if strcmp(built.parameters.soil_dom.sorption_model,'mend_dynamic')
        built.Y.soil_dom.doc_sorbed_g_cm3=total-free;
    end
    built=init.reserveInitialDOM(built,false);
    if isfield(built.Y,'soil_microbes')
        c=flip(total./rho)*1000;
        assert(all(c<=built.Y.soil_microbes.MOC), ...
            'reserveInitialDOM:carbonBudget','Initial DOC exceeds the accepted MOC donor.');
        built.Y.soil_microbes.MOC=built.Y.soil_microbes.MOC-c;
        built.Y.soil_microbes.Tot=built.Y.soil_microbes.Tot-c;
    end
    phase=soil.dom.phaseState(built.Y.soil_dom,'doc',theta,rho,built.parameters.soil_dom);
    errorMgkg=max(abs(phase.free_g_cm3_bulk./rho*1e6-target));
    assert(errorMgkg<=1e-9,'finalizeInitialFreeDOM:roundtrip', ...
        'Final accepted free DOC does not match the requested initial profile.');
    built.initialization.soil_dom.free_doc_pending=false;
    built.initialization.soil_dom.free_doc_finalized=true;
    built.initialization.soil_dom.free_doc_max_abs_error_mgkg=errorMgkg;
    built.initialization.soil_dom.free_doc_liquid_fraction_bottom_up=theta;
    built.initialization.soil_dom.physical_flux_added=false;
end
