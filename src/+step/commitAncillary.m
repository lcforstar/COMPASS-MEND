function model = commitAncillary(model, ledger, ancillary, kt)
%COMMITANCILLARY Commit accepted flux records, diagnostics, and management.
%
%   Inputs:
%     model   [inout] Model after all physical-domain state candidates commit.
%     ledger    [struct] Validated typed inter-domain exchange ledger.
%     ancillary [struct] Accepted diagnostic and management records.
%     kt      [-] Accepted one-based outer-step index.
%
%   Returns:
%     model [struct] Model with accepted write-only records and management.
%
%   This routine owns no prognostic physical state.  It writes accepted
%   diagnostics/output arrays from the transaction's explicit ancillary view;
%   it must never be called for a rejected outer transaction.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
%   Change: 2026-07-20 Jizixing Zhu -- remove generic flux input and make every
%     legacy workspace write a one-way projection of a typed accepted ledger.
    water = ledger.surface_to_soil.water_request;
    rootWater = ledger.vegetation_to_soil.water_request;
    rootNitrogen = ledger.vegetation_to_soil.nitrogen_request;
    model.built.workspace.flux.Throughfall = water.Throughfall;
    model.built.workspace.flux.CanopyInterception = water.CanopyInterception;
    model.built.workspace.flux.EVAP = water.EVAP;
    model.built.workspace.flux.Srt = rootWater.element_endpoint_rate_1_s;
    model.built.workspace.flux.SinkS = cat( ...
        3, rootNitrogen.sink_nh4, rootNitrogen.sink_no3);
    carbon = ledger.vegetation_to_atmosphere.carbon_result;
    if ~isempty(carbon.NEE)
        model.built.workspace.scope_step.fluxes.NEE = carbon.NEE;
    end
    if ~isempty(carbon.Resp)
        model.built.workspace.scope_step.fluxes.Resp = carbon.Resp;
    end
    if ~isempty(ancillary.diagnostics.v2_hydraulic)
        model.built.diagnostics.v2_hydraulic = ancillary.diagnostics.v2_hydraulic;
    end
    rootDemand = ancillary.diagnostics.root_demand;
    model.built.diagnostics.root_demand.nh4(kt) = rootDemand.nh4;
    model.built.diagnostics.root_demand.no3(kt) = rootDemand.no3;
    model.built.diagnostics.root_uptake_M1_diag = ...
        ancillary.diagnostics.root_uptake_M1_diag;
    if ~isempty(ancillary.diagnostics.canopy_water)
        historyRow = ancillary.diagnostics.canopy_water;
        canopyDiag = model.built.diagnostics.canopy_water;
        canopyDiag.h2ocan_history(kt, :) = historyRow.h2ocan;
        canopyDiag.h2ocan_sun_history(kt, :) = historyRow.h2ocan_sun;
        canopyDiag.h2ocan_sha_history(kt, :) = historyRow.h2ocan_sha;
        canopyDiag.fwet_history(kt, :) = historyRow.fwet;
        canopyDiag.fwet_sun_history(kt, :) = historyRow.fwet_sun;
        canopyDiag.fwet_sha_history(kt, :) = historyRow.fwet_sha;
        canopyDiag.lE_CI_history(kt, :) = historyRow.lE_CI;
        canopyDiag.canopy_water_residual(kt, :) = historyRow.residual;
        model.built.diagnostics.canopy_water = canopyDiag;
    end
    model.built.diagnostics.snow = ancillary.diagnostics.snow;
    model.runtime.management.fertilizer.applied = ...
        ancillary.management.fertilizer_applied;
end
