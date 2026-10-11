function model = initialize(siteDir)
%INITIALIZE  Top-level entry: build a ready-to-run model from a site folder.
%
%   model = init.initialize(siteDir) is the single facade for model setup. It
%   runs the full initialization chain and returns ONE task container (a model
%   struct), never a bare built struct and never via assignin('base',...):
%
%       config = io.loadSiteConfig(siteDir)  read site context and switches
%       plan   = init.resolveSimulationPlan(config) select active components
%       raw    = io.loadSite(siteDir, plan)  read only active component inputs
%       report = io.checkSite(raw, plan)     validate only active inputs
%       built  = init.buildState(raw) or init.buildWaterOnlyState(raw, plan)
%       built  = init.allocateWorkspace(built)  per-step scratch + diag buffers
%       built  = init.initializeSoilAux(built)  pre-step soil auxiliary refresh
%
%   checkSite errors are NOT swallowed: a missing required column/file aborts
%   initialize with the aggregated checkSite:invalid error, so model is only
%   returned for a site that passed validation (model.check_report.ok == true).
%
%   Inputs:
%     siteDir  [char]  path to a site folder, e.g. 'data/sites/changwu'
%
%   Returns model (one simulation task; the main loop consumes model.built):
%     model.siteDir       [char]      provenance (raw.siteDir)
%     model.raw           [struct]    raw site containers (kept for debug/audit)
%     model.check_report  [struct]    io.checkSite outcome (.ok == true)
%     model.built         [struct]    domain/control/boundary/sources/forcing/
%                                     parameters/vegetation_state/Y/p/workspace/diagnostics
%     model.status        [string]    "initialized"
%     model.created_at    [datetime]  initialization time
%
%   This facade only sequences single-writer builders. The pre-step soil
%   auxiliary refresh is intentionally separate from input loading and state
%   construction. See docs/init_refactor_progress.md.
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-11 Jizixing Zhu -- route field and laboratory sites through an explicit component input plan.
% Change: 2026-07-18 Jizixing Zhu -- align an initial bottom temperature node with
%         its same-time Dirichlet boundary and rebuild U-dependent auxiliaries.
% Change: 2026-09-29 Wenjie Liu -- condition warm field water and total DIN
%         with the actual accepted endpoint/control-volume initialization.

    config = io.loadSiteConfig(siteDir);
    plan = init.resolveSimulationPlan(config);
    raw = io.loadSite(siteDir, plan, config);
    report = io.checkSite(raw, plan);         % hard errors propagate (not swallowed)
    if plan.is_laboratory
        built = init.buildWaterOnlyState(raw, plan);
    else
        built = init.buildState(raw);
        built.control.simulation = plan;
    end
    built  = init.allocateWorkspace(built);
    built  = init.initializeSoilAux(built);
    built  = init.alignInitialThermalBoundary(built);
    [built,measuredAudit] = init.conditionMeasuredSoilInitialState(built,raw);
    built.p.initialization.measured_soil=measuredAudit;

    model              = struct();
    model.siteDir      = raw.siteDir;
    model.raw          = raw;
    model.check_report = report;
    model.built        = built;
    model.status       = "initialized";
    model.created_at   = datetime("now");
end
