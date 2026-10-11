function built = allocateWorkspace(built)
%ALLOCATEWORKSPACE  Pre-allocate the per-step solver workspace + diagnostic buffers.
%
%   built = init.allocateWorkspace(built) adds two single-writer containers to a
%   built state produced by init.buildState:
%       built.workspace    time-step scratch (trial/previous state, derived cache,
%                          solver matrices, flux scratch)
%       built.diagnostics  loose legacy output buffers (the registry owns the rest)
%
%   It is a single-writer for those two containers ONLY. It does NOT touch
%   built.Y (accepted prognostic state, owned by buildInitialState/buildPlant*),
%   built.parameters, or built.p. Everything here is zero/empty-initialised
%   scratch (recomputed before use), except the trial/previous step state, which
%   is warm-started from built.Y so the first Picard guess and the first temporal
%   term (hh-h)/dt are well posed.
%
%   Container layout (see docs/allocate_workspace_audit.md):
%     workspace.step.trial      {h,T,p_gas,U,nh4,no3,urea,h_frez,soil_gas}  warm copy of Y
%     workspace.step.previous   prior accepted (h^{n-1}, extrapolation) + backup{} (rollback)
%     workspace.derived         Se/KL_h/RHOV/Theta_LL/*BAR/D*DZ/Retard/... (recomputable)
%     workspace.solver.coeff    PDE coefficient bank Chh..Ccaa + N coeffs R1..G2
%     workspace.solver.{water,thermal,gas,nitrogen}  assembled tridiagonal systems
%     workspace.solver.SAVE     3x3x5 shared boundary-flux snapshot
%     workspace.flux            Srt/SinkNh4/Smin.../V_L/Dsh_* (per-step sink/velocity)
%     workspace.surface         DSTMAX/NBChh/IRPT1/IRPT2 + hcrit/h_target/hN (ponding cap + dynamic BC code + retry + top water-BC scratch)
%     diagnostics.{plantN_diag,nsc_diag,NH4_total,NO3_total,CHK2,CHK3}
%
%   Legacy provenance: legacy StateVariables / legacy SolverCoefficients /
%   legacy FluxVariables constructors + setInitialNitrogen blocks (2)(3).
%   The legacy h/hh/hOLD map to Y.h / workspace.step.trial.h /
%   workspace.step.previous.h (see audit doc, headline finding).
%
%   Inputs:
%     built  struct  init.buildState output; uses built.domain (NN/NL/nCell,
%                    grid.nD/grid.nS), built.control.time.tS, built.Y
%                    (warm start). For nCell>1, node fields are [NN x nCell]
%                    and element fields are [NL x nD x nCell].
%   Returns built with built.workspace + built.diagnostics added.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-06-26 Jizixing Zhu -- add freeze/N/gas derived scratch
%     (T_CRIT/TT_CRIT node; eT/eTheta/ePH/D_Vg edge) so workspace.derived
%     matches the legacy StateVariables fields the soil kernel recomputes
%     each step. Dead legacy fields (alpha_h/bx/b1/b11/ALPHA) intentionally
%     NOT migrated (no compute consumer; see docs/allocate_workspace_audit.md).
%   Change: 2026-06-28 Jizixing Zhu — add workspace.surface (DSTMAX/NBChh/IRPT1/IRPT2).
%   Change: 2026-07-01 Jizixing Zhu — add surface top water-BC scratch (hcrit/h_target/hN) to the fixed field set; runOneStep only UPDATES them (no dynamic field creation — Julia fixed-struct friendly).
%   Change: 2026-07-05 Jizixing Zhu — add built.diagnostics.canopy_water (8 tS x nCell history buffers, moved from runtime.canopy_water; state-ownership migration increment 3 stage 1).
%   Change: 2026-07-07 Jizixing Zhu — move DRHODAz from nodeDerived (NN x 1) to edgeDerived (NL x 1): it is a layer-pair gradient produced 1xNL by Density_DA, consumed by layer index (calculateEnergyParameters:198). B2 shape-correction.

%   Change: 2026-07-11 Jizixing Zhu -- preallocate named phase-storage and
%           freezing-head partial derivatives for the conservative water Jacobian.
%   Change: 2026-07-12 Jizixing Zhu -- add workspace.boundary.thermal for
%           recomputable Tbh/Tb runtime boundary values.
%   Change: 2026-07-13 Jizixing Zhu -- predeclare accepted hourly soil-heat
%           boundary fluxes for the new scalar heat solver.
%   Change: 2026-07-13 Jizixing Zhu -- rename the derived gas pressure
%           mobility field to Kg without changing its values.
%   Change: 2026-07-14 Jizixing Zhu -- predeclare wetting-heat derivatives
%           dW_dh and dW_dTheta in the recomputable soil auxiliary cache.
%   Change: 2026-07-17 Jizixing Zhu -- source the surface-water capacity from the
%           normalized site boundary instead of hard-coding zero.
%   Change: 2026-07-17 Jizixing Zhu -- predeclare the ClimaLand snow adapter
%           scratch, accepted snow diagnostics, and frozen-pond transfer.
%   Change: 2026-07-17 Jizixing Zhu -- document hcrit as compatibility reference;
%           residual-water state now owns automatic dry-set activation.
%   Change: 2026-07-20 Jizixing Zhu -- allocate soil-mode candidate and trial-only
%     dual-porosity exchange auxiliary outside accepted Y.

    Y  = built.Y;
    NN = built.domain.NN;
    NL = built.domain.NL;
    nD = getDim(built.domain.grid, 'nD', 2);
    nS = getDim(built.domain.grid, 'nS', 2);
    nCell = getDomainDim(built.domain, 'nCell', 1);
    tS = built.control.time.tS;
    validateStateShape(Y, NN, nCell);

    nodeSize = {NN, nCell};
    if nCell == 1
        elemSize   = {NL, nD};
        edgeSize   = {NL, 1};
        soluteSize = {NL, nD, nS};
    else
        elemSize   = {NL, nD, nCell};
        edgeSize   = {NL, nCell};
        soluteSize = {NL, nD, nS, nCell};
    end

    ws = struct();

    % ===================================================================
    %  step state: trial (t_{n+1} candidate) + previous (t_{n-1} + backups)
    %  warm-started from the accepted state in built.Y (Y owns t_n).
    % ===================================================================
    z = zeros(nodeSize{:});
    trial = struct('h', Y.h, 'T', Y.T, 'p_gas', Y.p_gas, 'U', Y.U, ...
                   'nh4', Y.nh4, 'no3', Y.no3, 'urea', z, 'h_frez', z, ...
                   'C_NH4_granule', Y.C_NH4_granule, ...
                   'canopy_air', Y.canopy_air, 'scope_soil', Y.scope_soil, ...
                   'soil_gas', Y.soil_gas, 'soil_modes', Y.soil_modes, ...
                   'surface_water', Y.surface_water);

    previous = struct('h', Y.h, 'T', Y.T, 'p_gas', Y.p_gas, 'U', Y.U, 'h_frez', z, ...
                      'nh4', Y.nh4, 'no3', Y.no3, 'urea', z, ...
                      'C_NH4_granule', Y.C_NH4_granule, ...
                      'canopy_air', Y.canopy_air, 'scope_soil', Y.scope_soil, ...
                      'soil_gas', Y.soil_gas, 'soil_modes', Y.soil_modes, ...
                      'surface_water', Y.surface_water);
    % rollback snapshots: filled on substep entry, empty until then
    previous.backup = struct('h', [], 'T', [], 'p_gas', [], ...
                             'nh4', [], 'no3', [], 'urea', [], 'C_NH4_granule', [], ...
                             'canopy_air', [], 'scope_soil', [], 'soil_gas', []);

    % L_f: per-step freeze latent heat [J/kg] (set by runOneStep.prepareStepState
    % from lFreezing + LFUSION_OVERRIDE); a fixed field so runOneStep only updates it.
    trialAux = struct( ...
        'dual_porosity', struct( ...
            'SinkIm_current', zeros(elemSize{:})), ...
        'surface_water', struct('active', false));
    ws.step = struct('trial', trial, 'trial_aux', trialAux, ...
        'previous', previous, 'L_f', 0);

    % scope_step: per-step SCOPE/ebal input bundle (rebuilt each outer step by
    % runOneStep.prepareCanopyStep via init.prepareScopeStep); empty until stepped.
    ws.scope_step = struct();

    % Runtime boundary auxiliaries are separate from static built.boundary
    % parameters and from equation coefficients/matrices.
    ws.boundary = struct('thermal', struct( ...
        'Tbh', zeros(tS, nCell), 'Tb', [], ...
        'Tb_jd_series', [], 'Tb_JD0', []));

    % ===================================================================
    %  derived cache: recomputable each step from state (Se, KL_h, *BAR, ...)
    % ===================================================================
    nodeDerived = {'RHODA','RHOV','RHOV_s','DRHOV_sT','DRHOVh','DRHOVT', ...
        'DRHODAt','Xaa','XaT','Xah','HR','D_A','L','LL','DhT', ...
        'Gamma_hh','COR','CORh','cc','DhFrezDh','DhFrezDT', ...
        'T_CRIT','TT_CRIT'};   % freezing-point depression (HT_frez -> KfL_T), NN x 1
    elemDerived = {'Se','KL_h','KfL_h','KfL_T','KL_T','DTheta_LLh','DTheta_UUh', ...
        'Ratio_ice','W','WW','dW_dh','dW_dTheta','Lambda_eff','c_unsat','EHCAP','ETCON','ZETA', ...
        'MU_W','f0','L_WT','D_Ta','Theta_g','Theta_V','D_V','k_g','Sa', ...
        'Alpha_Lg','POR_C','Eta','Kg','Theta_LL','Theta_UU','Theta_II', ...
        'Theta_U','Theta_I','Retard','fw','DshNh4','DshNo3','DifNh4','DifNo3', ...
        'Vpore','WFPS','CRootNh4','CRootNo3', ...
        'DallAmicoThermalTime','DallAmicoThetaLStar','DallAmicoPhaseSource', ...
        'DallAmicoLimiterCode', ...
        'DThetaLiquidDh','DThetaLiquidDT','DThetaWaterEquivDh','DThetaWaterEquivDT', ...
        'DThetaIceDh','DThetaIceDT','DThetaGasDh','DThetaGasDT'};
    edgeDerived = {'KLhBAR','KLTBAR','RHOVBAR','EtaBAR','Lambda_effBAR', ...
        'c_unsatBAR','Theta_LLBAR','DTDBAR','DETBAR','DEhBAR','DRHOVhDz', ...
        'DRHOVTDz','DIFBAR_NH4','DIFBAR_NO3','E1BAR','E2BAR','RetardBAR', ...
        'DhDZ','DTDZ','DPgDZ','DNH4DZ','DNO3DZ', ...
        'eT','eTheta','ePH', ...   % N rate temp/moisture/pH modifiers (solveNitrogen), NL x 1
        'D_Vg', ...                % gas-phase longitudinal dispersivity (matric coeff), NL x 1
        'DRHODAz'};                % dry-air density space gradient (Density_DA, layer-pair), NL x 1
    d = struct();
    d = allocZeros(d, nodeDerived, nodeSize{:});
    d = allocZeros(d, elemDerived, elemSize{:});
    d = allocZeros(d, edgeDerived, edgeSize{:});
    phys = equations.physicalConstants();
    d.Gamma_hh(:) = 1;                  % [-] no freezing partition before first update
    d.T_CRIT(:)   = phys.T_freeze;      % [K] legacy StateVariables ice-point default
    d.TT_CRIT(:)  = phys.T_freeze;      % [K] legacy StateVariables ice-point default
    ws.derived = d;

    % ===================================================================
    %  solver: PDE coefficient bank (.coeff) + per-equation assembled systems
    % ===================================================================
    coeffWater   = {'Chh','ChT','Khh','KhT','Kha','Vvh','VvT','Chg'};
    coeffThermal = {'CTh','CTT','CTa','KTh','KTT','KTa','VTh','VTT','VTa','CTg', ...
                    'Kcvh','KcvT','Kcva','Ccvh','CcvT'};
    coeffAir     = {'Cah','CaT','Caa','Kah','KaT','Kaa','Vah','VaT','Vaa','Cag', ...
                    'Kcah','KcaT','Kcaa','Ccah','CcaT','Ccaa'};
    coeffN       = {'R1','R1O','R2','R2O','E1','E2','B1','B2','F1','F2','G1','G2', ...
                    'C1O','C2O'};
    coeff = struct();
    coeff = allocZeros(coeff, [coeffWater, coeffThermal, coeffAir, coeffN], elemSize{:});

    solver = struct();
    solver.coeff    = coeff;
    solver.water    = assembledSystem(nodeSize, elemSize);
    solver.thermal  = assembledSystem(nodeSize, elemSize);
    solver.gas      = assembledSystem(nodeSize, elemSize);
    solver.nitrogen = assembledSystem(nodeSize, elemSize);
    solver.nitrogen.N_SAVE = zeros(3, 3, 5);   % N transport boundary-flux snapshot
    solver.SAVE = zeros(3, 3, 5);              % shared water/thermal/gas boundary-flux snapshot
    ws.solver = solver;

    % ===================================================================
    %  flux scratch: per-step sinks (root/N uptake), reaction rates, velocities
    % ===================================================================
    elemFlux  = {'Srt','SinkNh4','SinkNo3','AUptakeA','Smin','Snit','SVnit', ...
        'Sden','SVden','Svol','SVvol','Shys','V_L','V_A','Sdis'};
    soluteFlux = {'SinkS','Dsh_h','Dsh_T','Dsh_A'};   % NL x nD x nS
    fx = struct();
    fx = allocZeros(fx, elemFlux, elemSize{:});
    fx = allocZeros(fx, soluteFlux, soluteSize{:});
    % PR-E surface water fluxes (scalar per outer step; snow-adjusted throughfall + soil EVAP feed the
    % ④ soil-water top BC, CanopyInterception is diagnostic). Committed accept-only by stepCanopyWaterSnowEvap.
    fx.Throughfall        = zeros(1, nCell);   % [cm/s] snow-adjusted liquid water to ground (top water BC)
    fx.CanopyInterception = zeros(1, nCell);   % [cm/s] wet-canopy evaporation rate (diagnostic)
    fx.EVAP               = zeros(1, nCell);   % [cm/s] soil evaporation (top water BC)
    fx.heat_hourly = struct('QET', zeros(1, nCell), ...
        'QEB', zeros(1, nCell));               % [J cm-2 s-1] accepted top/bottom soil-heat flux
    fx.soil_energy_balance = soil.thermal.emptyColumnEnergyBalance();
    fx.soil_dry_air_mass_balance = soil.gas.emptyColumnMassBalance();
    fx.soil_boundary_hourly = soil.boundary.emptyFluxSummary();
    ws.flux = fx;

    % ===================================================================
    %  surface-water runtime scratch (legacy setBoundaryCondition.m:30-115):
    %  per-step ponding cap (DSTMAX), dynamic top water-head BC code (NBChh), and
    %  Picard substep retry flags. NBChh is flux(2) unless the top is Dirichlet(1).
    % ===================================================================
    NBChh0 = 2;
    if isfield(built, 'boundary') && isfield(built.boundary, 'NBCh_top') ...
            && built.boundary.NBCh_top == 1
        NBChh0 = 1;
    end
    % hcrit is a legacy finite dry-head reference, not an automatic switch.
    % The typed atmospheric active set enters its residual-dry branch only
    % after the converged top endpoint reaches residual water, and preserves
    % that state's actual head. h_target/hN remain fixed-shape kernel scratch.
    surfaceCapacity = 0;
    if isfield(built, 'boundary') && isfield(built.boundary, 'surface_water') ...
            && isfield(built.boundary.surface_water, 'capacity_cm')
        surfaceCapacity = built.boundary.surface_water.capacity_cm;
    end
    ws.surface = struct('DSTMAX', surfaceCapacity, 'NBChh', NBChh0, 'IRPT1', 0, 'IRPT2', 0, ...
                        'hcrit', -1.5e4, 'h_target', 0, 'hN', 0, ...
                        'frozen_transfer_cm', 0);
    ws.snow = struct('phase', struct(), 'flux', struct());

    built.workspace = ws;

    % ===================================================================
    %  diagnostics: only the LOOSE legacy buffers (registry owns the ~185 rest)
    %  init values match the legacy code (setInitialNitrogen blocks (2)(3)).
    % ===================================================================
    diag = struct();
    if nCell == 1
        diag.NH4_total = zeros(NN, tS);          % setInitialNitrogen.m:105
        diag.NO3_total = zeros(NN, tS);          % setInitialNitrogen.m:106
        diag.CHK2      = inf(NN, 1);             % NH4 Picard residual
        diag.CHK3      = inf(NN, 1);             % NO3 Picard residual
    else
        diag.NH4_total = zeros(NN, tS, nCell);
        diag.NO3_total = zeros(NN, tS, nCell);
        diag.CHK2      = inf(NN, nCell);
        diag.CHK3      = inf(NN, nCell);
    end
    diag.plantN_diag = plantNDiagTemplate(tS);   % setInitialNitrogen.m:193-200
    diag.nsc_diag    = nscDiagTemplate(tS);      % setInitialNitrogen.m:303-312
    % PR-D prereq (step 1b): plant C/N + root-uptake diagnostic slots pre-declared so the reused
    % runPlantCNStep / root-uptake glue never grow dynamic fields (Julia fixed-struct safety).
    diag.turnover_diag = turnoverDiagTemplate(tS);   % setInitialNitrogen.m:210-212 (tS x 1, 4 fields [])
    diag.frPHU_diag    = zeros(tS, 1);               % setInitialNitrogen.m:213 [-] frPHU snapshot per KT
    diag.root_demand   = struct('nh4', zeros(tS, 1), 'no3', zeros(tS, 1));   % per-KT root N demand cache (glue -> reg.write)
    diag.root_uptake_M1_diag = struct();             % current-KT DAISY+Nacry uptake diag (replaced each accept; empty when Nefc~=1)
    diag.v2_hydraulic = struct();                    % PR-A current-KT V2 dynamic-hydraulics diag (W/psi/RWC; replaced each accept; empty when scheme~=dynamic)
    % canopy-water diagnostic histories (moved from runtime.canopy_water 2026-07-05; row kt written
    % accept-only in commitStepPending). tS x nCell loose buffers; registry/flush wiring deferred (stage 5).
    diag.canopy_water = struct( ...
        'h2ocan_history',       zeros(tS, nCell), 'h2ocan_sun_history',    zeros(tS, nCell), ...
        'h2ocan_sha_history',   zeros(tS, nCell), 'fwet_history',          zeros(tS, nCell), ...
        'fwet_sun_history',     zeros(tS, nCell), 'fwet_sha_history',      zeros(tS, nCell), ...
        'lE_CI_history',        zeros(tS, nCell), 'canopy_water_residual', zeros(tS, nCell));
    diag.snow = struct();
    built.diagnostics = diag;
end

% ---------------------------------------------------------------------------
function s = allocZeros(s, names, varargin)
% Add each name in NAMES to struct S as zeros(varargin{:}).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    for i = 1:numel(names)
        s.(names{i}) = zeros(varargin{:});
    end
end

% ---------------------------------------------------------------------------
function sys = assembledSystem(nodeSize, elemSize)
% One equation's assembled tridiagonal workspace (legacy SolverCoefficients form):
%   element coefficients C1..C6 and assembled global C3_a..C9/RHS.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    sys = struct();
    sys = allocZeros(sys, {'C1','C2','C3','C4','C5','C6'}, elemSize{:});
    sys = allocZeros(sys, {'C3_a','C4_a','C5_a','C6_a','C7','C9','RHS','RHS1'}, nodeSize{:});
end

% ---------------------------------------------------------------------------
function v = getDim(gridCfg, name, default)
% Read a grid dimension with a fallback (nD/nS default to 2).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if isfield(gridCfg, name) && ~isempty(gridCfg.(name))
        v = gridCfg.(name);
    else
        v = default;
    end
end

% ---------------------------------------------------------------------------
function v = getDomainDim(domain, name, default)
% Read a domain dimension with a fallback (nCell defaults to one column).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if isfield(domain, name) && ~isempty(domain.(name))
        v = domain.(name);
    else
        v = default;
    end
end

% ---------------------------------------------------------------------------
function validateStateShape(Y, NN, nCell)
% Guard regional/global contracts before allocating matching workspace arrays.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    nodeFields = {'h','T','p_gas','U','nh4','no3'};
    expected = [NN, nCell];
    for i = 1:numel(nodeFields)
        name = nodeFields{i};
        if ~isfield(Y, name)
            error('allocateWorkspace:missingStateField', 'built.Y.%s is required.', name);
        end
        actual = size(Y.(name));
        if numel(actual) > 2 || ~isequal(actual(1:2), expected)
            error('allocateWorkspace:badStateShape', ...
                'built.Y.%s must be [%d x %d], got [%s].', ...
                name, expected(1), expected(2), strjoin(string(actual), ' x '));
        end
    end
    if ~isfield(Y, 'soil_gas') || ~isfield(Y.soil_gas, 'C_O2')
        error('allocateWorkspace:missingStateField', 'built.Y.soil_gas.C_O2 is required.');
    end
    actual = size(Y.soil_gas.C_O2);
    if numel(actual) > 2 || ~isequal(actual(1:2), expected)
        error('allocateWorkspace:badStateShape', ...
            'built.Y.soil_gas.C_O2 must be [%d x %d], got [%s].', ...
            expected(1), expected(2), strjoin(string(actual), ' x '));
    end
end

% ---------------------------------------------------------------------------
function p = plantNDiagTemplate(tS)
% tS x 1 struct array; all fields [] so the registry detects KT=1 via .NPP empty.
% Field list matches setInitialNitrogen.m:193-200 (computeNPP/allocatePlantN writers) PLUS A_croot:
% runPlantCNStep.m:109 writes plantN_diag(KT).A_croot, which the legacy template omits (grown
% dynamically). Pre-declared here (fixed-struct safety for the Julia port).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-04 Jizixing Zhu -- PR-D step 1b: pre-declare A_croot (was grown dynamically by runPlantCNStep, legacy template omits it).
    p = repmat(struct('NPP',[],'NPP_net',[],'R_a_plant',[],'R_m_leaf',[], ...
        'R_m_root',[],'R_m_stem',[],'R_m_livecroot',[],'R_m_fruit',[], ...
        'r_m_leaf_actual',[],'r_m_root_actual',[],'r_m_stem_actual',[], ...
        'r_m_livecroot_actual',[],'r_m_fruit_actual',[],'r_g_actual',[], ...
        'ra_above_actual',[],'ra_below_actual',[], ...
        'R_g',[],'f_N_growth',[],'A_leaf',[],'A_root',[],'A_stem',[],'A_croot',[],'A_fruit',[]), ...
        tS, 1);
end

% ---------------------------------------------------------------------------
function n = nscDiagTemplate(tS)
% tS x 1 struct array; numeric defaults (resp_beta_sub/resp_theta_actual = 1).
% Field list + defaults MUST match setInitialNitrogen.m:303-312 (applyNSCCoupling writer).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-04 Jizixing Zhu -- add NPP_actual_kg (applyNSCCoupling.m:304 writes it; template
%   had drifted, so runPlantCNStep's whole-struct nsc_diag(KT)=nsc_diag hit heterogeneousStrucAssignment).
    template = struct('R_a_potential_g', 0, 'R_a_actual_g', 0, ...
        'R_a_actual_kg', 0, 'NPP_actual_kg', 0, 'storage_drawdown_event', 0, 'nsc_total_g', 0, 'nsc_net_g', 0, ...
        'c_storage_drawdown_g', 0, ...
        'A_C_leaf', 0, 'A_C_root', 0, 'A_C_stem', 0, 'A_C_fruit', 0, ...
        'C_to_root_storage', 0, 'C_to_stem_storage', 0, ...
        'C_xsmr', 0, 'R_m_xs_g', 0, 'xsmr_recover_g', 0, ...
        'resp_beta_sub', 1, ...
        'resp_substrate_g', 0, 'C_labile_target_g', 0, 'C_storage_to_labile', 0, ...
        'resp_theta_actual', 1);
    n = repmat(template, tS, 1);
end

% ---------------------------------------------------------------------------
function t = turnoverDiagTemplate(tS)
% tS x 1 struct array; all fields [] (setInitialNitrogen.m:210-212). Written per KT by
% runPlantCNStep (turnover + harvest daily events).
% Editor: Jizixing Zhu; Date: 2026-07-04.
    t = repmat(struct('harvest_today', [], 'leaf_litter_today', [], ...
        'root_litter_today', [], 'storage_release_today', []), tS, 1);
end

% Change: 2026-09-30 Wenjie Liu -- predeclare annual stem WSC and
% standing-senescence applied diagnostics for transactional C/N updates.
% tS x 1 struct array; numeric defaults (resp_beta_sub/resp_theta_actual = 1).
% Field list + defaults MUST match setInitialNitrogen.m:303-312 (applyNSCCoupling writer).
% Change: 2026-09-30 Wenjie Liu -- withdraw this round's APSIM stem WSC
%   and standing-senescent-leaf extension; retain the original COMPASS litter
%   and CoLM crop C/N allocation/harvest contracts. Earlier Change records are historical.
