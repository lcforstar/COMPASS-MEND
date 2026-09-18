function [NitrogenVariables, SoilVariables, NitrogenFlux, NitrogenDiagnostics] = solveNitrogen( ...
        NitrogenVariables, NitrogenFlux, NitrogenBoundary, SoilVariables, ...
        AirVariables, GasVariables, cfg, KT, Delt_t, ExternalInputs)
%SOLVENITROGEN  Dispatcher for one outer-step N reaction + transport update.
%
% Top-level entry of +soil/+nitrogen/ package. Each Picard iteration calls
% +chemistry, +organic, +transport submodules; physics lives in submodules.
%
% Pipeline per Picard iter:
%   +chemistry/  calculateRetardationFactors → ureaHydrolysis → mineralization
%                ammoniumVolatilization → (O2 setup) → nitrification → denitrification
%   +transport/  soluteDispersion → solveTransportStep (assemble + BC + Thomas solve)
%
% Change: 2026-09-15 Jizixing Zhu -- bound denitrification by nodal FE
%   availability and propagate the reduced reaction to allocation/gas ledgers.
% Root N uptake handled by +vegetation/+root/applyRootUptake → +nitrogen/+root_uptake/nitrogenUptake.
%
% Inputs / Outputs: see body. Reference: HYDRUS-1D SOLUTE.f90,
% Zhu Jizixing thesis Ch 3.3, STEMMUS-N legacy.
%
% Refactored 2026-05-13: monolithic solveNitrogenEquations split into
% +chemistry/, +organic/, +transport/, +root_uptake/ subpackages
% (Skill §4 module organization, DAISY-inspired N decomposition).
% Editor: Jizixing Zhu; Date: 2026-07-20.
% Change: 2026-07-20 Jizixing Zhu — treat fertilizer as an explicit step
%   input while preserving the pre-event accepted storage time level.
% Change: 2026-07-20 Jizixing Zhu — remove temporary console ledger dumps; accepted
%   mass-balance diagnostics remain available through NitrogenDiagnostics.
% Change: 2026-08-16 Codex -- use signed DAISY NH4/NO3 exchange instead of
%   empirical mineralization when the native soil-carbon model is enabled.
% Change: 2026-09-07 Jizixing Zhu -- pass layer thickness to conservative
%   shared-node urea depletion and retain small positive urea stocks.
% Change: 2026-09-14 Codex -- hold reaction inventory and O2 old state at the
%   same physical step start; align DIN layer order and accepted gas products.
% Change: 2026-09-14 Codex -- reserve feasible soil respiration, limit actual
%   nitrification by remaining O2, and retain a conservative oxygen ledger.

if nargin < 10 || isempty(ExternalInputs)
    ExternalInputs = struct();
end
ExternalInputs = normalizeExternalInputs(ExternalInputs);

NL    = cfg.Grid.NL;
NN    = cfg.Grid.NN;
DeltZ = cfg.Grid.DeltZ;
nS    = cfg.Grid.nS;
if ~ExternalInputs.organic_enabled
    ExternalInputs.organic_nh4_rate = zeros(NL,1);
    ExternalInputs.organic_no3_rate = zeros(NL,1);
    ExternalInputs.organic_co2_c_rate = zeros(NL,1);
else
    assert(numel(ExternalInputs.organic_nh4_rate)==NL, ...
        'soil:nitrogen:organicExchangeSize', ...
        'Native organic exchange must contain one value per COMPASS element.');
end

NitrogenDiagnostics = struct();

% 2026-05-14 v11: N 模块自己维护 Theta_LL_for_N_prev 作为 ADE 的真旧步 θ
%
% 根因: 主循环 substep 末 SoilVariables.Theta_L = Theta_LL, 让 ADE 内部 θ_L = θ_LL
%   (两者数值相同). ADE 看不到 KT 内 θ 变化 → 不会 freeze-concentrate C → 解出
%   C_new ≈ C_old → cvno3 = θ_KT末 × C ≈ θ_KT末 × C_old, 而上 KT cvno3 用的是
%   θ_KT-1末 × C_old. 两者相比 cvno3 跳变 ~200 μg/cm² (春融期 θ 跳变大时).
%
% 修复: N 模块自己存上 KT 末 Theta_LL, 下个 KT 进来作 ADE 真旧步 θ. 这样 ADE
% 内部看到 θ_old ≠ θ_LL_new, 才能正确算 freeze-concentration.

% 临时覆盖 SoilVariables.Theta_L 为 N 模块自己维护的快照
if isfield(SoilVariables, 'Theta_LL_for_N_prev') && ~isempty(SoilVariables.Theta_LL_for_N_prev)
    SoilVariables.Theta_L = SoilVariables.Theta_LL_for_N_prev;  % 真"上 KT 末"
end
% (首次 KT=1 fallback: Theta_L 还是 substep 末值, IC 状态影响很小)

% P0a fix (2026-05-13): subtract ice volume in freezing periods so solute
% transport uses unfrozen water only. See effectiveLiquidWaterForSolute.m doc.
[Theta_Solute_LL, Theta_Solute_L] = ...
    soil.nitrogen.transport.effectiveLiquidWaterForSolute(SoilVariables, cfg);

% Preserve the management increment separately from accepted t^n.  The
% freeze-thaw partition operates on accepted storage; the explicit fertilizer
% source is installed only after that old storage has been measured.
fert_delta_nh4 = SoilVariables.C_NH4_new - SoilVariables.C_NH4;
fert_delta_no3 = SoilVariables.C_NO3_new - SoilVariables.C_NO3;
fert_delta_urea = SoilVariables.C_Urea_new - SoilVariables.C_Urea;
assert(all(fert_delta_nh4 >= -1e-12) && all(fert_delta_no3 >= -1e-12) ...
        && all(fert_delta_urea >= -1e-12), ...
    'soil:nitrogen:negativeManagementIncrement', ...
    'Staged fertilizer increments must not remove accepted nitrogen mass.');

% P0a.5: freeze-thaw N partitioning between liquid + ice.
%   Ice grows (Δθ_ice > 0) → entrap K_segr × C_liquid into ice
%   Ice melts (Δθ_ice < 0) → release ice-N back to liquid
% Mass-conservative; ice-phase N is prognostic state (C_NH4_ice/C_NO3_ice).
%
% 2026-05-14 bug fix history:
%   v1 (failed): used SoilVariables.Theta_I — but Theta_I never updated in main loop.
%   v2 (failed): used SoilVariables.Theta_II_prev — but main loop overwrites it
%       (Theta_II_prev ← Theta_II) inside the water-substep hysteresis block at
%       COMPASS.m:904, BEFORE N module fires at COMPASS.m:1712. So by the time N
%       gets it, Theta_II_prev = Theta_II → Δθ_ice = 0.
%   v3 (current): N module maintains its OWN prev-snapshot field Theta_II_for_N_prev.
%       Read it at top of solveNitrogen (= last-KT Theta_II), write it at bottom
%       (= this-KT Theta_II for next call). Decoupled from hysteresis pipeline.
% 2026-05-14 v6 setup: 冻融 N partition 已移到 Picard 后 (cv 计算前)
% 这里不调用 freezeThawPartitioning, 只初始化占位 partition_diag
% (NitrogenDiagnostics.partition 由 ADE 后那次 partition 调用真正赋值)

PH1 = NitrogenVariables.pH;
if cfg.lNtotal == 1
    for ML = 1:NL
        MN = ML;
        PH1(ML) = 0.202 * exp(-(NitrogenVariables.NO3_total(MN, KT) ./ NitrogenVariables.RHO_bulk(MN)) / 165.591) + ...
            0.208 * exp(-(NitrogenVariables.NO3_total(MN, KT) ./ NitrogenVariables.RHO_bulk(MN)) / 165.562) + 7.872;
    end
end

error_thresholds = [10, 200];
error_tolerances = [0.1, 0.1, 0.1];
max_iter = 100;
tolerance_factor = 1;

% 2026-05-14 v10: 强制 cvno3OLD(KT) = cvno3(KT-1)
% 关键洞察: 单步 dcv = cvno3(KT) - cvno3OLD(KT) ≈ 0 (test16 数据), 但累计 res = -303.
%   说明 cvno3OLD(KT) ≠ cvno3(KT-1), 之间损失"幽灵 N" 累积成 303.
%
% 物理原因: cvno3OLD 用 Theta_L × C_NH4_old, 而 Theta_L 在水量 substep 末被覆盖
% 成 Theta_LL_KT (本 KT 末水量, 不是"上 KT 末"). 同时 C_NH4 = 上 KT 末 C_new.
% 两者时步不匹配 → cvno3OLD 跟上 KT cvno3 不连续, 累计成残差.
%
% 修复: 不重算 cvno3OLD, 直接用上 KT 末记录的 cvno3 (KT=1 例外, 用 IC 一次性算).
if KT == 1 || isempty(NitrogenFlux.cvno3) || NitrogenFlux.cvno3(KT-1) == 0
    % 首次: 用 IC 算 — 用 Theta_Solute_L (= "上一步"θ, KT=1 时即 t=0 初值 θ)
    % 不能用 Theta_Solute_LL: 后者是 KT=1 水量解后的θ, 与初值 C_NO3 时步不匹配,
    % 会把 KT=1 内水量变化造成的"θ×C"跳跃记为 cv_init,产生 ~158 μg/cm² 幽灵残差
    % 冰相: IC C_NH4_ice/C_NO3_ice = 0, Theta_II_old = 0, 项为零 (KT=1 也调一次保完整)
    Theta_II_old_for_cv = zeros(NL, 2);
    if isfield(SoilVariables, 'Theta_II_for_N_prev') && ~isempty(SoilVariables.Theta_II_for_N_prev)
        Theta_II_old_for_cv = SoilVariables.Theta_II_for_N_prev;
    end
    % [2026-05-17 dissolution pool] 加 C_NH4_granule 到 NH4 储量 (KT=1 时 = 0)
    granule_KT1 = SoilVariables.C_NH4_granule;
    if isempty(granule_KT1), granule_KT1 = zeros(NN, 1); end
    [NitrogenFlux.cvnh4OLD(KT), NitrogenFlux.cvno3OLD(KT)] = ...
        soil.nitrogen.transport.calculateStorage( ...
            Theta_Solute_L, SoilVariables.C_NH4, SoilVariables.C_NO3, ...
            NitrogenVariables.RHOKG, DeltZ, NL, ...
            Theta_II_old_for_cv, SoilVariables.C_NH4_ice, SoilVariables.C_NO3_ice, ...
            granule_KT1);
else
    % 后续 KT: cvno3OLD = 上 KT 末 cvno3 (强制守恒账本连续, 包括冰相)
    NitrogenFlux.cvnh4OLD(KT) = NitrogenFlux.cvnh4(KT-1);
    NitrogenFlux.cvno3OLD(KT) = NitrogenFlux.cvno3(KT-1);
end

phys = equations.physicalConstants();
R = phys.R;
Tr = phys.Tr;
PeCr = phys.PeCr;

% === 冻融 N 分配 (Killawee 1998 K_segr=0.3) ===
% 在反应+ADE 之前调一次, 让 N 在液相/冰相之间重分配, 然后 ADE 在分配后的
% C_NH4/C_NO3 上跑. 这样 ADE 看到的"上一步液相 N"已经扣除了被冰包合的部分,
% 不会双重浓缩.
% 输入: 上 KT 末冰量 (Theta_II_for_N_prev, 跟 Theta_LL_for_N_prev 同源时步)
%       本 KT 末冰量 (SoilVariables.Theta_II, 水量解之后)
% 输出: C_NH4/C_NO3 液相被扣 (或释放), C_NH4_ice/C_NO3_ice 冰相同步更新
%       partition_diag.delta_*_to_ice / from_ice [μg/cm²] 进守恒账本
if isfield(SoilVariables, 'Theta_II_for_N_prev') && ~isempty(SoilVariables.Theta_II_for_N_prev)
    Theta_II_old_for_partition = SoilVariables.Theta_II_for_N_prev;
else
    Theta_II_old_for_partition = SoilVariables.Theta_II;   % KT=1: Δ=0, partition 不动
end
[SoilVariables.C_NH4, SoilVariables.C_NO3, SoilVariables.C_NH4_ice, SoilVariables.C_NO3_ice, partition_diag] = ...
    soil.nitrogen.chemistry.freezeThawPartitioning( ...
        SoilVariables.C_NH4, SoilVariables.C_NO3, ...
        SoilVariables.C_NH4_ice, SoilVariables.C_NO3_ice, ...
        SoilVariables.Theta_II, Theta_II_old_for_partition, ...
        SoilVariables.Theta_LL, SoilVariables.Theta_L, ...
        NL, NN, cfg, DeltZ, NitrogenVariables.RHOKG);
NitrogenDiagnostics.partition = partition_diag;
% 同步 _new: ADE Picard 起点用 C_*_new = C_* (post-partition)
SoilVariables.C_NH4_new = SoilVariables.C_NH4;
SoilVariables.C_NO3_new = SoilVariables.C_NO3;

% 2026-05-14 partition 后重算 cvno3OLD/cvnh4OLD: 让账本起点反映 partition 后状态
% (partition 改了 SoilVariables.C_NH4 和 C_NH4_ice, 但 cvno3(KT-1) 是 partition 前状态;
%  ADE 接下来 march forward 用 partition 后 C, 所以账本对账起点必须是 partition 后 cv)
%
% [2026-05-15 N-closure fix] 冰相 θ 必须用 Theta_II_new (= SoilVariables.Theta_II) 而非
% _old_for_partition. 原因：partition.m line 162-163 反算 C_NO3_ice 时除以 theta_ice_new
% (C_ice_new = N_ice_new / θ_ice_new). 因此 calculateStorage(C_ice, θ_ice_old) 会把
% N 错乘 θ_old/θ_new 倍，冬季 θ_ice 涨期单步漏 ~67 μg/cm² NO3 → 72 KT 累积 -559.
% 用 Theta_II_new 跟 partition 内部反算一致, partition 出口数学守恒 (验证: cvno3OLD == cvno3(KT-1)).
% [2026-05-17 dissolution pool] 加 C_NH4_granule 到 NH4 储量 (partition 后再算)
gran_old = SoilVariables.C_NH4_granule;
if isfield(SoilVariables, 'C_NH4_granule_old')
    gran_old = SoilVariables.C_NH4_granule_old;
end
if isempty(gran_old), gran_old = zeros(NN, 1); end
[NitrogenFlux.cvnh4OLD(KT), NitrogenFlux.cvno3OLD(KT)] = ...
    soil.nitrogen.transport.calculateStorage( ...
        Theta_Solute_L, SoilVariables.C_NH4, SoilVariables.C_NO3, ...
        NitrogenVariables.RHOKG, DeltZ, NL, ...
        SoilVariables.Theta_II, SoilVariables.C_NH4_ice, SoilVariables.C_NO3_ice, ...
        gran_old);

% Apply the event at the start of the accepted N interval.  Both the old-side
% coefficient field and Picard initial candidate see the source, while
% cv*OLD above remains the genuine pre-event accepted inventory.
SoilVariables.C_NH4 = SoilVariables.C_NH4 + fert_delta_nh4;
SoilVariables.C_NO3 = SoilVariables.C_NO3 + fert_delta_no3;
SoilVariables.C_Urea = SoilVariables.C_Urea + fert_delta_urea;
SoilVariables.C_NH4_new = SoilVariables.C_NH4;
SoilVariables.C_NO3_new = SoilVariables.C_NO3;
SoilVariables.C_Urea_new = SoilVariables.C_Urea;

% The reservation and chemistry both belong to this one physical interval.
% Freeze partition and the accepted management event are already applied;
% Picard candidates must not become a new inventory or O2 time origin.
[inventoryNh4Top,inventoryNo3Top]=soil.organic_n.availableDinMgG( ...
    SoilVariables.C_NH4,SoilVariables.C_NO3,Theta_Solute_LL, ...
    NitrogenVariables.RHOKG,NitrogenVariables.RHO_bulk);
reactionInventoryBottom=[flip(inventoryNh4Top),flip(inventoryNo3Top)];
oxygenAtStepStart=NitrogenVariables.C_O2;
oxygenPreviousCandidate=oxygenAtStepStart;
aerobicPreviousCandidate=NitrogenVariables.f_aerobic;
oxygenConverged=false;
oxygenNitCap=inf(NL,1);

CHK2 = inf(NN, 1);
CHK3 = inf(NN, 1);
for iter = 1:max_iter
    [NitrogenVariables.Retard, NitrogenFlux.RetardBAR] = soil.nitrogen.chemistry.calculateRetardationFactors( ...
        Theta_Solute_LL, NitrogenVariables.RHOKG, NitrogenVariables.Retard, NitrogenFlux.RetardBAR, NL);
    [NitrogenVariables.Dsh, NitrogenVariables.Dif, NitrogenVariables.Disp, NitrogenVariables.Dsh_h, NitrogenVariables.Dsh_T, NitrogenVariables.Dsh_A] = ...
        soil.nitrogen.transport.soluteDispersion(NitrogenVariables.DIF1111, nS, NL, Theta_Solute_LL, SoilVariables.Theta_s(:).', ...
        AirVariables.QL, AirVariables.QL_h, AirVariables.QL_T, AirVariables.QL_A, NitrogenVariables.DISP, NN, DeltZ, ...
        R, Tr, SoilVariables.TT, NitrogenVariables.eT, Delt_t, SoilVariables.Srt, ...
        Theta_Solute_L, cfg.lArtD, cfg.lUpW, NitrogenVariables.RHOKG, ...
        PeCr, NitrogenVariables.Retard);
    if ExternalInputs.organic_enabled
        NitrogenVariables.Smin = repmat(ExternalInputs.organic_nh4_rate(:),1,2);
        NitrogenVariables.SorgNO3 = repmat(ExternalInputs.organic_no3_rate(:),1,2);
    else
        NitrogenVariables.Smin = soil.nitrogen.organic.mineralization(SoilVariables.TT, NitrogenVariables.K_mineralization, ...
            Theta_Solute_LL, NitrogenVariables.Theta_f, NL, Tr, R, ...
            NitrogenVariables.eT, NitrogenVariables.eTheta, NitrogenVariables.MSOC_volume, NitrogenVariables.Smin, PH1, NitrogenVariables.ePH);
        NitrogenVariables.SorgNO3 = zeros(NL,2);
    end
    NitrogenVariables.Svol = soil.nitrogen.chemistry.ammoniumVolatilization(SoilVariables.TT, SoilVariables.C_NH4, NL, PH1, KT, ...
        NitrogenVariables.Svol, NitrogenVariables.Tm, NitrogenVariables.u_2, NN, NitrogenBoundary.NoTime);

    if isfield(GasVariables, 'O2a_molcm3')
        C_O2_top_bc = pickGasDriver(GasVariables.O2a_molcm3, KT);
    else
        Oa_bc = pickOptionalGasDriver(GasVariables, 'Oa', KT, []);
        ea_bc = pickOptionalGasDriver(GasVariables, 'ea', KT, 0);
        C_O2_top_bc = soil.gas.calcAtmO2(SoilVariables.P_g(NN) / 10, SoilVariables.TT(NN), Oa_bc, ea_bc);
    end
    if ExternalInputs.organic_enabled
        organicCarbonRate = repmat(ExternalInputs.organic_co2_c_rate(:),1,2);
    else
        organicCarbonRate = NitrogenVariables.Smin;
    end
    R_O2_sink_bc = soil.gas.calcO2Sink(NitrogenVariables.Snit, organicCarbonRate, NL);
    NitrogenVariables.R_O2_sink = R_O2_sink_bc;   % [diag 2026-05-16] expose to main loop for reg.write
    oxygenPhase=SoilVariables.oxygen_phase;
    oxygenPhase.reserved_sink_mol_cm3_s=soil.gas.calcO2Sink( ...
        zeros(NL,2),organicCarbonRate,NL);
    oxygenSnitUsed=NitrogenVariables.Snit;
    [NitrogenVariables.C_O2,~,NitrogenVariables.f_aerobic, ...
        NitrogenVariables.f_anaerobic,oxygenLedger]=soil.gas.solveO2Transport( ...
        oxygenAtStepStart,SoilVariables.Theta_g,SoilVariables.POR, ...
        SoilVariables.TT,GasVariables.V_A,R_O2_sink_bc,DeltZ,Delt_t,NL,NN, ...
        C_O2_top_bc,[],oxygenPhase);
    limitedOxygen=oxygenLedger.additional_scale<1-1e-10;
    allowed=mean(oxygenSnitUsed,2).*oxygenLedger.additional_scale;
    oxygenNitCap(limitedOxygen)=min(oxygenNitCap(limitedOxygen),allowed(limitedOxygen));
    NitrogenDiagnostics.o2_fallback=false;
    oxygenChange=max(abs(NitrogenVariables.C_O2(:)-oxygenPreviousCandidate(:)));
    aerobicChange=max(abs(NitrogenVariables.f_aerobic(:)-aerobicPreviousCandidate(:)));
    oxygenConverged=iter>1 && (NitrogenDiagnostics.o2_fallback || ...
        (oxygenChange<=1e-3*max(1e-10,max(abs(oxygenAtStepStart(:)))) ...
        &&aerobicChange<=1e-3));
    oxygenPreviousCandidate=NitrogenVariables.C_O2;
    aerobicPreviousCandidate=NitrogenVariables.f_aerobic;
    % [2026-05-16] N 化学改 Michaelis-Menten (DAISY-style)，废弃线性 K_nitrification
    % Snit 显式进 G1 (NH4 汇) + G2 (NO3 源), SVnit 闲置 = 0
    % Vmax_nit / Km_nit 从 Excel 'Vmax_nit' / 'Km_nit' 列读 (见 SoilProperties.m)
    [NitrogenVariables.Snit, NitrogenVariables.SVnit, NitrogenFlux.Cvnit_N2o, ...
        NitrogenVariables.nit_mm_factor, NitrogenVariables.nit_M_total, ...
        NitrogenVariables.eT_nit, NitrogenVariables.Snit_N2o] = ...
        soil.nitrogen.chemistry.nitrification( ...
            SoilVariables.TT, SoilVariables.C_NH4, NL, Theta_Solute_LL, ...
            NitrogenVariables.Theta_f, SoilVariables.Theta_s, ...
            NitrogenVariables.RHOKG, NitrogenVariables.Vmax_nit, NitrogenVariables.Km_nit, ...
            NitrogenVariables.Snit, R, Tr, NitrogenVariables.ePH, NitrogenFlux.Cvnit_N2o, ...
            PH1, NitrogenVariables.SVnit, KT, DeltZ, NitrogenVariables.f_aerobic);
    [NitrogenVariables.Sden, NitrogenVariables.SVden, NitrogenFlux.Cvdenit_N2o, NitrogenFlux.Cvdenit_N2] = ...
        soil.nitrogen.chemistry.denitrification(SoilVariables.TT, NitrogenVariables.Vmax_denit, Theta_Solute_LL, ...
        SoilVariables.C_NO3_new, NL, NitrogenVariables.eT, NitrogenVariables.eTheta, NitrogenVariables.Sden, ...
        R, Tr, NitrogenFlux.Cvdenit_N2o, NitrogenFlux.Cvdenit_N2, NitrogenVariables.SVden, ...
        KT, DeltZ, NitrogenVariables.Km_denit, NitrogenVariables.f_anaerobic);   % [2026-05-17 A2 土基对齐]
    NitrogenDiagnostics.reaction_inventory = reactionInventoryLedger( ...
        NitrogenVariables.Snit,NitrogenVariables.Sden,NitrogenVariables.Svol, ...
        NitrogenVariables.SinkS,NitrogenVariables.Smin,NitrogenVariables.SorgNO3,Delt_t);
    reactionReservationTop=zeros(NL,2);
    if isfield(ExternalInputs,'pre_reaction_reservation_mg_g') ...
            &&~isempty(ExternalInputs.pre_reaction_reservation_mg_g)
        reactionReservationTop=ExternalInputs.pre_reaction_reservation_mg_g;
    elseif isfield(ExternalInputs,'root_reservation_mg_g') ...
            &&~isempty(ExternalInputs.root_reservation_mg_g)
        reactionReservationTop=ExternalInputs.root_reservation_mg_g;
    end
    reactionReservationBottom=flip(reactionReservationTop);
    [NitrogenVariables.Snit,NitrogenVariables.Sden,NitrogenVariables.SinkS, ...
        NitrogenDiagnostics.reaction_allocation]=soil.nitrogen.limitReactionRates( ...
        NitrogenVariables.Snit,NitrogenVariables.Sden,NitrogenVariables.SinkS, ...
        reactionInventoryBottom,NitrogenVariables.RHO_bulk,Delt_t, ...
        reactionReservationBottom);
    oxygenFactor=min(1,oxygenNitCap./max(mean(NitrogenVariables.Snit,2),realmin));
    oxygenNUnmet=mean(NitrogenVariables.Snit,2).*(1-oxygenFactor) ...
        .*Delt_t./NitrogenVariables.RHO_bulk(:)*1e-3;
    NitrogenVariables.Snit=NitrogenVariables.Snit.*oxygenFactor;
    NitrogenDiagnostics.reaction_allocation.actual_mg_g(:,1)= ...
        NitrogenDiagnostics.reaction_allocation.actual_mg_g(:,1)-oxygenNUnmet;
    NitrogenDiagnostics.reaction_allocation.unmet_mg_g(:,1)= ...
        NitrogenDiagnostics.reaction_allocation.unmet_mg_g(:,1)+oxygenNUnmet;
    NitrogenDiagnostics.reaction_allocation.remaining_mg_g(:,1)= ...
        NitrogenDiagnostics.reaction_allocation.remaining_mg_g(:,1)+oxygenNUnmet;
    NitrogenDiagnostics.reaction_allocation.scale(:,1)= ...
        NitrogenDiagnostics.reaction_allocation.scale(:,1).*oxygenFactor;
    oxygenConverged=oxygenConverged&&all(abs(NitrogenVariables.Snit(:)-oxygenSnitUsed(:)) ...
        <=1e-5*max(1e-12,max(abs(NitrogenVariables.Snit(:)))));
    NitrogenVariables.Snit_N2o=0.02*NitrogenVariables.Snit;
    NitrogenFlux.Cvnit_N2o(KT)=DeltZ(:)'*mean(NitrogenVariables.Snit_N2o,2);
    beforeNodalDenit=NitrogenVariables.Sden;
    [NitrogenVariables.Sden,nodalCap]=soil.nitrogen.limitDenitrificationAtNodes( ...
        beforeNodalDenit,.98*NitrogenVariables.Snit-NitrogenVariables.SinkS(:,:,2)+NitrogenVariables.SorgNO3, ...
        SoilVariables.C_NO3,Theta_Solute_L,DeltZ,Delt_t);
    nodalUnmet=mean(beforeNodalDenit-NitrogenVariables.Sden,2).*Delt_t./NitrogenVariables.RHO_bulk(:)*1e-3;
    NitrogenDiagnostics.reaction_allocation.actual_mg_g(:,2)=NitrogenDiagnostics.reaction_allocation.actual_mg_g(:,2)-nodalUnmet;
    NitrogenDiagnostics.reaction_allocation.unmet_mg_g(:,2)=NitrogenDiagnostics.reaction_allocation.unmet_mg_g(:,2)+nodalUnmet;
    NitrogenDiagnostics.reaction_allocation.remaining_mg_g(:,2)=NitrogenDiagnostics.reaction_allocation.remaining_mg_g(:,2)+nodalUnmet;
    NitrogenDiagnostics.reaction_allocation.scale(:,2)=NitrogenDiagnostics.reaction_allocation.scale(:,2).*nodalCap.element_scale;
    NitrogenDiagnostics.nodal_denitrification=nodalCap;
    NitrogenFlux.Cvdenit_N2o(KT)=0.11*DeltZ(:)'*mean(NitrogenVariables.Sden,2);
    NitrogenFlux.Cvdenit_N2(KT)=0.89*DeltZ(:)'*mean(NitrogenVariables.Sden,2);
    NitrogenDiagnostics.reaction_inventory_basis='post_partition_post_event_step_start_bottom_up';
    [NitrogenVariables.Shys, NitrogenVariables.Vpore, NitrogenVariables.WFPS, DeltaUrea] = soil.nitrogen.chemistry.ureaHydrolysis( ...
        NitrogenVariables.Kurea, NL, SoilVariables.C_Urea, NitrogenVariables.Vpore, NitrogenVariables.WFPS, ...
        SoilVariables.Theta_g, Theta_Solute_LL, NitrogenVariables.Shys, Delt_t, SoilVariables.TT, DeltZ);

    % [2026-05-17 dissolution pool] NH4 颗粒肥溶解 (mirrors urea hydrolysis structure).
    %   Picard iter 内每次重算 Sdis 但**用同一个 granule 池**, ADE G1 += Sdis 累积.
    %   granule 池在 iter 收敛后才扣 (line 270, 与 C_Urea_new 同时更新), 这与
    %   ureaHydrolysis 完全对称: ureaHydrolysis 也是 iter 内反复算 Shys 而 Urea 池
    %   在 iter 末才扣 DeltaUrea.
    %
    %   守恒性: 假设 iter=K 次收敛, 每次 Sdis 都用同一个 granule_old, dG 不变,
    %   ADE 每次 G1 加 Sdis 注入 dissolved. iter 收敛后 granule -= dG, 此时
    %   dissolved 已经累计 K 次 Sdis 注入. **但只该扣 1 次 granule!**
    %   实际上, dissolveSdis(granule_old, dt) 已经算了"过 Δt 该释放多少", iter 内
    %   重复算 Sdis 是数值上 Picard 迭代收敛 dissolved 用的常量, 不重复物理释放.
    %   关键: 实际 dissolved 增量 = G1 净源 × Δt (ADE 求解一次, 不是 K 次), 守恒 OK.
    if isfield(SoilVariables, 'C_NH4_granule')
        [NitrogenVariables.Sdis, DeltaGranule] = soil.nitrogen.fertilizer.nh4FertDissolution( ...
            NitrogenVariables.K_nh4_dissolve, NL, SoilVariables.C_NH4_granule, ...
            NitrogenVariables.WFPS, Delt_t, SoilVariables.TT);
    else
        NitrogenVariables.Sdis = zeros(NL, 2);
        DeltaGranule = zeros(NN, 1);
    end

    [NitrogenVariables, SoilVariables, NitrogenFlux] = soil.nitrogen.transport.solveTransportStep( ...
        NitrogenVariables, SoilVariables, NitrogenFlux, NitrogenBoundary, AirVariables, ...
        Theta_Solute_L, Theta_Solute_LL, cfg, PeCr, KT, Delt_t);

    SoilVariables.C_Urea_new = SoilVariables.C_Urea - DeltaUrea;
    % Retain small positive stocks: threshold deletion would lose nitrogen.

    error_toleranceNno3 = solver.calculateErrorTolerance(SoilVariables.C_NO3_new(1:NN), error_thresholds, error_tolerances);
    error_toleranceNn = solver.calculateErrorTolerance(SoilVariables.C_NH4_new(1:NN), error_thresholds, error_tolerances);
    CHK2 = NitrogenVariables.CHK2;
    CHK3 = NitrogenVariables.CHK3;
    if all(CHK3(:) < tolerance_factor * error_toleranceNno3(:)) && ...
            all(CHK2(:) < tolerance_factor * error_toleranceNn(:)) ...
            &&oxygenConverged
        break
    end
end
assert(iter<max_iter || (all(CHK3(:) < tolerance_factor * error_toleranceNno3(:)) ...
    &&all(CHK2(:) < tolerance_factor * error_toleranceNn(:))&&oxygenConverged), ...
    'soil:nitrogen:picardNonconvergence', ...
    'Nitrogen/O2 Picard iteration did not converge for one physical interval.');

% [2026-05-17 dissolution pool] Update granule pool ONCE per outer step (after Picard converged).
finalDemand=soil.gas.calcO2Sink(NitrogenVariables.Snit,organicCarbonRate,NL);
[NitrogenVariables.C_O2,~,NitrogenVariables.f_aerobic, ...
    NitrogenVariables.f_anaerobic,NitrogenDiagnostics.oxygen]=soil.gas.solveO2Transport( ...
    oxygenAtStepStart,SoilVariables.Theta_g,SoilVariables.POR, ...
    SoilVariables.TT,GasVariables.V_A,finalDemand,DeltZ,Delt_t,NL,NN, ...
    C_O2_top_bc,[],oxygenPhase);
assert(all(NitrogenDiagnostics.oxygen.additional_scale>=1-1e-9), ...
    'soil:gas:oxygenReactionMismatch','Accepted nitrogen demand exceeds actual O2 consumption.');
NitrogenVariables.R_O2_sink=NitrogenDiagnostics.oxygen.actual_sink_mol_cm3_s;
NitrogenDiagnostics.oxygen_nitrification_cap=oxygenNitCap;
%   Why outside iter loop: Picard iter 不该改 prognostic state, 每次 iter 算 Sdis 用同一 granule_old.
%   ADE 内部稳态求解 dissolved 一次, granule 也只该扣 1 次 ΔG (= ∫Sdis dt). 这样 mass conservation:
%       Δ(dissolved + granule) = External (Fert + Smin + Bot - Snit - Svol - SinkS)
%   (granule → dissolved 是内部转移, Sdis 项在 ADE G1 已加, granule 扣 DeltaGranule, 自动平衡)
if isfield(SoilVariables, 'C_NH4_granule')
    SoilVariables.C_NH4_granule = SoilVariables.C_NH4_granule - DeltaGranule;
    SoilVariables.C_NH4_granule(SoilVariables.C_NH4_granule <= 1e-6) = 0;
end

% 储量包含冰相 (P0a.5 Killawee 1998 冻融 N 分配)
% [2026-05-17 dissolution pool] 加 NH4 granule 池 (内部转移到 dissolved 后剩余)
gran_end = SoilVariables.C_NH4_granule;
if isempty(gran_end), gran_end = zeros(NN, 1); end
[NitrogenFlux.cvnh4(KT), NitrogenFlux.cvno3(KT)] = ...
    soil.nitrogen.transport.calculateStorage( ...
        Theta_Solute_LL, SoilVariables.C_NH4_new, SoilVariables.C_NO3_new, ...
        NitrogenVariables.RHOKG, DeltZ, NL, ...
        SoilVariables.Theta_II, SoilVariables.C_NH4_ice, SoilVariables.C_NO3_ice, ...
        gran_end);

% 保存本 KT 末 θ_LL / θ_II 快照, 下 KT 进 N 模块时作 ADE 真"上 KT 末" θ.
% [2026-05-15 N-closure fix] 保存 effective Theta_Solute_LL 而不是 raw Theta_LL —
% cvno3(KT) 在 line 247-249 用 effective (ponding-aware + POR-cap), 下 KT 算
% cvno3OLD 必须用同一个 effective θ, 否则 ponding 节点 ΔStorage 跟 ledger 对不上,
% 跨 KT 累积 ~7000 μg/cm² 残差 (冬季冰相剧烈期暴露). raw Theta_LL 在 ponding 节点
% 比 effective 小 hh/DeltZ 倍 → 单步漏 ~67 μg/cm² NO3.
SoilVariables.Theta_LL_for_N_prev = Theta_Solute_LL;
SoilVariables.Theta_II_for_N_prev = SoilVariables.Theta_II;

NitrogenDiagnostics.iterations = iter;
NitrogenDiagnostics.Theta_Solute_LL = Theta_Solute_LL;
NitrogenDiagnostics.Theta_Solute_L = Theta_Solute_L;
cvOrganicNO3 = DeltZ(:)' * ExternalInputs.organic_no3_rate(:);
NitrogenDiagnostics.organic_enabled = ExternalInputs.organic_enabled;
NitrogenDiagnostics.organic_nh4_rate = ExternalInputs.organic_nh4_rate;
NitrogenDiagnostics.organic_no3_rate = ExternalInputs.organic_no3_rate;
NitrogenDiagnostics.organic_co2_c_rate = ExternalInputs.organic_co2_c_rate;
NitrogenDiagnostics.cvOrganicNO3 = cvOrganicNO3;
NitrogenDiagnostics.CHK2 = CHK2;
NitrogenDiagnostics.CHK3 = CHK3;
% [2026-05-16] cvSnit 已是 NH4 100% 总消耗 (含 NO3+N2O 两条路径), Cvnit_N2o 是 cvSnit 子集, 不再独立计入 NH4 ledger
NitrogenDiagnostics.Massbalance_absolute_error_nh4 = ...
    (NitrogenFlux.cvnh4(KT) - NitrogenFlux.cvnh4OLD(KT)) / Delt_t - ...
    (NitrogenFlux.cvTop(KT) + NitrogenFlux.cvBot(KT) - NitrogenFlux.cvSnit(KT) + ...
    NitrogenFlux.cvSmin(KT) - NitrogenFlux.cvSinkSnh4(KT) - ...
    NitrogenFlux.cvSvol(KT) + ExternalInputs.nh4_added / Delt_t);
denom_nh4 = max(abs((NitrogenFlux.cvnh4(KT) - NitrogenFlux.cvnh4OLD(KT)) / Delt_t), ...
    abs(NitrogenFlux.cvTop(KT)) + abs(NitrogenFlux.cvBot(KT)) + abs(NitrogenFlux.cvSnit(KT)) + ...
    abs(NitrogenFlux.cvSmin(KT)) + abs(NitrogenFlux.cvSinkSnh4(KT)) + ...
    abs(NitrogenFlux.cvSvol(KT)) + ExternalInputs.nh4_added / Delt_t);
if denom_nh4 == 0, denom_nh4 = 1e-10; end
NitrogenDiagnostics.Massbalance_relative_error_nh4 = 100 * abs(NitrogenDiagnostics.Massbalance_absolute_error_nh4) / denom_nh4;

% [2026-05-16] NO3 in 是 0.98×cvSnit (硝化产物 98% 入 NO3, 2% → N2O 不入), 跟 G2 = 0.98×Snit ... 一致
NitrogenDiagnostics.Massbalance_absolute_error_no3 = ...
    (NitrogenFlux.cvno3(KT) - NitrogenFlux.cvno3OLD(KT)) / Delt_t - ...
    (NitrogenFlux.cvTop1(KT) + NitrogenFlux.cvBot1(KT) - NitrogenFlux.cvSinkSno3(KT) - ...
    NitrogenFlux.cvSden(KT) + 0.98 * NitrogenFlux.cvSnit(KT) + cvOrganicNO3 + ...
    ExternalInputs.no3_added / Delt_t);
denom_no3 = max(abs((NitrogenFlux.cvno3(KT) - NitrogenFlux.cvno3OLD(KT)) / Delt_t), ...
    abs(NitrogenFlux.cvTop1(KT)) + abs(NitrogenFlux.cvBot1(KT)) + abs(NitrogenFlux.cvSinkSno3(KT)) + ...
    abs(NitrogenFlux.cvSden(KT)) + abs(0.98 * NitrogenFlux.cvSnit(KT)) + abs(cvOrganicNO3) + ...
    ExternalInputs.no3_added / Delt_t);
if denom_no3 == 0, denom_no3 = 1e-10; end
NitrogenDiagnostics.Massbalance_relative_error_no3 = 100 * abs(NitrogenDiagnostics.Massbalance_absolute_error_no3) / denom_no3;

% [ClimaLand-style cumulative balance, see SKILL.md §6]
% Per-step accumulators for whole-run N conservation check.
% Δstorage − Σ(inputs − outputs) over [1, KT] must remain near zero.
if ~isfield(NitrogenFlux, 'cum_nh4_in_mg_cm2')
    NitrogenFlux.cum_nh4_in_mg_cm2  = 0;
    NitrogenFlux.cum_nh4_out_mg_cm2 = 0;
    NitrogenFlux.cum_no3_in_mg_cm2  = 0;
    NitrogenFlux.cum_no3_out_mg_cm2 = 0;
    NitrogenFlux.cv_nh4_init_mg_cm2 = NitrogenFlux.cvnh4OLD(KT);
    NitrogenFlux.cv_no3_init_mg_cm2 = NitrogenFlux.cvno3OLD(KT);
    % Per-process cumulative tracking for debug (μg/cm²)
    NitrogenFlux.cum_cvTop_NH4    = 0;
    NitrogenFlux.cum_cvBot_NH4    = 0;
    NitrogenFlux.cum_cvSmin       = 0;
    NitrogenFlux.cum_cvSnit       = 0;
    NitrogenFlux.cum_cvSinkSnh4   = 0;
    NitrogenFlux.cum_cvSvol       = 0;
    NitrogenFlux.cum_Cvnit_N2o    = 0;
    NitrogenFlux.cum_cvTop_NO3    = 0;
    NitrogenFlux.cum_cvBot_NO3    = 0;
    NitrogenFlux.cum_cvSinkSno3   = 0;
    NitrogenFlux.cum_cvSden       = 0;
    NitrogenFlux.cum_iceNH4_in    = 0;
    NitrogenFlux.cum_iceNH4_out   = 0;
    NitrogenFlux.cum_iceNO3_in    = 0;
    NitrogenFlux.cum_iceNO3_out   = 0;
    NitrogenFlux.cum_cvShys       = 0;   % urea → NH4 hydrolysis (in for NH4)
end
% Depth-integrate urea hydrolysis source: Shys is [μg/(cm³·s), NL×nD]
% Same trapezoidal rule the transport assemble uses (Smin / Snit ...) for unit consistency.
Shys_avg_now = (NitrogenVariables.Shys(:,1) + NitrogenVariables.Shys(:,2)) * 0.5;
cvShys_now   = DeltZ(:)' * Shys_avg_now;   % [μg/(cm²·s)]
NitrogenFlux.cum_cvShys = NitrogenFlux.cum_cvShys + Delt_t * cvShys_now;
% Per-process cumulative integral (μg/cm² over [1..KT])
NitrogenFlux.cum_cvTop_NH4  = NitrogenFlux.cum_cvTop_NH4  + Delt_t * NitrogenFlux.cvTop(KT);
NitrogenFlux.cum_cvBot_NH4  = NitrogenFlux.cum_cvBot_NH4  + Delt_t * NitrogenFlux.cvBot(KT);
NitrogenFlux.cum_cvSmin     = NitrogenFlux.cum_cvSmin     + Delt_t * NitrogenFlux.cvSmin(KT);
NitrogenFlux.cum_cvSnit     = NitrogenFlux.cum_cvSnit     + Delt_t * NitrogenFlux.cvSnit(KT);
NitrogenFlux.cum_cvSinkSnh4 = NitrogenFlux.cum_cvSinkSnh4 + Delt_t * NitrogenFlux.cvSinkSnh4(KT);
NitrogenFlux.cum_cvSvol     = NitrogenFlux.cum_cvSvol     + Delt_t * NitrogenFlux.cvSvol(KT);
NitrogenFlux.cum_Cvnit_N2o  = NitrogenFlux.cum_Cvnit_N2o  + Delt_t * NitrogenFlux.Cvnit_N2o(KT);
NitrogenFlux.cum_cvTop_NO3  = NitrogenFlux.cum_cvTop_NO3  + Delt_t * NitrogenFlux.cvTop1(KT);
NitrogenFlux.cum_cvBot_NO3  = NitrogenFlux.cum_cvBot_NO3  + Delt_t * NitrogenFlux.cvBot1(KT);
NitrogenFlux.cum_cvSinkSno3 = NitrogenFlux.cum_cvSinkSno3 + Delt_t * NitrogenFlux.cvSinkSno3(KT);
NitrogenFlux.cum_cvSden     = NitrogenFlux.cum_cvSden     + Delt_t * NitrogenFlux.cvSden(KT);
% [2026-05-16] NH4 ledger: +Top +Bot +Smin +Shys, −Snit −SinkSnh4 −Svol
% (Snit 已含 NO3+N2O 全部 NH4 消耗, Cvnit_N2o 是 Snit 子集不再独立)
NitrogenFlux.cum_nh4_in_mg_cm2 = NitrogenFlux.cum_nh4_in_mg_cm2 + Delt_t * ( ...
    NitrogenFlux.cvTop(KT) + NitrogenFlux.cvBot(KT) + NitrogenFlux.cvSmin(KT) + cvShys_now);
NitrogenFlux.cum_nh4_out_mg_cm2 = NitrogenFlux.cum_nh4_out_mg_cm2 + Delt_t * ( ...
    NitrogenFlux.cvSnit(KT) + NitrogenFlux.cvSinkSnh4(KT) + NitrogenFlux.cvSvol(KT));
% [2026-05-16] NO3 ledger: +Top +Bot +0.98×Snit (硝化 98% 入 NO3, 2%→N2O 不入), −SinkSno3 −Sden
NitrogenFlux.cum_no3_in_mg_cm2 = NitrogenFlux.cum_no3_in_mg_cm2 + Delt_t * ( ...
    NitrogenFlux.cvTop1(KT) + NitrogenFlux.cvBot1(KT) + 0.98 * NitrogenFlux.cvSnit(KT));
NitrogenFlux.cum_no3_out_mg_cm2 = NitrogenFlux.cum_no3_out_mg_cm2 + Delt_t * ( ...
    NitrogenFlux.cvSinkSno3(KT) + NitrogenFlux.cvSden(KT));
NitrogenFlux.cum_no3_in_mg_cm2 = NitrogenFlux.cum_no3_in_mg_cm2 + ...
    Delt_t * max(cvOrganicNO3,0);
NitrogenFlux.cum_no3_out_mg_cm2 = NitrogenFlux.cum_no3_out_mg_cm2 + ...
    Delt_t * max(-cvOrganicNO3,0);
NitrogenFlux.cum_nh4_in_mg_cm2 = NitrogenFlux.cum_nh4_in_mg_cm2 + ...
    ExternalInputs.nh4_added;
NitrogenFlux.cum_no3_in_mg_cm2 = NitrogenFlux.cum_no3_in_mg_cm2 + ...
    ExternalInputs.no3_added;
if ~isfield(NitrogenFlux, 'cum_Fert_NH4'), NitrogenFlux.cum_Fert_NH4 = 0; end
if ~isfield(NitrogenFlux, 'cum_Fert_NO3'), NitrogenFlux.cum_Fert_NO3 = 0; end
if ~isfield(NitrogenFlux, 'cum_Fert_Urea'), NitrogenFlux.cum_Fert_Urea = 0; end
NitrogenFlux.cum_Fert_NH4 = NitrogenFlux.cum_Fert_NH4 + ExternalInputs.nh4_added;
NitrogenFlux.cum_Fert_NO3 = NitrogenFlux.cum_Fert_NO3 + ExternalInputs.no3_added;
NitrogenFlux.cum_Fert_Urea = NitrogenFlux.cum_Fert_Urea + ExternalInputs.urea_added;

% --- P0a.5 freeze-thaw partitioning: ice-phase N diagnostic tracking only ---
% 不加进 cum_nh4_in/out (cvnh4 = 液相+冰相总和, partition 在总账本中是内部转移, 中性)
% 只保留 cum_iceN*_in/out 作为单冰相诊断 (输出报告参考)
if isfield(NitrogenDiagnostics, 'partition') && ~isempty(NitrogenDiagnostics.partition)
    pd = NitrogenDiagnostics.partition;
    NitrogenFlux.cum_iceNH4_in  = NitrogenFlux.cum_iceNH4_in  + pd.delta_NH4_from_ice;
    NitrogenFlux.cum_iceNH4_out = NitrogenFlux.cum_iceNH4_out + pd.delta_NH4_to_ice;
    NitrogenFlux.cum_iceNO3_in  = NitrogenFlux.cum_iceNO3_in  + pd.delta_NO3_from_ice;
    NitrogenFlux.cum_iceNO3_out = NitrogenFlux.cum_iceNO3_out + pd.delta_NO3_to_ice;
end

NitrogenDiagnostics.cum_residual_nh4 = ...
    (NitrogenFlux.cvnh4(KT) - NitrogenFlux.cv_nh4_init_mg_cm2) ...
    - (NitrogenFlux.cum_nh4_in_mg_cm2 - NitrogenFlux.cum_nh4_out_mg_cm2);
NitrogenDiagnostics.cum_residual_no3 = ...
    (NitrogenFlux.cvno3(KT) - NitrogenFlux.cv_no3_init_mg_cm2) ...
    - (NitrogenFlux.cum_no3_in_mg_cm2 - NitrogenFlux.cum_no3_out_mg_cm2);
end

function inputs = normalizeExternalInputs(candidate)
% Keep the solver adapter thin; independently tested organic/event gate.
% Editor: Jizixing Zhu; Date: 2026-09-11.
    inputs=soil.nitrogen.normalizeOrganicExternalInputs(candidate);
end

function ledger = reactionInventoryLedger(snit,sden,svol,sink,smin,sorgNo3,dt)
%REACTIONINVENTORYLEDGER Record native-rate request/actual/unmet fields.
% Products are explicitly kept for a later named stage.
% Editor: Jizixing Zhu; Date: 2026-09-12.
    z=size(snit); zero=zeros(z);
    ledger=struct('dt_s',dt,'product_policy', ...
        'products become available only in their next named stage');
    ledger.nh4=struct('nitrification_requested',snit,'nitrification_actual',snit, ...
        'nitrification_unmet',zero,'volatilization_requested',svol, ...
        'volatilization_actual',svol,'volatilization_unmet',zero, ...
        'root_requested',sink(:,:,1),'root_actual',sink(:,:,1),'root_unmet',zero, ...
        'organic_source',max(smin,0),'organic_sink_requested',max(-smin,0), ...
        'organic_sink_actual',max(-smin,0),'organic_sink_unmet',zero);
    ledger.no3=struct('denitrification_requested',sden,'denitrification_actual',sden, ...
        'denitrification_unmet',zero,'root_requested',sink(:,:,2), ...
        'root_actual',sink(:,:,2),'root_unmet',zero,'organic_source',max(sorgNo3,0), ...
        'organic_sink_requested',max(-sorgNo3,0),'organic_sink_actual',max(-sorgNo3,0), ...
        'organic_sink_unmet',zero);
end

function v = pickOptionalGasDriver(GasVariables, fieldName, KT, defaultValue)
%PICKOPTIONALGASDRIVER  Read scalar/time-series gas driver with a fallback.
    if isfield(GasVariables, fieldName) && ~isempty(GasVariables.(fieldName))
        v = pickGasDriver(GasVariables.(fieldName), KT);
    else
        v = defaultValue;
    end
end

function v = pickGasDriver(x, KT)
%PICKGASDRIVER  Return a scalar gas driver at KT from scalar or time series.
    if isscalar(x)
        v = x;
    else
        idx = min(KT, numel(x));
        v = x(idx);
    end
end
