function NitrogenVariables = setInitialNitrogen(SoilProperties, VanGenuchten, MeteoData, cfg, state)
% setInitialNitrogen - 初始化氮素模块的静态参数与动态求解状态
%                      Initialize static N parameters and dynamic solver state
%
% 输入 / Inputs:
%   SoilProperties - legacy SoilProperties 对象
%   VanGenuchten   - VG 参数结构体（含 Theta_f）
%   MeteoData      - 气象数据对象（含 Ta/U 时间序列）
%   cfg            - legacy ModelSettings 对象
%   state          - StateVariables 对象（含初始动态状态字段）
%
% 输出 / Outputs:
%   NitrogenVariables - 氮素平坦结构体，含静态参数 + 动态求解状态
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-27

    if cfg.lNefc ~= 1
        NitrogenVariables = struct();
        return;
    end

    NL = cfg.Grid.NL;
    nD = cfg.Grid.nD;
    nS = cfg.Grid.nS;
    NN = cfg.Grid.NN;
    tS = cfg.Time.tS;

    % --- 静态运移系数 / Static transport coefficients ---
    NitrogenVariables.RHO_bulk = SoilProperties.Bulk_density(:);
    NitrogenVariables.RHOKG    = SoilProperties.KG(:) .* NitrogenVariables.RHO_bulk;

    NitrogenVariables.DIF1111          = zeros(1, NL, nS);
    NitrogenVariables.DIF1111(1, :, 1) = SoilProperties.DIF_NH4(:).';
    NitrogenVariables.DIF1111(1, :, 2) = SoilProperties.DIF_NO3(:).';
    NitrogenVariables.DISP             = SoilProperties.DISP(:);

    % --- Static biochemical rates ---
    NitrogenVariables.K_mineralization  = SoilProperties.K_mineralization(:);
    NitrogenVariables.Vmax_nit          = SoilProperties.Vmax_nit(:);            % [2026-05-16 MM]  硝化 V_max [1/s]
    NitrogenVariables.Km_nit            = SoilProperties.Km_nit(:);              % [2026-05-16 MM]  硝化 K_m [μg/cm³ soil]
    NitrogenVariables.Vmax_denit = SoilProperties.Vmax_denit(:);   % [2026-05-17 A2] 反硝化 V_max [μg/(cm³ 土·s)] (Excel 'Vmax_denit' /86400)
    NitrogenVariables.Km_denit   = SoilProperties.Km_denit(:);     % [2026-05-17 A2] 反硝化 K_m  [μg N/cm³ 土] (Excel 'Km_denit', 土基)
    NitrogenVariables.K_nh4_dissolve = SoilProperties.K_nh4_dissolve(:);   % [2026-05-17 dissolution pool] NH4 颗粒肥溶解 [1/s]
    NitrogenVariables.Kurea = SoilProperties.Kurea(:);
    NitrogenVariables.pH    = SoilProperties.pH(:);

    % SOM volumetric fraction [-], substrate basis for mineralization
    NitrogenVariables.MSOC_volume = SoilProperties.MSOC_volume(:);

    % --- VG 导出 / From VanGenuchten ---
    NitrogenVariables.Theta_f = VanGenuchten.Theta_f(:);

    % --- 气象导出 / From MeteoData ---
    if isprop(MeteoData, 'Ta_msr') && ~isempty(MeteoData.Ta_msr)
        NitrogenVariables.Tm = MeteoData.Ta_msr(:);
    else
        NitrogenVariables.Tm = MeteoData.Ta(:);
    end
    if isprop(MeteoData, 'U_msr') && ~isempty(MeteoData.U_msr)
        NitrogenVariables.u_2 = MeteoData.U_msr(:);
    else
        NitrogenVariables.u_2 = MeteoData.U(:);
    end

    % --- 动态求解状态 / Dynamic solver state ---
    NitrogenVariables.Retard = state.Retard;
    NitrogenVariables.Dsh    = state.DshNh4;
    NitrogenVariables.Dif    = state.DifNh4;
    NitrogenVariables.Disp   = state.Disp;

    NitrogenVariables.Dsh_h = zeros(NL, nD, nS);
    NitrogenVariables.Dsh_T = zeros(NL, nD, nS);
    NitrogenVariables.Dsh_A = zeros(NL, nD, nS);

    NitrogenVariables.Smin  = state.Smin;
    NitrogenVariables.Sden  = state.Sden;
    NitrogenVariables.SVden = state.SVden;
    NitrogenVariables.Snit  = state.Snit;
    NitrogenVariables.SVnit = state.SVnit;
    NitrogenVariables.Svol  = state.Svol;
    NitrogenVariables.SVvol = state.SVvol;
    NitrogenVariables.Shys  = state.Shys;
    NitrogenVariables.Vpore = state.Vpore;
    NitrogenVariables.WFPS  = state.WFPS;
    NitrogenVariables.eT    = state.eT;
    NitrogenVariables.eTheta = state.eTheta;
    NitrogenVariables.ePH   = state.ePH;

    NitrogenVariables.SinkS         = zeros(NL, nD, nS);
    NitrogenVariables.SinkS(:, :, 1) = state.SinkNh4;
    NitrogenVariables.SinkS(:, :, 2) = state.SinkNo3;

    % [2026-05-17 dissolution pool] Sdis lives in NitrogenVariables (rate, NL×nD),
    % whereas C_NH4_granule (NN×1 prognostic state) is owned by SoilVariables
    % (mirrors C_Urea convention; see StateVariables.m line 254-255).
    NitrogenVariables.Sdis = zeros(NL, nD);            % [μg/(cm³·s), NL×nD]  dissolution rate (granule→NH4)

    NitrogenVariables.C_O2       = [];      % 由 StartInit 填充初值
    NitrogenVariables.f_aerobic  = ones(NL, 2);
    NitrogenVariables.f_anaerobic = zeros(NL, 2);

    % --- 累计输出与收敛检查 / Accumulation and convergence arrays ---
    NitrogenVariables.N_SAVE    = zeros(3, 3, 5);
    NitrogenVariables.NH4_total = zeros(NN, tS);
    NitrogenVariables.NO3_total = zeros(NN, tS);
    NitrogenVariables.CHK2      = inf(NN, 1);
    NitrogenVariables.CHK3      = inf(NN, 1);

    % --- Phase C M2-A (2026-05-18): plant N labile pool ---
    % Single-pool plant N (DAISY-style luxury threshold); state at end of KT.
    % cn_coupling §12.10. C-primary: N_labile(1) = C_labile_0 / CN_L, assigned in the
    % C block below once CN_L is available (here just pre-allocate the array).
    NitrogenVariables.N_labile           = zeros(tS, 1);    % [g N/m²]  prognostic, end-of-KT

    % --- Phase C M2-B.1b/spin-up + C-primary (v12 2026-05-29): 4 organ pools, C 主 N 派生 ---
    % [v12 C 主导] 4 结构器官 (leaf/root/stem/fruit) 初值以 C 为主:
    %   config 给 C_*_init [g C/m²] (实测生物量, 或 N×CN 反算的等价值), N_X_0 = C_X_0 / CN_target_X 派生。
    %   理由: C 主导模型里 C_leaf→LAI→GPP 是驱动链, 被信任量应是 C 而非派生 C; 生物量比组织 N 好测;
    %         对标 FruitTree (Dabbaghi 2022 GMD) 以 gC (transplant) 初始化。
    % 例外: N_labile (NSC 供给池, 非器官) 仍以 N 为主 (上方读入), C_labile 派生 = N_labile × CN_L (见下 C 块);
    %       N_storage / C_storage (跨年储备) 无单一 C:N, 两者由 config 各自独立给。
    % CN_target 来自 PFT 表 (cfg.Veg.PFT_type); spin-up 末值写回 config 时直接更新 C_*_init。
    if cfg.lNefc == 1 && isprop(cfg.Veg, 'PFT_type')
        CN_L_      = vegetation.getPFTValue(cfg.Veg.PFT_type, 'CN_L');
        CN_root_   = vegetation.getPFTValue(cfg.Veg.PFT_type, 'CN_target_root');
        CN_fruit_  = vegetation.getPFTValue(cfg.Veg.PFT_type, 'CN_target_fruit');
        % [13c dual stem pool] live (CN_lw, respires) + dead (CN_dw, ~N-free) C:N
        CN_lw_     = vegetation.getPFTValue(cfg.Veg.PFT_type, 'CN_lw');
        CN_dw_     = vegetation.getPFTValue(cfg.Veg.PFT_type, 'CN_dw');
        flivewd_   = vegetation.getPFTValue(cfg.Veg.PFT_type, 'flivewd');
        % organs: C-primary (config C_*_init), N derived = C / CN_target
        C_leaf_0       = max(getInitVal(cfg, 'C_leaf_init',      0), 0);
        C_root_fine_0  = max(getInitVal(cfg, 'C_root_fine_init', 0), 0);
        C_fruit_0      = max(getInitVal(cfg, 'C_fruit_init',     0), 0);
        N_leaf_0       = C_leaf_0      / max(CN_L_,        1);
        N_root_fine_0  = C_root_fine_0 / max(CN_root_,     1);
        N_fruit_0      = C_fruit_0     / max(CN_fruit_,    1);
        % [13c] split total stem C by flivewd; N derived per-pool (live CN_lw / dead CN_dw)
        C_stem_total_0 = max(getInitVal(cfg, 'C_stem_wood_init', 0), 0);
        C_livestem_0   = C_stem_total_0 * flivewd_;
        C_deadstem_0   = C_stem_total_0 * (1 - flivewd_);
        N_livestem_0   = C_livestem_0 / max(CN_lw_, 1);
        N_deadstem_0   = C_deadstem_0 / max(CN_dw_, 1);
        % [coarse-root 2026-06-12] explicit live/dead coarse-root pools = croot_stem x stem,
        %   the belowground mirror of the 13c stem dual pool. croot_stem=0.225 (Zanotelli 2013
        %   Biogeosciences 10:3089 woodBG/woodAG). Faithful to CLM5 livecroot/deadcroot.
        croot_stem_    = vegetation.getPFTValue(cfg.Veg.PFT_type, 'croot_stem', 0);
        C_livecroot_0  = croot_stem_ * C_livestem_0;
        C_deadcroot_0  = croot_stem_ * C_deadstem_0;
        N_livecroot_0  = C_livecroot_0 / max(CN_lw_, 1);   % live coarse root respires (CN_lw)
        N_deadcroot_0  = C_deadcroot_0 / max(CN_dw_, 1);   % dead coarse root ~N-free (CN_dw)
        % 跨年储备 N: config 独立给 (C 侧在下方 C 块)。
        % [dual-storage 2026-06-06] Split single N_storage into coarse-root + stem reserve.
        % f_storage_root (PFT param) is the SINGLE source for the split fraction across all
        % storage sites; see getPFTParameterSet.m PFT 23 for source/units.
        f_sto_root       = vegetation.getPFTValue(cfg.Veg.PFT_type, 'f_storage_root', 0.55);
        N_storage_tot_0  = max(getInitVal(cfg, 'N_storage_init',  0), 0);
        N_root_storage_0 = N_storage_tot_0 * f_sto_root;          % [g N/m2] coarse-root reserve
        N_stem_storage_0 = N_storage_tot_0 * (1 - f_sto_root);    % [g N/m2] stem sapwood reserve
    else
        CN_L_ = 1;   % 仅 C_labile fallback (下方 C 块) 在 lNefc==1 时用; 此处保 defined 防 lint undefined
        C_leaf_0 = 0; C_root_fine_0 = 0; C_fruit_0 = 0;
        C_livestem_0 = 0; C_deadstem_0 = 0;
        N_leaf_0 = 0; N_root_fine_0 = 0; N_fruit_0 = 0;
        N_root_storage_0 = 0; N_stem_storage_0 = 0;   % [dual-storage]
        N_livestem_0 = 0; N_deadstem_0 = 0;
        C_livecroot_0 = 0; C_deadcroot_0 = 0; N_livecroot_0 = 0; N_deadcroot_0 = 0;   % [coarse-root]
    end
    NitrogenVariables.N_leaf         = zeros(tS, 1);  % [g N/m²]  LUNA 反馈用 (M2-B.3)
    NitrogenVariables.N_root_fine    = zeros(tS, 1);  % [g N/m²]  细根 N
    NitrogenVariables.N_livestem     = zeros(tS, 1);  % [g N/m²]  [13c] 活木 N (CN_lw, 呼吸)
    NitrogenVariables.N_deadstem     = zeros(tS, 1);  % [g N/m²]  [13c] 死木 N (CN_dw, 不呼吸)
    NitrogenVariables.N_livecroot    = zeros(tS, 1);  % [g N/m²]  [coarse-root] live coarse root N (respires)
    NitrogenVariables.N_deadcroot    = zeros(tS, 1);  % [g N/m²]  [coarse-root] dead coarse root N (~N-free)
    NitrogenVariables.N_fruit        = zeros(tS, 1);  % [g N/m²]  果实/种子 N (M2-B.5 harvest)
    % [dual-storage 2026-06-06] root + stem reserve pools replace single N_storage
    NitrogenVariables.N_root_storage = zeros(tS, 1);  % [g N/m²]  coarse-root seasonal reserve
    NitrogenVariables.N_stem_storage = zeros(tS, 1);  % [g N/m²]  stem sapwood seasonal reserve
    NitrogenVariables.N_leaf(1)      = N_leaf_0;       % init (cfg.Veg fallback PFT 表)
    NitrogenVariables.N_root_fine(1) = N_root_fine_0;
    NitrogenVariables.N_livestem(1)  = N_livestem_0;   % [13c]
    NitrogenVariables.N_deadstem(1)  = N_deadstem_0;   % [13c]
    NitrogenVariables.N_livecroot(1) = N_livecroot_0;  % [coarse-root]
    NitrogenVariables.N_deadcroot(1) = N_deadcroot_0;  % [coarse-root]
    NitrogenVariables.N_fruit(1)     = N_fruit_0;
    NitrogenVariables.N_root_storage(1) = N_root_storage_0;   % [spin-up] 跨年储备初值 (root)
    NitrogenVariables.N_stem_storage(1) = N_stem_storage_0;   % [spin-up] (stem)

    % --- Phase C M2-B.1b (2026-05-19): pre-allocate plantN_diag struct array ---
    % 12 fields × tS rows; lets registry reg.write detect KT=1 via .NPP empty check
    % NPP_net added 2026-06-02: unclamped GPP-R_a for true annual CUE (see computeNPP.m).
    NitrogenVariables.plantN_diag = repmat( ...
        struct('NPP',[],'NPP_net',[],'R_a_plant',[],'R_m_leaf',[],'R_m_root',[],'R_m_stem',[], ...
               'R_m_livecroot',[],'R_m_fruit',[], ...   % [Ra completeness 2026-06-14] now registered/output
               'r_m_leaf_actual',[],'r_m_root_actual',[],'r_m_stem_actual',[], ...   % [actual organ resp 2026-06-15] = theta * potential organ Rm
               'r_m_livecroot_actual',[],'r_m_fruit_actual',[],'r_g_actual',[], ...
               'ra_above_actual',[],'ra_below_actual',[], ...
               'R_g',[],'f_N_growth',[],'A_leaf',[],'A_root',[],'A_stem',[],'A_fruit',[]), ...
        tS, 1);

    % --- Phase C M2-B.4 (2026-05-21): turnover + harvest accumulators ---
    NitrogenVariables.cum_N_harvest      = 0;   % [g N/m² ground] cumulative fruit removed
    NitrogenVariables.cum_N_leaf_litter  = 0;   % [g N/m² ground] → M3 SOM 入口
    NitrogenVariables.cum_N_root_litter  = 0;   % [g N/m² ground] → M3 SOM 入口
    % [N 守恒 2026-05-30] cum_U_N: 累积根吸 N (植物 N 唯一进项, 对标 C 侧 cum_GPP)。
    %   植物 N 守恒残差 = cum_U_N − cum_N_litter − cum_N_harvest − (ΣN_pool − N_pool_init_const)。
    NitrogenVariables.cum_U_N            = 0;   % [g N/m²] cumulative root N uptake (进项)
    NitrogenVariables.turnover_state = struct('harvested_this_year', 0, 'year_now', [], 'frPHU_prev', []);
    NitrogenVariables.turnover_diag  = repmat( ...
        struct('harvest_today',[],'leaf_litter_today',[],'root_litter_today',[],'storage_release_today',[]), ...
        tS, 1);
    NitrogenVariables.frPHU_diag = zeros(tS, 1);   % [-] phenology frPHU snapshot

    % --- Phase C M2-B.8/spin-up + C-primary (v12 2026-05-29): Plant C 6 pools ---
    % 4 结构器官 C (C_leaf/root/stem/fruit) 已在上方"N 块"以 C 为主从 config 读入 (C_*_init);
    % 此处补 C_labile (NSC) 与 C_storage (跨年储备, config 独立给, 非 N 派生)。
    % [一致性 2026-05-30] C_labile 也改 C 主导 (与 4 器官同范式):
    %   C_labile_0 = config C_labile_init (主); N_labile_0 = C_labile_0 / CN_L (派生, 回填上方 N_labile)。
    %   运行时两 labile 池本就各自独立演化 (C←GPP−R_a, N←根吸 U_N), 初值统一仅为范式一致, 不动守恒/物理。
    % spin-up 末值写回 config: 4 器官更 C_*_init, 储备更 C_storage_init/N_storage_init, labile 更 C_labile_init。
    if cfg.lNefc == 1 && isprop(cfg.Veg, 'PFT_type')
        C_labile_0     = max(getInitVal(cfg, 'C_labile_init',  0), 0);   % NSC: C 主 (与 4 器官同范式)
        N_labile_0     = C_labile_0 / max(CN_L_, 1);                     % N 派生 (C 主导)
        % [dual-storage 2026-06-06] split C_storage by f_storage_root (same source as N side).
        C_storage_tot_0  = max(getInitVal(cfg, 'C_storage_init', 0), 0); % 跨年储备 (非 N 派生)
        C_root_storage_0 = C_storage_tot_0 * f_sto_root;          % [g C/m2] coarse-root reserve
        C_stem_storage_0 = C_storage_tot_0 * (1 - f_sto_root);    % [g C/m2] stem sapwood reserve
    else
        C_labile_0 = 0; N_labile_0 = 0;
        C_root_storage_0 = 0; C_stem_storage_0 = 0;   % [dual-storage]
    end
    NitrogenVariables.N_labile(1)   = N_labile_0;   % [g N/m²] C 主导反算后回填 (见上方预分配)
    NitrogenVariables.C_labile      = zeros(tS, 1);   % [g C/m²] NSC 临时仓库 (跟 N_labile 平行)
    NitrogenVariables.C_leaf        = zeros(tS, 1);   % [g C/m²] 叶 C
    NitrogenVariables.C_root_fine   = zeros(tS, 1);   % [g C/m²] 细根 C
    NitrogenVariables.C_livestem    = zeros(tS, 1);   % [g C/m²] [13c] 活木 C
    NitrogenVariables.C_deadstem    = zeros(tS, 1);   % [g C/m²] [13c] 死木 C (multi-year)
    NitrogenVariables.C_livecroot   = zeros(tS, 1);   % [g C/m²] [coarse-root] live coarse root C (respires)
    NitrogenVariables.C_deadcroot   = zeros(tS, 1);   % [g C/m²] [coarse-root] dead coarse root C (multi-year)
    NitrogenVariables.C_fruit       = zeros(tS, 1);   % [g C/m²] 果实 C
    % [dual-storage 2026-06-06] root + stem C reserve pools replace single C_storage
    NitrogenVariables.C_root_storage = zeros(tS, 1);  % [g C/m²] coarse-root seasonal reserve
    NitrogenVariables.C_stem_storage = zeros(tS, 1);  % [g C/m²] stem sapwood seasonal reserve
    % [xsmrpool 2026-06-11] maintenance-respiration deficit pool (<=0 = debt). Full R_m always paid;
    %   the shortfall is carried here and repaid from future GPP. CLM5 xsmrpool_patch
    %   (CNVegCarbonStateType.F90:60, D:/OneDrive/应用/Model/CTSM-release-clm5.0.37/src/).
    NitrogenVariables.C_xsmr = zeros(tS, 1);          % [g C/m²] xsmrpool (may be negative)
    % [E onset-xfer 2026-06-05] Transfer pools: storage mobilised at budburst, drained
    %   to display over ndays_stor. Faithful to CLM5 CNPhenologyMod storage->xfer->display.
    %   Init 0; included in the C/N conservation sums (storage->xfer->display is internal).
    NitrogenVariables.C_leaf_xfer     = zeros(tS, 1);   % [g C/m²] leaf transfer pool
    NitrogenVariables.C_root_xfer     = zeros(tS, 1);   % [g C/m²] fine-root transfer pool
    NitrogenVariables.C_livestem_xfer = zeros(tS, 1);   % [g C/m²] live-stem transfer pool
    NitrogenVariables.N_leaf_xfer     = zeros(tS, 1);   % [g N/m²] leaf N transfer pool
    NitrogenVariables.N_root_xfer     = zeros(tS, 1);   % [g N/m²] fine-root N transfer pool
    NitrogenVariables.N_livestem_xfer = zeros(tS, 1);   % [g N/m²] live-stem N transfer pool
    NitrogenVariables.onset_mobilized_flag = false;     % [-] [E] one-shot: storage already mobilised this onset (reset pre-budburst)
    NitrogenVariables.resp_night_hours = 0;             % [-] [resp-substrate EXP 2026-06-12] consecutive dark hours since sunset (bruhn scheme h)
    % [G pruning 2026-06-05] deadstem at start of year (for annual-increment pruning) + one-shot flag
    NitrogenVariables.C_deadstem_soy   = C_deadstem_0;  % [g C/m²] deadstem C at year start (Jan 1 capture)
    NitrogenVariables.N_deadstem_soy   = N_deadstem_0;  % [g N/m²] deadstem N at year start
    NitrogenVariables.pruned_this_year = false;         % [-] one-shot: pruned this autumn (reset at year start)
    NitrogenVariables.cum_C_prune_litter = 0;           % [g C/m²] cumulative pruned wood C -> litter (conservation)
    NitrogenVariables.cum_N_prune_litter = 0;           % [g N/m²] cumulative pruned wood N -> litter
    NitrogenVariables.cum_C_mortality    = 0;           % [g C/m²] [gap-mortality 2026-06-13] cumulative woody C -> litter/CWD (CNGapMortalityMod)
    NitrogenVariables.cum_N_mortality    = 0;           % [g N/m²] cumulative woody N -> litter/CWD
    NitrogenVariables.C_labile(1)     = C_labile_0;
    NitrogenVariables.C_leaf(1)       = C_leaf_0;
    NitrogenVariables.C_root_fine(1)  = C_root_fine_0;
    NitrogenVariables.C_livestem(1)   = C_livestem_0;   % [13c]
    NitrogenVariables.C_deadstem(1)   = C_deadstem_0;   % [13c]
    NitrogenVariables.C_livecroot(1)  = C_livecroot_0;  % [coarse-root]
    NitrogenVariables.C_deadcroot(1)  = C_deadcroot_0;  % [coarse-root]
    NitrogenVariables.C_fruit(1)      = C_fruit_0;
    NitrogenVariables.C_root_storage(1) = C_root_storage_0;   % [dual-storage] root reserve init
    NitrogenVariables.C_stem_storage(1) = C_stem_storage_0;   % [dual-storage] stem reserve init
    NitrogenVariables.C_xsmr(1)         = 0;                  % [xsmrpool] starts with no maintenance-resp debt (CLM5 init 0)
    % Frozen sum of initial C pools, 用于守恒残差 (不被 KT=1 步覆写)
    % [dual-storage 2026-06-06] C_storage_0 -> C_root_storage_0 + C_stem_storage_0
    NitrogenVariables.C_pool_init_const = C_labile_0 + C_leaf_0 + C_root_fine_0 + ...
                                          C_livestem_0 + C_deadstem_0 + C_fruit_0 + ...
                                          C_livecroot_0 + C_deadcroot_0 + ...   % [coarse-root]
                                          C_root_storage_0 + C_stem_storage_0 + ...
                                          0;   % [xsmrpool] C_xsmr init = 0 (included so ΔC_pool tracks the debt excursion)
    % [N 守恒 2026-05-30] Frozen sum of initial N pools, 对标 C_pool_init_const
    % [dual-storage 2026-06-06] N_storage_0 -> N_root_storage_0 + N_stem_storage_0
    NitrogenVariables.N_pool_init_const = N_labile_0 + N_leaf_0 + N_root_fine_0 + ...
                                          N_livestem_0 + N_deadstem_0 + N_fruit_0 + ...
                                          N_livecroot_0 + N_deadcroot_0 + ...   % [coarse-root]
                                          N_root_storage_0 + N_stem_storage_0;

    % M2-B.8 cumulative C transfers (送师弟 SOM 用)
    NitrogenVariables.cum_C_harvest      = 0;   % [g C/m²] 果实出地
    NitrogenVariables.cum_C_leaf_litter  = 0;   % [g C/m²] → M3 SOM
    NitrogenVariables.cum_C_root_litter  = 0;   % [g C/m²] → M3 SOM
    NitrogenVariables.cum_GPP            = 0;   % [g C/m²] 累积 GPP (守恒检查用)
    NitrogenVariables.cum_R_a            = 0;   % [g C/m²] 累积自养呼吸 actual (cap 后, 守恒严格)
    NitrogenVariables.cum_C_xsmr_writeoff = 0; % [g C/m²] crop debt reset; bookkeeping source, not CO2

    % M2-B.8b NSC 联动诊断累积器 (applyNSCCoupling 调用后填)
    NitrogenVariables.cum_R_a_potential    = 0; % [g C/m²] 累积理论自养呼吸 R_a=R_m+R_g (cap 前)
    NitrogenVariables.cum_storage_drawdown_count = 0; % [-] 全年动用 NSC 储备的 KT 步数 (非致死)
    % [nsc_diag drift fix 2026-07-04] NPP_actual_kg added: applyNSCCoupling.m:304 writes it, but this
    %   template had drifted (missing it), so runPlantCNStep.m:217 whole-struct nsc_diag(KT)=nsc_diag
    %   would hit MATLAB:heterogeneousStrucAssignment whenever nsc_diag is pre-seeded from this template.
    nsc_diag_template = struct('R_a_potential_g', 0, 'R_a_actual_g', 0, ...
        'R_a_actual_kg', 0, 'NPP_actual_kg', 0, 'storage_drawdown_event', 0, 'nsc_total_g', 0, 'nsc_net_g', 0, ...
        'c_storage_drawdown_g', 0, ...
        'A_C_leaf', 0, 'A_C_root', 0, 'A_C_stem', 0, 'A_C_fruit', 0, ...
        'C_to_root_storage', 0, 'C_to_stem_storage', 0, ...   % [N-consistency 2026-06-11] C refill (N follows C)
        'C_xsmr', 0, 'R_m_xs_g', 0, 'xsmr_recover_g', 0, ...   % [xsmrpool 2026-06-11] debt pool + fluxes (template MUST match applyNSCCoupling nsc_diag)
        'resp_beta_sub', 1, ...   % [resp-substrate EXP 2026-06-12] effective R_m turn-down (1=none)
        'resp_substrate_g', 0, 'C_labile_target_g', 0, 'C_storage_to_labile', 0, ...   % [labile sugar floor 2026-06-14] resp substrate + sugar-floor diags
        'resp_theta_actual', 1);   % [actual organ resp 2026-06-15] effective maintenance throttle = actual R_m/potential R_m (1=none)
    NitrogenVariables.nsc_diag = repmat(nsc_diag_template, tS, 1);

end

% Change: 2026-09-02 Jizixing Zhu -- add the crop xsmr debt-writeoff C ledger.

% ============================================================================
% Helper: 优先读 cfg.Veg.(name), fallback PFT 表 vegetation.getPFTValue, 再 fallback default.
%
%   v = getInitVal(cfg, 'C_leaf_init', 0)
%
%   解 spin-up 工作流: 写回 config_vegetation.txt 不动 PFT 表 (物种字典).
%   cfg.Veg.(name) 在 VegetationData properties 中默认 [] = 未设, 走 PFT 表 fallback.
% ============================================================================
function v = getInitVal(cfg, name, default)
    v = default;
    if isprop(cfg.Veg, name)
        val = cfg.Veg.(name);
        if ~isempty(val) && isnumeric(val) && isfinite(val)
            v = val;
            return;
        end
    end
    if isprop(cfg.Veg, 'PFT_type')
        try
            v = vegetation.getPFTValue(cfg.Veg.PFT_type, name);
            if isempty(v) || ~isfinite(v), v = default; end
        catch
            % PFT 表没这个字段 (如 C_storage_init), 用 default
            v = default;
        end
    end
end

% Change: 2026-09-30 Wenjie Liu -- keep the compatibility NSC diagnostic
% template consistent with annual WSC and standing-senescence diagnostics.
% Change: 2026-09-30 Wenjie Liu -- withdraw this round's APSIM stem WSC
%   and standing-senescent-leaf extension; retain the original COMPASS litter
%   and CoLM crop C/N allocation/harvest contracts. Earlier Change records are historical.
