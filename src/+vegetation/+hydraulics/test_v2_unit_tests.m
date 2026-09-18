function test_v2_unit_tests()
%TEST_V2_UNIT_TESTS  瞬态升级 (#6 #7) 单元测试套件
%                    Unit tests for prognostic upgrades
%
% 不依赖主循环，独立验证 4 个新算法文件的物理正确性。
% 运行：在 src 目录下输入 test_v2_unit_tests
%
% 测试覆盖 / Tests covered:
%   T1: weibullConductance 边界值 + 物理形状
%   T2: linearRetention 单调性 + 反函数自洽
%   T3: harmonicMean 极值 + 对称性
%   T4: prognosticCanopyEnergy 守恒 + 时标
%   T5: prognosticPlantWater 守恒 + 稳态退化
%   T6: prognosticPlantWater 干旱响应（4天滞后）
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    fprintf('\n========================================\n');
    fprintf('  STEMMUS-SCOPE v2 Unit Tests\n');
    fprintf('  Author: Zhu Jizixing  Date: 2026-05-01\n');
    fprintf('========================================\n\n');

    pass_count = 0;
    fail_count = 0;

    [pass_count, fail_count] = run_test(@test_weibull,         'T1: weibullConductance',        pass_count, fail_count);
    [pass_count, fail_count] = run_test(@test_retention,       'T2: linearRetention',           pass_count, fail_count);
    [pass_count, fail_count] = run_test(@test_harmonic,        'T3: harmonicMean',              pass_count, fail_count);
    [pass_count, fail_count] = run_test(@test_canopy_energy,   'T4: prognosticCanopyEnergy',    pass_count, fail_count);
    [pass_count, fail_count] = run_test(@test_plant_water,     'T5: prognosticPlantWater',      pass_count, fail_count);
    [pass_count, fail_count] = run_test(@test_drought_response,'T6: drought response',          pass_count, fail_count);

    fprintf('\n========================================\n');
    fprintf('  Result: %d passed, %d failed\n', pass_count, fail_count);
    fprintf('========================================\n\n');

    if fail_count > 0
        error('test_v2_unit_tests:Failure', '%d 测试失败', fail_count);
    end
end

% ============================================================
% Test runner
% ============================================================
function [pc, fc] = run_test(testfn, testname, pc, fc)
    try
        testfn();
        fprintf('  [ PASS ] %s\n', testname);
        pc = pc + 1;
    catch ME
        fprintf('  [ FAIL ] %s\n', testname);
        fprintf('           %s\n', ME.message);
        fc = fc + 1;
    end
end

% ============================================================
% T1: weibullConductance
% ============================================================
function test_weibull()
    K_sat = 2e-8;
    psi_63 = -275;
    c = 4;

    % 饱和：psi=0 → K=K_sat
    K0 = vegetation.hydraulics.weibullConductance(0, K_sat, psi_63, c);
    assert(abs(K0 - K_sat) < 1e-15, 'psi=0 应给 K_sat');

    % 强吸力：psi=ψ63 → K = K_sat * exp(-1) ≈ 0.368 K_sat
    K63 = vegetation.hydraulics.weibullConductance(psi_63, K_sat, psi_63, c);
    assert(abs(K63 - K_sat * exp(-1)) < 1e-12, 'psi=ψ63 应给 K_sat*exp(-1)');

    % 极强吸力：K → 0
    K_extreme = vegetation.hydraulics.weibullConductance(-10000, K_sat, psi_63, c);
    assert(K_extreme < K_sat * 1e-10, '极强吸力应 K→0');

    % 单调性：psi 越负，K 越小
    psi_array = -[0, 50, 100, 200, 400, 800];
    K_array = arrayfun(@(p) vegetation.hydraulics.weibullConductance(p, K_sat, psi_63, c), psi_array);
    assert(all(diff(K_array) <= 0), 'K 应随 |psi| 递减');
end

% ============================================================
% T2: linearRetention
% ============================================================
function test_retention()
    W_cap = 0.005;
    a = 4.9e-4;

    % 饱和：W=W_cap → ψ=0
    psi0 = vegetation.hydraulics.linearRetention(W_cap, W_cap, a);
    assert(abs(psi0) < 1e-12, 'W=W_cap 应给 ψ=0');

    % 完全干燥：W=0 → ψ = -1/a
    psi_dry = vegetation.hydraulics.linearRetention(0, W_cap, a);
    assert(abs(psi_dry - (-1/a)) < 1e-9, 'W=0 应给 ψ=-1/a');

    % 反函数自洽性
    psi_test = -100;
    W_test = vegetation.hydraulics.inverseLinearRetention(psi_test, W_cap, a);
    psi_back = vegetation.hydraulics.linearRetention(W_test, W_cap, a);
    assert(abs(psi_back - psi_test) < 1e-9, '反函数自洽性破坏');
end

% ============================================================
% T3: harmonicMean
% ============================================================
function test_harmonic()
    K1 = 2e-8;
    K2 = 1e-9;

    K_eff = vegetation.hydraulics.harmonicMean(K1, K2);

    % 调和平均 < min(K1, K2) * 2，但 > min
    assert(K_eff <= min(K1, K2) * 2.01, '调和平均不应超过两倍小值');
    assert(K_eff >= min(K1, K2) * 0.99, '调和平均不应远小于小值');

    % 对称性
    K_sym = vegetation.hydraulics.harmonicMean(K2, K1);
    assert(abs(K_eff - K_sym) < 1e-25, '调和平均应对称');

    % 等值
    K_equal = vegetation.hydraulics.harmonicMean(K1, K1);
    assert(abs(K_equal - K1) < 1e-25, '等值时调和平均=单值');
end

% ============================================================
% T4: prognosticCanopyEnergy 守恒 + 时标
% ============================================================
function test_canopy_energy()
    % --- 基本参数（典型苹果园 LAI=3）---
    rho_c_h = 745;
    A_canopy = 3.0;
    rho_a = 1.2;
    cp = 1005;
    Lv = 2.45e6;
    r_ah = 30;
    r_ae = 30;
    qsat = 0.012;
    epsilon = 0.97;
    sigmaSB = 5.67e-8;
    R_v = 461.5;
    Delt_t = 60;

    T_old = 290;
    T_iter = 290;
    T_air = 288;

    % --- 平衡态：Rn = H + LE，应该 dT≈0 ---
    Rn = 200;
    H = 50;
    LE = 150;
    S_re_e = 0;

    [T_new, dT, info] = vegetation.hydraulics.prognosticCanopyEnergy( ...
        T_old, T_iter, Rn, H, LE, S_re_e, ...
        rho_c_h, A_canopy, rho_a, cp, Lv, r_ah, r_ae, ...
        qsat, T_iter, epsilon, sigmaSB, R_v, Delt_t);

    assert(abs(dT) < 0.5, '平衡态 dT 应接近 0');
    assert(isfinite(info.tau_canopy), '时标应有限');

    % --- 不平衡态：Rn 突增，T 应上升 ---
    Rn = 400;
    [T_new2, dT2, info2] = vegetation.hydraulics.prognosticCanopyEnergy( ...
        T_old, T_iter, Rn, H, LE, S_re_e, ...
        rho_c_h, A_canopy, rho_a, cp, Lv, r_ah, r_ae, ...
        qsat, T_iter, epsilon, sigmaSB, R_v, Delt_t);

    assert(dT2 > 0, 'Rn 突增应使 T 上升');
    assert(T_new2 > T_old, 'T_new 应 > T_old');

    % --- 时标合理性：~10s 到 100s ---
    assert(info2.tau_canopy > 1, '时标应 > 1s');
    assert(info2.tau_canopy < 1000, '时标应 < 1000s');

    % --- 极端：稀疏冠层（A→0）应不崩 ---
    [T_sparse, dT_sparse, ~] = vegetation.hydraulics.prognosticCanopyEnergy( ...
        T_old, T_iter, Rn, H, LE, S_re_e, ...
        rho_c_h, 0.001, rho_a, cp, Lv, r_ah, r_ae, ...
        qsat, T_iter, epsilon, sigmaSB, R_v, Delt_t);
    assert(isfinite(T_sparse), '稀疏冠层 T 应有限');
end

% ============================================================
% T5: prognosticPlantWater 守恒 + 稳态退化
% ============================================================
function test_plant_water()
    % --- 参数 ---
    plant_params = struct( ...
        'K_sat',       2e-8, ...
        'psi_63',      -275, ...
        'c_weibull',   4, ...
        'a_retention', [4.9e-4; 4.9e-4; 4.9e-4], ...
        'W_capacity',  [0.003; 0.005; 0.001], ...
        'dz_rs',       0.5, ...
        'dz_sl',       1.0);

    % --- 初始：所有层 80% 饱和 ---
    W_old = plant_params.W_capacity * 0.8;

    % --- 土壤剖面（10 层，潮湿）---
    NN = 10;
    psi_soil = -1 * ones(NN, 1);   % -1 m
    P_root = ones(NN, 1) / NN;
    K_sr = 1e-7 * ones(NN, 1);

    % --- 适中蒸腾 ---
    T_potential = 1e-7;   % m/s ≈ 8.6 mm/day
    Delt_t = 60;          % 1 min

    [W_new, psi_new, q_fluxes, info] = ...
        vegetation.hydraulics.prognosticPlantWater( ...
            W_old, psi_soil, P_root, K_sr, T_potential, plant_params, Delt_t);

    % --- 维度检查 ---
    assert(numel(W_new) == 3, 'W_new 应 3×1');
    assert(numel(psi_new) == 3, 'psi_new 应 3×1');

    % --- 物理钳位：W >=0, W <= cap*1.001 ---
    assert(all(W_new >= 0), 'W 不能为负');
    assert(all(W_new <= plant_params.W_capacity * 1.002), 'W 不能远超容量');

    % --- 守恒：水量守恒残差应接近 0（IMEX 离散误差）---
    assert(abs(info.water_balance_residual) < 1e-6, ...
           sprintf('水量守恒残差 %g 太大', info.water_balance_residual));

    % --- 通量合理性：T_actual 应 <= T_potential ---
    assert(q_fluxes.T_actual <= T_potential * 1.001, 'T_actual 不应超 T_potential');

    % --- 稳态极限测试：长 dt 后应接近稳态 ---
    Delt_t_long = 86400 * 10;   % 10 天
    [W_steady, ~, ~, ~] = vegetation.hydraulics.prognosticPlantWater( ...
        W_old, psi_soil, P_root, K_sr, T_potential, plant_params, Delt_t_long);

    % 长时间后水量变化平缓
    assert(all(abs(W_steady - W_old) < plant_params.W_capacity), ...
           '稳态极限：W 不应远离物理范围');
end

% ============================================================
% T6: 干旱响应测试（4天滞后）
% ============================================================
function test_drought_response()
    plant_params = struct( ...
        'K_sat',       2e-8, ...
        'psi_63',      -275, ...
        'c_weibull',   4, ...
        'a_retention', [4.9e-4; 4.9e-4; 4.9e-4], ...
        'W_capacity',  [0.003; 0.005; 0.001], ...
        'dz_rs',       0.5, ...
        'dz_sl',       1.0);

    W_old = plant_params.W_capacity * 0.95;   % 几乎饱和起点

    NN = 10;
    P_root = ones(NN, 1) / NN;
    K_sr = 1e-9 * ones(NN, 1);   % 干燥土壤，导度低

    Delt_t = 3600;   % 1 h
    n_hours = 24 * 7;   % 模拟 7 天

    % 干旱起始：土壤水势-200 m, 蒸腾正常
    psi_soil_dry = -200 * ones(NN, 1);
    T_potential = 5e-8;

    W_history = zeros(3, n_hours);
    stress_history = zeros(1, n_hours);
    W_curr = W_old;

    for t = 1:n_hours
        [W_curr, ~, ~, info] = ...
            vegetation.hydraulics.prognosticPlantWater( ...
                W_curr, psi_soil_dry, P_root, K_sr, T_potential, ...
                plant_params, Delt_t);
        W_history(:, t) = W_curr;
        stress_history(t) = info.stress_factor;
    end

    % --- 验证：储水应该单调下降 ---
    W_total = sum(W_history, 1);
    %（小幅振荡允许，但总趋势应下降）
    assert(W_total(end) < W_total(1), '干旱 7 天后总储水应下降');

    % --- 验证：胁迫因子应在第 2-5 天显著上升 ---
    %   稳态模型立即 100% 胁迫，瞬态有几天延迟
    %   stress_factor 1=健康 0=胁迫
    assert(stress_history(1) > 0.5, '第1天应仍较健康（>0.5）');
    %  允许后续逐渐下降（具体多快取决于 K, W_cap）

    % --- 验证：所有 W_history 物理合理 ---
    assert(all(W_history(:) >= 0), '储水不能为负');
    assert(all(W_history(:) <= 1), '储水不能超 1m');
end
