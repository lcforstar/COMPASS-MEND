function demo_v2_step_response()
%DEMO_V2_STEP_RESPONSE  瞬态版本阶跃响应演示 / Step-response demo
%
% 不依赖任何主循环，独立演示 #6 (冠层热容) 和 #7 (植物水力) 的物理行为。
% 输出 6 张图：
%   1. 冠层温度对辐射阶跃的响应
%   2. 冠层温度对 VPD 阶跃的响应
%   3. 不同冠层热容下的时标对比
%   4. 植物水量对干旱的响应
%   5. 叶水势昼夜节律
%   6. 守恒诊断：累积通量 vs 储量变化
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    fprintf('\n=== Demo 1: 冠层温度对辐射阶跃 ===\n');
    demo_canopy_radiation_step();

    fprintf('\n=== Demo 2: 冠层热容对比 ===\n');
    demo_canopy_thermal_inertia();

    fprintf('\n=== Demo 3: 植物水力干旱响应 ===\n');
    demo_plant_drought();

    fprintf('\n=== Demo 4: 守恒验证 ===\n');
    demo_conservation();

    fprintf('\n所有演示完成 / All demos completed.\n');
end

% ============================================================
% Demo 1: 冠层温度对辐射阶跃
% ============================================================
function demo_canopy_radiation_step()
    % 参数
    params = default_canopy_params();
    Delt_t = 60;   % 1 min
    n_steps = 60;  % 60 min 总时长

    % 初始平衡态：T = T_air, Rn = H + LE
    T_old = 290;
    T_air = 290;
    Rn_baseline = 200;
    H_baseline = 60;
    LE_baseline = 140;

    % 阶跃：t=10 min 时 Rn 突增到 500 W/m²
    T_history = zeros(n_steps, 1);
    T_history(1) = T_old;
    for k = 2:n_steps
        if k < 10
            Rn = Rn_baseline;
        else
            Rn = 500;   % 阶跃
        end

        % 显热/潜热（简化：跟随 T - T_air）
        H = (T_old - T_air) * 1005 * 1.2 / 30 + H_baseline;
        LE = LE_baseline;   % 暂不让 LE 跟随，看纯热惯性

        [T_new, dT, info] = vegetation.hydraulics.prognosticCanopyEnergy( ...
            T_old, T_old, Rn, H, LE, 0, ...
            params.rho_c_h, params.A_canopy, ...
            params.rho_a, params.cp, params.Lv, ...
            params.r_ah, params.r_ae, ...
            params.qsat, T_old, params.epsilon, params.sigmaSB, params.R_v, ...
            Delt_t);

        T_history(k) = T_new;
        T_old = T_new;
    end

    fprintf('  时标 τ_canopy = %.1f s (理论 %.1f s)\n', ...
        info.tau_canopy, params.rho_c_h * params.A_canopy / 40);
    fprintf('  阶跃前 T = %.2f K, 阶跃后 1 h 末 T = %.2f K\n', ...
        T_history(1), T_history(end));
    fprintf('  温度上升量 ΔT = %.2f K\n', T_history(end) - T_history(1));

    % 简单文本可视化
    print_simple_plot(T_history, 'Canopy T (K)', n_steps);
end

% ============================================================
% Demo 2: 冠层热容对比 (不同 ρch 时标差异)
% ============================================================
function demo_canopy_thermal_inertia()
    rho_c_h_list = [100, 745, 5000, 65000];   % 极轻 / 默认 / 中等 / 苹果园(含果)
    labels = {'轻量(100)', '默认(745)', '中等(5000)', '苹果园(65000)'};

    Delt_t = 60;
    n_steps = 120;

    for i = 1:length(rho_c_h_list)
        params = default_canopy_params();
        params.rho_c_h = rho_c_h_list(i);

        T_old = 290;
        T_history = zeros(n_steps, 1);
        T_history(1) = T_old;

        for k = 2:n_steps
            if k < 10, Rn = 200; else, Rn = 500; end

            [T_new, ~, info] = vegetation.hydraulics.prognosticCanopyEnergy( ...
                T_old, T_old, Rn, 60, 140, 0, ...
                params.rho_c_h, params.A_canopy, ...
                params.rho_a, params.cp, params.Lv, ...
                params.r_ah, params.r_ae, ...
                params.qsat, T_old, params.epsilon, params.sigmaSB, params.R_v, ...
                Delt_t);
            T_history(k) = T_new;
            T_old = T_new;
        end

        % 找到达到 63% 平衡的时间
        T_steady = T_history(end);
        T_init = T_history(1);
        T_63 = T_init + 0.63 * (T_steady - T_init);
        idx_63 = find(T_history >= T_63, 1);
        if isempty(idx_63), idx_63 = NaN; end

        fprintf('  ρch=%-5d %s: τ=%.1f s, 实测 63%% 时间 ≈ %d min\n', ...
            rho_c_h_list(i), labels{i}, info.tau_canopy, ...
            (idx_63 - 9));
    end
end

% ============================================================
% Demo 3: 植物水力干旱响应 (展示储水缓冲)
% ============================================================
function demo_plant_drought()
    plant_params = struct( ...
        'K_sat',       2e-8, ...
        'psi_63',      -275, ...
        'c_weibull',   4, ...
        'a_retention', [4.9e-4; 4.9e-4; 4.9e-4], ...
        'W_capacity',  [0.003; 0.005; 0.001], ...
        'dz_rs',       0.5, ...
        'dz_sl',       1.0);

    % 初始：完全饱和
    W_old = plant_params.W_capacity * 0.95;

    NN = 10;
    P_root = ones(NN, 1) / NN;
    K_sr = 1e-8 * ones(NN, 1);   % 较干燥
    psi_soil_dry = -300 * ones(NN, 1);   % 强吸力
    T_pot = 5e-8;   % 持续蒸腾压力

    Delt_t = 3600;
    n_days = 14;
    n_steps = n_days * 24;

    W_hist = zeros(3, n_steps);
    psi_hist = zeros(3, n_steps);
    stress_hist = zeros(1, n_steps);
    T_actual_hist = zeros(1, n_steps);

    for t = 1:n_steps
        [W_old, psi_new, q_fluxes, info] = ...
            vegetation.hydraulics.prognosticPlantWater( ...
                W_old, psi_soil_dry, P_root, K_sr, T_pot, ...
                plant_params, Delt_t);
        W_hist(:, t) = W_old;
        psi_hist(:, t) = psi_new;
        stress_hist(t) = info.stress_factor;
        T_actual_hist(t) = q_fluxes.T_actual;
    end

    % 找胁迫开始时间
    idx_stress = find(stress_hist < 0.5, 1);
    if isempty(idx_stress), idx_stress = n_steps; end
    days_to_stress = idx_stress / 24;

    fprintf('  起始 W_total = %.3f mm\n', sum(W_hist(:,1))*1000);
    fprintf('  %d 天后 W_total = %.3f mm\n', n_days, sum(W_hist(:,end))*1000);
    fprintf('  达到胁迫(stress<0.5) 时间: %.1f 天\n', days_to_stress);
    fprintf('  实际蒸腾比 (T_actual/T_pot) 末态: %.1f%%\n', ...
        T_actual_hist(end)/T_pot*100);

    print_simple_plot(stress_hist, '水力胁迫因子', n_steps);
end

% ============================================================
% Demo 4: 守恒验证
% ============================================================
function demo_conservation()
    % --- 冠层能量守恒 ---
    params = default_canopy_params();
    Delt_t = 60;
    n_steps = 1440;   % 1 天

    T_old = 290;
    T_init = T_old;

    cum_energy_in = 0;
    for k = 1:n_steps
        % 模拟日变化辐射
        Rn = 200 + 300 * sin(2*pi*k/1440);
        H = 60;
        LE = 140;

        [T_new, dT, info] = vegetation.hydraulics.prognosticCanopyEnergy( ...
            T_old, T_old, Rn, H, LE, 0, ...
            params.rho_c_h, params.A_canopy, ...
            params.rho_a, params.cp, params.Lv, ...
            params.r_ah, params.r_ae, ...
            params.qsat, T_old, params.epsilon, params.sigmaSB, params.R_v, ...
            Delt_t);

        cum_energy_in = cum_energy_in + (Rn - H - LE) * Delt_t;
        T_old = T_new;
    end
    storage_change = params.rho_c_h * params.A_canopy * (T_old - T_init);
    err_canopy = (storage_change - cum_energy_in) / max(abs(cum_energy_in), 1);
    fprintf('  冠层能量守恒：储能变化 = %.2f J/m², 累积通量 = %.2f J/m²\n', ...
        storage_change, cum_energy_in);
    fprintf('  相对误差: %.4f%%\n', err_canopy * 100);

    % --- 植物水量守恒 ---
    plant_params = struct( ...
        'K_sat',       2e-8, ...
        'psi_63',      -275, ...
        'c_weibull',   4, ...
        'a_retention', [4.9e-4; 4.9e-4; 4.9e-4], ...
        'W_capacity',  [0.003; 0.005; 0.001], ...
        'dz_rs',       0.5, ...
        'dz_sl',       1.0);

    W_old = plant_params.W_capacity * 0.8;
    W_init = W_old;

    NN = 10;
    P_root = ones(NN, 1) / NN;
    K_sr = 1e-7 * ones(NN, 1);
    psi_soil = -1 * ones(NN, 1);
    T_pot = 5e-8;
    Delt_t = 3600;
    n_steps = 24;   % 1 天

    cum_q_in = 0;
    cum_T = 0;
    for t = 1:n_steps
        [W_old, ~, q, ~] = vegetation.hydraulics.prognosticPlantWater( ...
            W_old, psi_soil, P_root, K_sr, T_pot, plant_params, Delt_t);
        cum_q_in = cum_q_in + q.q_soil_to_root * Delt_t;
        cum_T = cum_T + q.T_actual * Delt_t;
    end
    storage_change_w = sum(W_old - W_init);
    expected_change = cum_q_in - cum_T;
    err_water = (storage_change_w - expected_change);
    fprintf('  植物水量守恒：储水变化 = %.6f m, 累积通量 = %.6f m\n', ...
        storage_change_w, expected_change);
    fprintf('  绝对误差: %.2e m\n', err_water);
    if abs(err_water) < 1e-6
        fprintf('  ✓ 守恒通过\n');
    else
        fprintf('  ⚠ 守恒误差超阈\n');
    end
end

% ============================================================
% 辅助：默认冠层参数
% ============================================================
function params = default_canopy_params()
    params = struct( ...
        'rho_c_h',   745, ...
        'A_canopy',  3.0, ...
        'rho_a',     1.2, ...
        'cp',        1005, ...
        'Lv',        2.45e6, ...
        'r_ah',      30, ...
        'r_ae',      30, ...
        'qsat',      0.012, ...
        'epsilon',   0.97, ...
        'sigmaSB',   5.67e-8, ...
        'R_v',       461.5);
end

% ============================================================
% 辅助：简单文本绘图
% ============================================================
function print_simple_plot(data, label, n_max)
    if length(data) < 2
        return;
    end
    n_show = min(50, length(data));
    step = max(1, floor(length(data) / n_show));
    sub = data(1:step:end);

    minv = min(sub);
    maxv = max(sub);
    if maxv == minv, return; end

    fprintf('  %s 演化（前→后）:\n', label);
    width = 50;
    for i = 1:length(sub)
        bar_len = round((sub(i) - minv) / (maxv - minv) * width);
        bar = repmat('█', 1, bar_len);
        if mod(i, max(1, floor(length(sub)/10))) == 1 || i == length(sub)
            fprintf('    %5.2f |%s\n', sub(i), bar);
        end
    end
end
