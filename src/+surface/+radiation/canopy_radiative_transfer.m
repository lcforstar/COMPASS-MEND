function [I_up, I_down, absorbed, reflected, transmitted] = canopy_radiative_transfer(L, S, mu, G_mu, omega, beta, beta0, alpha_g_mu, alpha_g, I_sky_diffuse, varargin)
%CANOPY_RADIATIVE_TRANSFER_COMPLETE 完整的植被冠层辐射传输两流模型
%   包含正确的大气漫射辐射上边界条件和近似值计算
%
%   输入参数:
%     L: 叶面积指数 [无量纲]
%     S: 茎面积指数 [无量纲] 
%     mu: 太阳天顶角余弦 [无量纲]
%     G_mu: 叶片投影函数 [无量纲]
%     omega: 散射系数 [无量纲]
%     beta: 漫射光向上散射参数 [无量纲]
%     beta0: 直射光向上散射参数 [无量纲]
%     alpha_g_mu: 地表直射反照率 [无量纲]
%     alpha_g: 地表漫射反照率 [无量纲]
%     I_sky_diffuse: 大气漫射辐射通量 [单位入射通量]
%     varargin: 可选参数
%
%   输出参数:
%     I_up: 向上漫射辐射通量（每层）[无量纲]
%     I_down: 向下漫射辐射通量（每层）[无量纲]
%     absorbed: 植被吸收的总辐射 [无量纲]
%     reflected: 冠层反射的总辐射 [无量纲]
%     transmitted: 冠层透射的总辐射 [无量纲]

% ======================== 参数验证和默认值 ========================
if nargin < 10
    % 如果没有提供大气漫射辐射，使用近似值
    I_sky_diffuse = estimate_diffuse_sky_radiation(mu);
end

% 验证输入参数
validateattributes(L, {'numeric'}, {'scalar', '>=', 0});
validateattributes(S, {'numeric'}, {'scalar', '>=', 0});
validateattributes(mu, {'numeric'}, {'scalar', '>', 0, '<=', 1});
validateattributes(omega, {'numeric'}, {'scalar', '>=', 0, '<=', 1});
validateattributes(beta, {'numeric'}, {'scalar', '>=', 0, '<=', 1});
validateattributes(beta0, {'numeric'}, {'scalar', '>=', 0, '<=', 1});
validateattributes(I_sky_diffuse, {'numeric'}, {'scalar', '>=', 0});

% 设置默认参数
N_layers = 100; % 分层数
tau_total = L + S; % 总光学深度
K = G_mu / mu; % 直射光束的光学深度
mu_bar = 0.5; % 平均漫射光学深度倒数（典型值）

% ======================== 网格设置 ========================
tau = linspace(0, tau_total, N_layers); % 光学深度坐标
d_tau = tau_total / (N_layers - 1); % 层间光学深度增量

% ======================== 初始化通量数组 ========================
I_up = zeros(1, N_layers); % 向上漫射通量
I_down = zeros(1, N_layers); % 向下漫射通量

% ======================== 正确的边界条件设置 ========================
% 上边界条件 (tau = 0): 来自大气的漫射辐射
I_down(1) = I_sky_diffuse; % 关键修正：大气漫射辐射输入

% 下边界条件 (tau = tau_total): 地表反射
% 将在迭代过程中动态更新

% ======================== 有限差分求解 ========================
max_iter = 1000;
tolerance = 1e-6;

fprintf('开始求解两流方程...\n');
fprintf('上边界条件: I_sky = %.4f\n', I_sky_diffuse);

for iter = 1:max_iter
    I_up_old = I_up;
    I_down_old = I_down;
    
    % 从底部向上更新向上通量 (i = N-1 到 2)
    for i = N_layers-1:-1:2
        % 方程(3.1)离散化: -μ̄ dI↑/dτ + [1-(1-β)ω]I↑ - ωβI↓ = ωμ̄Kβ₀e^{-Kτ}
        
        % 中心差分近似导数
        if i == 1
            dI_up_dtau = (I_up(i+1) - I_up(i)) / d_tau; % 前向差分（上边界）
        elseif i == N_layers
            dI_up_dtau = (I_up(i) - I_up(i-1)) / d_tau; % 后向差分（下边界）
        else
            dI_up_dtau = (I_up(i+1) - I_up(i-1)) / (2 * d_tau); % 中心差分
        end
        
        % 源项：直射光散射贡献
        source_term = omega * mu_bar * K * beta0 * exp(-K * tau(i));
        
        % 更新I_up (重新整理方程3.1)
        numerator = source_term + mu_bar * dI_up_dtau + omega * beta * I_down(i);
        denominator = 1 - (1 - beta) * omega;
        
        if denominator ~= 0
            I_up(i) = numerator / denominator;
        end
    end
    
    % 从顶部向下更新向下通量 (i = 2 到 N-1)
    for i = 2:N_layers-1
        % 方程(3.2)离散化: μ̄ dI↓/dτ + [1-(1-β)ω]I↓ - ωβI↑ = ωμ̄K(1-β₀)e^{-Kτ}
        
        % 中心差分近似导数
        if i == 1
            dI_down_dtau = (I_down(i+1) - I_down(i)) / d_tau; % 前向差分
        elseif i == N_layers
            dI_down_dtau = (I_down(i) - I_down(i-1)) / d_tau; % 后向差分
        else
            dI_down_dtau = (I_down(i+1) - I_down(i-1)) / (2 * d_tau); % 中心差分
        end
        
        % 源项：直射光散射贡献
        source_term = omega * mu_bar * K * (1 - beta0) * exp(-K * tau(i));
        
        % 更新I_down (重新整理方程3.2)
        numerator = source_term - mu_bar * dI_down_dtau + omega * beta * I_up(i);
        denominator = 1 - (1 - beta) * omega;
        
        if denominator ~= 0
            I_down(i) = numerator / denominator;
        end
    end
    
    % 应用边界条件
    % 上边界保持: I_down(1) = I_sky_diffuse (不变)
    
    % 下边界条件: I↑(τ_total) = α_g^μ e^{-Kτ_total} + α_g I↓(τ_total)
    I_up(end) = alpha_g_mu * exp(-K * tau_total) + alpha_g * I_down(end);
    
    % 检查收敛
    max_change_up = max(abs(I_up - I_up_old));
    max_change_down = max(abs(I_down - I_down_old));
    
    if max_change_up < tolerance && max_change_down < tolerance
        fprintf('收敛于第 %d 次迭代, 最大变化: %.2e\n', iter, max(max_change_up, max_change_down));
        break;
    end
    
    if iter == max_iter
        warning('未达到收敛容差 (最大变化: %.2e)', max(max_change_up, max_change_down));
    end
end

% ======================== 计算辐射平衡 ========================
% 直射光贡献
direct_beam_top = 1.0; % 单位入射直射通量
direct_beam_bottom = direct_beam_top * exp(-K * tau_total);

% 植被吸收的辐射
% 直射光吸收 + 漫射光吸收
absorbed_direct = direct_beam_top - direct_beam_bottom;
absorbed_diffuse = trapz(tau, (1 - omega) * (I_up + I_down));
absorbed = absorbed_direct + absorbed_diffuse;

% 冠层反射的辐射（顶部向上通量）
reflected = I_up(1);

% 冠层透射的辐射
% 直射光透射 + 漫射光透射
transmitted_direct = direct_beam_bottom;
transmitted_diffuse = I_down(end);
transmitted = transmitted_direct + transmitted_diffuse;

% ======================== 结果显示 ========================
display_results_complete(tau, I_up, I_down, tau_total, I_sky_diffuse, ...
    absorbed, reflected, transmitted, direct_beam_bottom);

end

function I_sky = estimate_diffuse_sky_radiation(mu)
%ESTIMATE_DIFFUSE_SKY_RADIATION 估算大气漫射辐射（近似值）
%   基于经验公式估算大气漫射辐射与直射辐射的比例
%   参考: Liu and Jordan (1960), Erbs et al. (1982) 等经验模型
%
%   输入:
%     mu: 太阳天顶角余弦
%
%   输出:
%     I_sky: 大气漫射辐射（单位直射辐射）

    % 晴空指数近似 (Kt = 漫射/总辐射)
    % 对于晴空条件，典型值范围
    if mu > 0.7 % 太阳高度角较大
        I_sky = 0.1 + 0.1 * (1 - mu); % 较小比例的漫射辐射
    elseif mu > 0.3
        I_sky = 0.2 + 0.2 * (0.7 - mu);
    else % 太阳高度角较小
        I_sky = 0.3 + 0.3 * (0.3 - mu); % 较大比例的漫射辐射
    end
    
    % 确保合理性
    I_sky = max(0.05, min(0.5, I_sky));
end

function display_results_complete(tau, I_up, I_down, tau_total, I_sky, ...
    absorbed, reflected, transmitted, direct_trans,G_mu,mu)
%DISPLAY_RESULTS_COMPLETE 显示完整结果
    figure('Position', [100, 100, 1200, 500]);
    
    % 子图1: 漫射辐射通量分布
    subplot(1,3,1);
    plot(I_up, tau, 'r-', 'LineWidth', 2); hold on;
    plot(I_down, tau, 'b-', 'LineWidth', 2);
    plot([0, max(I_down)*1.1], [0, 0], 'k--', 'LineWidth', 1); % 顶部边界
    plot([0, max(I_up)*1.1], [tau_total, tau_total], 'k--', 'LineWidth', 1); % 底部边界
    
    % 标记边界条件
    plot(I_sky, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(I_sky, -0.02, sprintf('I_{sky}=%.3f', I_sky), 'HorizontalAlignment', 'center');
    
    set(gca, 'YDir', 'reverse');
    ylabel('光学深度 \tau');
    xlabel('漫射辐射通量');
    title('冠层内漫射辐射分布');
    legend('向上通量 I↑', '向下通量 I↓', '冠层边界', 'Location', 'best');
    grid on;
    ylim([-0.05, tau_total+0.05]);
    
    % 子图2: 直射光衰减
    subplot(1,3,2);
    direct_beam = exp(-(G_mu/mu) * tau);
    plot(direct_beam, tau, 'g-', 'LineWidth', 2);
    set(gca, 'YDir', 'reverse');
    ylabel('光学深度 \tau');
    xlabel('直射辐射通量');
    title('直射光衰减');
    text(0.5, tau_total/2, sprintf('透射率: %.3f', direct_trans), 'FontSize', 12);
    grid on;
    
    % 子图3: 辐射平衡
    subplot(1,3,3);
    components = [absorbed, reflected, transmitted];
    labels = {'吸收', '反射', '透射'};
    pie(components, labels);
    title('辐射平衡');
    
    % 数值结果
    fprintf('\n=== 完整辐射传输结果 ===\n');
    fprintf('大气漫射输入 (I_sky): %.4f\n', I_sky);
    fprintf('植被吸收: %.4f (%.1f%%)\n', absorbed, absorbed*100);
    fprintf('冠层反射: %.4f (%.1f%%)\n', reflected, reflected*100);
    fprintf('冠层透射: %.4f (%.1f%%)\n', transmitted, transmitted*100);
    fprintf('   - 直射透射: %.4f\n', direct_trans);
    fprintf('   - 漫射透射: %.4f\n', transmitted - direct_trans);
    fprintf('能量平衡: 吸收 + 反射 + 透射 = %.4f\n', absorbed + reflected + transmitted);
    fprintf('平衡误差: %.6f\n', 1 - (absorbed + reflected + transmitted));
end

