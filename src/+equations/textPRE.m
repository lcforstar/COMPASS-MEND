% --- 测试脚本：1000次【完全随机】降水模拟分析 ---
clear; clc;

N_sim = 1000;         % 模拟次数
R_day = 50;           % 日降雨量 50 mm

% 预分配矩阵：1000 行 (每次模拟)，24 列 (0-23 小时)
all_hourly_rain = zeros(N_sim, 24);

% 1. 跑 1000 次循环
for i = 1:N_sim
    % 【核心修改】：历时、峰值位置、开始时间 全部传入 []
    % 让函数内部自由发挥，模拟最真实的自然概率
    all_hourly_rain(i, :) = equations.downscale_rain_swat(R_day, [], [], []);
end

% 2. 数据统计分析
mean_rain = mean(all_hourly_rain, 1);        
p95_rain  = prctile(all_hourly_rain, 95, 1); 

fprintf('\n✓ 1000 次【完全随机】模拟完成！\n');
fprintf('设定单次降雨量   = %.2f mm\n', R_day);
fprintf('千次平均总降雨量 = %.2f mm (严格守恒)\n', sum(mean_rain));

% 3. 绘图展示
figure('Color', 'w', 'Position', [150, 150, 1000, 500]);

boxplot(all_hourly_rain, 'Colors', [0.6 0.6 0.6], 'Symbol', '.');
hold on;
plot(1:24, p95_rain, 'r--', 'LineWidth', 1.5, 'DisplayName', '95% 极端降雨上限');
plot(1:24, mean_rain, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', '千次平均降雨曲线');

title(sprintf('SWAT双指数降尺度: 1000次【完全解封】模拟 (日降雨 %d mm, 历时 2~12 h 随机)', R_day), 'FontSize', 14);
xlabel('时间 (小时, 0-23)', 'FontSize', 12);
ylabel('小时降雨量 (mm/h)', 'FontSize', 12);
xticklabels(0:23);
legend('Location', 'northeast', 'AutoUpdate', 'off'); 
grid on;