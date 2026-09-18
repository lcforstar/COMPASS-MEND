function sunTrajectory()
% 定义参数
R = 1;          % 半径
beta_deg = 70;  % 角度参数 (度)
theta_deg = 40; % 角度参数 (度)

% 将角度转换为弧度
beta = beta_deg * pi / 180;
theta = theta_deg * pi / 180;

d = R * cos(beta); % 计算d的值

% 创建φ值数组 (0 到 2π)
phi = linspace(0, 2 * pi, 100);

% 计算轨迹的x, y, z坐标 (核心投影公式)
x = R * sin(beta) * sin(phi);
y = R * sin(beta) * cos(phi) * sin(theta) + d * cos(theta);
z = R * sin(beta) * cos(phi) * cos(theta) - d * sin(theta);

% 创建3D图形
figure;
ax = axes('Parent', gcf);
hold(ax, 'on');
grid(ax, 'on');
view(ax, 3); % 设置3D视角

% 绘制原点
scatter3(ax, 0, 0, 0, 'green', 'filled', 'DisplayName', 'Origin');

% 绘制基础参考平面（例如：地平圈）
circle_x = R * cos(phi);
circle_y = R * sin(phi);
plot3(ax, circle_x, circle_y, zeros(size(phi)), 'Color', '#87CEFA', 'LineWidth', 1);

% 绘制轨迹，颜色根据z的值确定 (分段绘制以实现颜色区分)
for i = 1:(length(z)-1)
    if z(i) > 0
        plot3(ax, x(i:i+1), y(i:i+1), z(i:i+1), 'Color', 'red', 'LineWidth', 2);
    else
        plot3(ax, x(i:i+1), y(i:i+1), z(i:i+1), 'Color', 'blue', 'LineWidth', 2);
    end
end

% 创建并显示图例
% 创建虚拟绘图对象用于图例
h1 = plot3(ax, NaN, NaN, NaN, 'r-', 'LineWidth', 2, 'DisplayName', 'Day (z>0)');
h2 = plot3(ax, NaN, NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName', 'Night (z<0)');
legend(ax, [h1, h2], 'Location', 'best');

% 设置坐标轴标签
xlabel(ax, 'X');
ylabel(ax, 'Y');
zlabel(ax, 'Z');
title(ax, 'Sun Trajectory Model with Day/Night Distinction');

% 保持图形显示
hold(ax, 'off');
end