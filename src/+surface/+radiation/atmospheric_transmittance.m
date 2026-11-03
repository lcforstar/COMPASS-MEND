function [T, tau_total, tau_components] = atmospheric_transmittance(lambda, h, theta, p, MR0, X, w, varargin)
%ATMOSPHERIC_TRANSMITTANCE 计算太阳辐射的大气透射率
%   基于Elterman（1968,1970）和Gates and Harrop（1963）模型
%   整合了瑞利散射、气溶胶消光、臭氧吸收和水汽吸收的贡献
%
%   输入参数:
%     lambda: 波长（微米），范围0.27-4.0 μm [μm]
%     h: 海拔高度（千米） [km]
%     theta: 太阳天顶角（度） [度]
%     p: 海平面气压（毫巴），默认1013 mb [mb]
%     MR0: 海平面气象视程（千米），默认25 km [km]
%     X: 总臭氧含量（厘米NTP），默认0.34 cm [cm NTP]
%     w: 总可降水水汽（厘米），默认1.0 cm [cm]
%     varargin: 可选参数，用于水汽剖面比率A4(h)或其他设置
%
%   输出参数:
%     T: 大气透射率（0-1） [无量纲]
%     tau_total: 总光学厚度 [无量纲]
%     tau_components: 各成分光学厚度结构体 [无量纲]
%
%   参考文献:
%     Elterman, L. (1968). UV, visible, and IR attenuation for altitudes to 50 km
%     Gates, D. M., & Harrop, W. J. (1963). Infrared transmission of the atmosphere
%     Robinson, N. (1966). Solar Radiation

% ======================== 参数验证和默认值 ========================
if nargin < 7
    w = 1.0; % 默认可降水水汽1 cm
end
if nargin < 6
    X = 0.34; % 默认臭氧含量0.34 cm NTP
end
if nargin < 5
    MR0 = 25; % 默认气象视程25 km
end
if nargin < 4
    p = 1013; % 默认海平面气压1013 mb
end

% 验证输入参数
validateattributes(lambda, {'numeric'}, {'scalar', '>=', 0.27, '<=', 4.0});
validateattributes(h, {'numeric'}, {'scalar', '>=', 0, '<=', 50});
validateattributes(theta, {'numeric'}, {'scalar', '>=', 0, '<=', 90});
validateattributes(p, {'numeric'}, {'scalar', '>', 0});
validateattributes(MR0, {'numeric'}, {'scalar', '>', 0});
validateattributes(X, {'numeric'}, {'scalar', '>', 0});
validateattributes(w, {'numeric'}, {'scalar', '>', 0});

% ======================== 大气质量计算 ========================
m = calculate_air_mass(theta); % 计算大气质量

% ======================== 海平面消光系数获取 ========================
% 从Elterman和Gates & Harrop表获取消光系数（需要实际数据）
[tau_R0, tau_A0, tau_O0, tau_W0] = get_sea_level_extinction(lambda);

% ======================== 海拔减少因子获取 ========================
% 从表3获取A1, A2, A3（瑞利、气溶胶、臭氧）
A1 = get_altitude_factor(h, 1); % 瑞利减少因子
A2 = get_altitude_factor(h, 2); % 气溶胶减少因子
A3 = get_altitude_factor(h, 3); % 臭氧减少因子

% 水汽减少因子A4（需要特殊处理，见图片4说明）
A4 = get_water_vapor_factor(h, w, varargin{:});

% ======================== 光学厚度计算（公式14） ========================
% 瑞利散射光学厚度（公式14第一部分）
tau_R = m * tau_R0 * (p / 1013) * A1;

% 气溶胶消光光学厚度（公式14第二部分）
tau_A = m * tau_A0 * (25 / MR0) * A2;

% 臭氧吸收光学厚度（公式14第三部分）
tau_O = m * tau_O0 * (X / 0.34) * A3;

% 水汽吸收光学厚度（公式13和14第四部分）
tau_W = tau_W0 * sqrt(w * m * A4);

% 总光学厚度
tau_total = tau_R + tau_A + tau_O + tau_W;

% ======================== 透射率计算 ========================
T = exp(-tau_total);

% ======================== 输出组件光学厚度 ========================
tau_components.tau_R = tau_R;
tau_components.tau_A = tau_A;
tau_components.tau_O = tau_O;
tau_components.tau_W = tau_W;

% ======================== 结果显示 ========================
displayTransmittanceResults(lambda, h, theta, T, tau_total, tau_components);

end

%% ======================== 辅助函数 ========================
function m = calculate_air_mass(theta)
%CALCULATE_AIR_MASS 计算大气质量
%   基于太阳天顶角计算相对大气质量
%
%   输入参数:
%     theta: 太阳天顶角（度） [度]
%
%   输出参数:
%     m: 大气质量（相对值） [无量纲]

    theta_rad = theta * pi / 180; % 转换为弧度
    
    % 简单的大气质量公式（当theta < 75°时精度较好）
    if theta < 75
        m = 1 / cos(theta_rad);
    else
        % 对于大天顶角使用更精确的公式
        m = 1 / (cos(theta_rad) + 0.50572 * (96.07995 - theta)^(-1.6364));
    end
end

function [tau_R0, tau_A0, tau_O0, tau_W0] = get_sea_level_extinction(lambda)
%GET_SEA_LEVEL_EXTINCTION 获取海平面消光系数
%   从Elterman和Gates & Harrop表获取消光系数
%   注意：这里需要实际的表格数据，目前使用近似值
%
%   输入参数:
%     lambda: 波长（微米） [μm]
%
%   输出参数:
%     tau_R0: 瑞利散射海平面消光系数 [无量纲]
%     tau_A0: 气溶胶海平面消光系数 [无量纲]
%     tau_O0: 臭氧海平面消光系数 [无量纲]
%     tau_W0: 水汽海平面消光系数（1 cm可降水） [无量纲]

    % 这里应该是从表格中插值的代码
    % 由于缺乏完整数据表，这里使用近似公式
    
    % 瑞利散射（与λ^(-4)成正比）
    tau_R0 = 0.0086 * (0.55 / lambda)^4;
    
    % 气溶胶（与λ^(-1.3)成正比，Angstrom定律）
    tau_A0 = 0.1 * (0.55 / lambda)^1.3;
    
    % 臭氧（在紫外和可见光波段有吸收带）
    % 简化处理：在哈特莱-哈金斯带（0.2-0.3μm）和查普uis带（0.5-0.7μm）
    if lambda >= 0.2 && lambda <= 0.3
        tau_O0 = 0.05 * exp(-10 * (lambda - 0.25)^2);
    elseif lambda >= 0.5 && lambda <= 0.7
        tau_O0 = 0.02 * exp(-20 * (lambda - 0.6)^2);
    else
        tau_O0 = 0.001;
    end
    
    % 水汽（在红外波段有吸收带）
    % 主要吸收带：0.94, 1.1, 1.38, 1.87, 2.7, 3.2 μm
    if lambda >= 0.8 && lambda <= 1.0
        tau_W0 = 0.1 * exp(-50 * (lambda - 0.94)^2);
    elseif lambda >= 1.0 && lambda <= 1.2
        tau_W0 = 0.08 * exp(-40 * (lambda - 1.1)^2);
    elseif lambda >= 1.3 && lambda <= 1.5
        tau_W0 = 0.15 * exp(-60 * (lambda - 1.38)^2);
    elseif lambda >= 1.7 && lambda <= 2.0
        tau_W0 = 0.2 * exp(-50 * (lambda - 1.87)^2);
    elseif lambda >= 2.5 && lambda <= 3.0
        tau_W0 = 0.3 * exp(-40 * (lambda - 2.7)^2);
    elseif lambda >= 3.0 && lambda <= 3.5
        tau_W0 = 0.25 * exp(-50 * (lambda - 3.2)^2);
    else
        tau_W0 = 0.01;
    end
    
    % 注意：这些是近似值，实际应用需要完整的Elterman和Gates & Harrop数据表
end

function A = get_altitude_factor(h, component_type)
%GET_ALTITUDE_FACTOR 获取海拔减少因子（表3数据）
%   基于表3的插值计算
%
%   输入参数:
%     h: 海拔高度（km） [km]
%     component_type: 成分类型（1=瑞利, 2=气溶胶, 3=臭氧）
%
%   输出参数:
%     A: 海拔减少因子 [无量纲]

    % 表3数据（海拔km vs A1, A2, A3）
    altitude_data = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    
    switch component_type
        case 1 % 瑞利散射（A1）
            A_data = [1.000, 0.887, 0.784, 0.692, 0.609, 0.534, 0.467, 0.406, 0.353, 0.304, 0.262];
        case 2 % 气溶胶（A2）
            A_data = [1.000, 0.546, 0.353, 0.260, 0.223, 0.199, 0.181, 0.168, 0.155, 0.142, 0.129];
        case 3 % 臭氧（A3）
            A_data = [1.000, 0.999, 0.981, 0.973, 0.966, 0.959, 0.952, 0.946, 0.939, 0.932, 0.923];
        otherwise
            error('无效的成分类型。必须是1（瑞利）、2（气溶胶）或3（臭氧）');
    end
    
    % 线性插值
    if h <= 0
        A = A_data(1);
    elseif h >= 10
        A = A_data(end);
    else
        A = interp1(altitude_data, A_data, h, 'linear');
    end
end

function A4 = get_water_vapor_factor(h, w, varargin)
%GET_WATER_VAPOR_FACTOR 获取水汽海拔减少因子A4(h)
%   水汽垂直分布变化大，需要特殊处理（图片4说明）
%
%   输入参数:
%     h: 海拔高度（km） [km]
%     w: 总可降水水汽（cm） [cm]
%     varargin: 可选参数，可提供实际A4值或分布模型
%
%   输出参数:
%     A4: 水汽海拔减少因子 [无量纲]

    % 如果有用户提供的A4值，使用它
    if nargin > 2 && ~isempty(varargin)
        if isnumeric(varargin{1})
            A4 = varargin{1};
            return;
        end
    end
    
    % 默认方法：使用指数衰减模型近似水汽分布
    % 水汽随高度近似指数减少：A4(h) = exp(-h / H_w)
    % H_w = 尺度高度，典型值2-3 km
    H_w = 2.5; % km，水汽尺度高度
    
    A4 = exp(-h / H_w);
    
    % 注意：这是近似值。实际应用中应使用无线电探空仪数据（图片4说明）
end

function displayTransmittanceResults(lambda, h, theta, T, tau_total, tau_components)
%DISPLAYTRANSMITTANCERESULTS 显示透射率计算结果
%   格式化显示大气透射率和各成分光学厚度
%
%   输入参数:
%     lambda: 波长（μm） [μm]
%     h: 海拔高度（km） [km]
%     theta: 太阳天顶角（度） [度]
%     T: 透射率 [无量纲]
%     tau_total: 总光学厚度 [无量纲]
%     tau_components: 各成分光学厚度结构体 [无量纲]

    fprintf('\n=== 大气透射率计算结果 ===\n');
    fprintf('输入参数:\n');
    fprintf('   波长: %.2f μm\n', lambda);
    fprintf('   海拔高度: %.1f km\n', h);
    fprintf('   太阳天顶角: %.1f°\n', theta);
    
    fprintf('\n光学厚度分析:\n');
    fprintf('   总光学厚度: %.4f\n', tau_total);
    fprintf('   瑞利散射光学厚度: %.4f (占%.1f%%)\n', ...
        tau_components.tau_R, 100*tau_components.tau_R/tau_total);
    fprintf('   气溶胶光学厚度: %.4f (占%.1f%%)\n', ...
        tau_components.tau_A, 100*tau_components.tau_A/tau_total);
    fprintf('   臭氧吸收光学厚度: %.4f (占%.1f%%)\n', ...
        tau_components.tau_O, 100*tau_components.tau_O/tau_total);
    fprintf('   水汽吸收光学厚度: %.4f (占%.1f%%)\n', ...
        tau_components.tau_W, 100*tau_components.tau_W/tau_total);
    
    fprintf('\n透射率结果:\n');
    fprintf('   大气透射率: %.4f (%.2f%%)\n', T, T*100);
    
    % 根据透射率分类
    if T >= 0.8
        fprintf('   大气条件: 非常清澈\n');
    elseif T >= 0.6
        fprintf('   大气条件: 清澈\n');
    elseif T >= 0.4
        fprintf('   大气条件: 中等\n');
    elseif T >= 0.2
        fprintf('   大气条件: 浑浊\n');
    else
        fprintf('   大气条件: 非常浑浊/高吸收\n');
    end
    
    fprintf('\n模型说明:\n');
    fprintf('   基于Elterman (1968,1970) 和 Gates & Harrop (1963) 模型\n');
    fprintf('   使用了表3的海拔减少因子和公式(12)-(14)\n');
end