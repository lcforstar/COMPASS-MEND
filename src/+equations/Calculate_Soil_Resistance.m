function r_ss = Calculate_Soil_Resistance(Theta_surf, Theta_r, Theta_s,Theta_f)
% Calculate_Soil_Resistance 计算土壤表面阻力
%
% [改进 2026-04-18 朱寄子星] 从 +CalculateEvapotranspiration/ 迁移至 +equations/
%   按物理方程归类，与 zo_and_d.m 等阻力类方程并置。
%   调用改为 equations.Calculate_Soil_Resistance(...)
%
% CLM4.5-style improvement:
%   Jizixing Zhu (朱寄子星), 18 Apr 2026.
%
% 方法: S4 - Mahfouf and Noilhan (1991) / Sellers-style soil evaporation resistance
%   ebal_clm45 每个时间步用表层土壤含水量 Theta_surf 动态更新 rss，
%   用于替代原始固定 soil.rss。该公式使干土的水汽扩散阻力快速增大，
%   湿土的蒸发阻力降低，从而把土壤水分状态反馈到潜热通量。
%
% Method document:
%   docs/clm45_aerodynamic_resistance_method.md
 is=1;
if is ==1 
% 1. 获取输入限制
% 防止数值越界
Theta_surf = max(0, Theta_surf);
Theta_s    = max(1e-3, Theta_s);

% 2. 计算相对饱和度 (Ratio)
% Mahfouf 公式使用的是 theta / theta_sat
ratio = Theta_surf / Theta_s;
ratio = max(0, min(1.0, ratio));

% 3. 计算阻力 (Mahfouf and Noilhan, 1991)
% 公式: r_ss = 38113 * exp(-13.515 * ratio)
% 物理含义: 
% - 极干 (ratio=0): r_ss = 38113 (极大)
% - 极湿 (ratio=1): r_ss = 38113 * 1.3e-6 = 0.05 (极小)
% - 半干 (ratio=0.3): r_ss ≈ 660 (适中，正好修正你的误差)

coeff_A = 38113; % 或者 3.8e4
coeff_B = 12.0;

r_ss = coeff_A * exp(-coeff_B * ratio);

% 4. 必要的边界限制 (Clamping)
% 测试记录 2026-05-10: 尝试 rss_max 2000→5000 想抑制冻土LE，结果 lEs 反升
%   (12.1→15.0 W/m²)。ebal 非线性反馈让 Ts 重新平衡，rss 单点修改改变不了 LE
%   分配。维持 2000。
r_ss_min = 50;   % 最小阻力
r_ss_max = 2000; % 最大阻力

r_ss = max(r_ss_min, min(r_ss_max, r_ss));
else   
% Calculate_Soil_Resistance 计算土壤表面阻力
% 方法: Sellers et al. (1992) / SCOPE 模型策略
% 逻辑: 使用田持水量 (Field Capacity) 作为参考点，指数衰减
%
% 输入:
%   Theta_surf : 表层土壤含水量 (cm3/cm3) -> 对应 SP.SMC
%   Theta_r    : 土壤残余含水量 (cm3/cm3) -> 对应 SP.ResidualMC
%   Theta_f    : 田持水量 (cm3/cm3)       -> 对应 SP.fieldMC
%
% 输出:
%   r_ss       : 土壤表面阻力 [s/m]

    % 1. 定义系数 (来自你提供的代码片段)
    aa = 3.8; 
    
    % 截距项 (aa + 4.1 = 7.9)
    % 物理含义: 当水分接近残余值时，exp(7.9) ≈ 2700 s/m (最大阻力)
    intercept = aa + 4.1; 
    
    % 2. 准备分母 (田持水量 - 残余含水量)
    denom = Theta_f - Theta_r;
    if denom < 1e-6
        denom = 1e-6; % 防止除以0
    end
    
    % 3. 计算相对含水量权重 (Weight)
    % 公式: (SMC - Res) / (Field - Res)
    weight = (Theta_surf - Theta_r) ./ denom;
    
    % 4. 边界限制 (很重要)
    % 物理上，如果 Theta_surf > Theta_f (过饱和/积水)，weight > 1
    % 此时 exp(7.9 - 3.8 * >1) 会变得很小 (<60)，这是合理的。
    % 但是，如果 Theta_surf < Theta_r，weight < 0，exp 会爆炸。
    % 所以至少要限制下限为 0 (即 Theta_surf 不低于 Theta_r)
    weight = max(0.0, weight);
    
    % 5. 计算阻力
    % rss = exp((aa + 4.1) - aa * weight)
    r_ss = exp(intercept - aa .* weight);
    
    % 6. 安全边界限制 (Clamping)
    % 理论值参考:
    % - At Field Cap (w=1): exp(7.9 - 3.8) = exp(4.1) ≈ 60 s/m (湿润)
    % - At Residual (w=0):  exp(7.9) ≈ 2700 s/m (干燥)
    
    r_ss_min = 50;   % 最小阻力
    r_ss_max = 3000; % 最大阻力 (防止数值溢出，Sellers模型也就是2700左右)
    
    r_ss = max(r_ss_min, min(r_ss_max, r_ss));

end
end
