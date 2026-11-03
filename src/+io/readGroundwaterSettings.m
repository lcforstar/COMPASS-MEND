function GroundwaterSettings = readGroundwaterSettings(cfg)
% readGroundwaterSettings - 初始化地下水耦合设置
%                           Initialize groundwater coupling settings
%
% 函数定义 / Function Definition:
%   GroundwaterSettings = readGroundwaterSettings(cfg)
%
% 输入 / Inputs:
%   cfg - 全局模型配置对象 / Global model settings (legacy ModelSettings)
%         含 cfg.Grid.NN, cfg.Grid.NL, cfg.Grid.DeltZ, cfg.Grid.Tot_Depth
%
% 输出 / Outputs:
%   GroundwaterSettings - 地下水耦合设置结构体 / Groundwater coupling struct
%     .GroundwaterCoupling  - 耦合开关 0=关闭,1=开启 (默认0，通过BMI更新)
%     .headBotmLayer        - 底层水头初值 (cm)，由MODFLOW通过BMI传入
%     .indexBotmLayer       - 含当前底层水头的层索引，由MODFLOW通过BMI传入
%     .soilLayerThickness   - 累积层厚 (NN×1, cm)，从地表向下累加
%     .topLevel             - 含水层顶部高程 (cm)，由MODFLOW通过BMI传入
%     .tempBotm             - 底层地下水温度 (°C)，由MODFLOW通过BMI传入 (NaN=缺测)
%
% 参考 / References:
%   Mostafa et al. (2022) GMD preprint — MODFLOW-STEMMUS BMI coupling
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-13

NN = cfg.Grid.NN;
NL = cfg.Grid.NL;
TDeltZ = flip(cfg.Grid.DeltZ);   % DeltZ 从底到顶，flip后从顶到底

% --- [1] 耦合开关与初始水头 / Coupling switch and initial head ---
% Default off; set GroundwaterCoupling=1 via BMI to enable MODFLOW coupling
GroundwaterSettings.GroundwaterCoupling = 0;
GroundwaterSettings.headBotmLayer   = 1750.0;  % groundwater head at bottom layer [cm]
GroundwaterSettings.topLevel        = 2000.0;  % elevation of top surface aquifer layer [cm]
GroundwaterSettings.tempBotm        = NaN;     % groundwater temperature at bottom layer [degC]
GroundwaterSettings.indexBotmLayer  = 40;

% --- [2] 累积层厚 / Cumulative soil layer thickness ---
soilLayerThickness = zeros(NN, 1);
for ML = 2:NL
    soilLayerThickness(ML) = soilLayerThickness(ML - 1) + TDeltZ(ML - 1);
end
soilLayerThickness(NN) = cfg.Grid.Tot_Depth;
GroundwaterSettings.soilLayerThickness = soilLayerThickness;

end
