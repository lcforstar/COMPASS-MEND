function [ScopeParameters, xyt, canopy] = loadTimeSeries(ScopeParameters, ~, Soil, canopy, meteo, ~, xyt, ModelSettings, options)
%LOADTIMESERIES 将 MeteoData 逐小时数组映射到 SCOPE 所需的时间序列结构体
%
%   MeteoData 已在 run_stemmus_initialize 中完成：
%     load → to_hourly → compute_derived → compute_solar
%   本函数直接读取其属性，不再从文件加载。
%
%   输出不变：[ScopeParameters, xyt, canopy]

Nh   = meteo.n_hours;
DELT = ModelSettings.Time.Delt_t;
t0   = ModelSettings.Time.start_time;   % datetime，本地标准时

%% 1. 时间轴（SCOPE 约定：t = DOY-1 + hour/24，1月1日00:00 = 0）
timestamps = t0 + seconds((0:Nh-1)' * DELT);
doy_vec    = day(timestamps, 'dayofyear');
hr_vec     = hour(timestamps) + minute(timestamps)/60 + second(timestamps)/3600;
t_         = (doy_vec - 1) + hr_vec / 24;

xyt.t    = t_;
xyt.year = year(timestamps);

%% 2. 太阳天顶角
% MeteoData.compute_solar 已按经纬度逐小时计算，直接用
ScopeParameters.tts = meteo.tts;   % 太阳天顶角 (°), 0~90 已钳位

%% 3. 辐射
ScopeParameters.Rin = meteo.Rin;   % 入射短波辐射 (W/m²)
ScopeParameters.Rli = meteo.Rli;   % 入射长波辐射 (W/m²)

%% 4. 风速、气温、湿度、气压
ScopeParameters.u  = meteo.U;      % 风速 (m/s)
ScopeParameters.Ta = meteo.Ta;     % 气温 (°C)
ScopeParameters.Ts = meteo.Ts;     % 实测地表温度 (°C)，用于 Soil.Ts 初始化
ScopeParameters.ea = meteo.ea / 100; % Pa → hPa，SCOPE 单位
ScopeParameters.p  = meteo.Pg;      % 气压 (Pa)，SCOPE 单位

%% 5. 冠层结构（测量高度、冠层高度、LAI）
ScopeParameters.z   = ModelSettings.Env.Z_measure * ones(Nh, 1);  % 测风高度 (m)
% --- [改进] w/Cd/rwc 从统一新输入系统覆盖 xlsx / Override xlsx with new config ---
%   改进者 / Modified by:   朱寄子星 (Jizixing Zhu)
%   改进日期 / Modified date: 2026-04-18
ScopeParameters.leafwidth = ModelSettings.Veg.w_leaf_max * ones(Nh, 1);  % 叶宽 (m)
ScopeParameters.Cd        = ModelSettings.Veg.Cd  * ones(Nh, 1);                % 拖曳系数 [-]
ScopeParameters.rwc       = ModelSettings.Veg.rwc * ones(Nh, 1);                % 冠层内空气阻力 [s/m]
% Verhoef 方法专属参数（Method=verhoef 时用）
ScopeParameters.CR        = ModelSettings.Veg.CR     * ones(Nh, 1);             % 孤立树拖曳系数 [-]
ScopeParameters.CD1       = ModelSettings.Veg.CD1    * ones(Nh, 1);             % 拟合参数 [-]
ScopeParameters.Psicor    = ModelSettings.Veg.Psicor * ones(Nh, 1);             % 粗糙层修正 [-]
ScopeParameters.CSSOIL    = ModelSettings.Veg.CSSOIL * ones(Nh, 1);             % 土壤拖曳系数 [-]
ScopeParameters.LAI = ModelSettings.Veg.LAI;   % Nh×1，由 VegetationData.to_hourly 准备
ScopeParameters.hc  = ModelSettings.Veg.hc;    % Nh×1
ScopeParameters.SAI = ModelSettings.Veg.SAI;   % Nh×1
canopy.LAI         = ScopeParameters.LAI;
canopy.hc          = ScopeParameters.hc;
canopy.SAI         = ScopeParameters.SAI;
canopy.Method      = ModelSettings.Veg.Method;      % standard|table|optimal|verhoef
canopy.SurfaceType = ModelSettings.Veg.SurfaceType; % table 模式地表类型
if options.calc_zo && ~strcmpi(ModelSettings.Veg.LAI_mode, 'model')
    % constant/data 模式：LAI 已知，一次性计算全时间步 zo/d
    [ScopeParameters.zo, ScopeParameters.d] = equations.zo_and_d(Soil, canopy);
else
    ScopeParameters.zo = canopy.zo * ones(Nh, 1);
    ScopeParameters.d  = canopy.d  * ones(Nh, 1);
end

%% 6. CO₂ 浓度
% MeteoData 已处理：有实测列用实测，无则用 cfg.Env.Ca_default
ScopeParameters.Ca = meteo.Ca;     % CO₂浓度 (ppm)

%% 7. Leaf Vcmo, Tparam, m, Type, Rdparam, leafwidth（保留原注释供后续接入）
    % landcovers = unique(ModelSettings.landcoverClass, 'stable');
    % for timeIndex = 1:length(xyt.t)
    %     landcoverIndex = find(strcmp(ModelSettings.landcoverClass(timeIndex), landcovers));
    %     ScopeParameters.Vcmo(timeIndex, 1) = ScopeParameters.lcVcmo(landcoverIndex, 1);
    %     ScopeParameters.Tparam(timeIndex, :) = ScopeParameters.lcTparam(landcoverIndex, :);
    %     ScopeParameters.m(timeIndex, 1) = ScopeParameters.lcm(landcoverIndex, 1);
    %     ScopeParameters.Type(timeIndex, 1) = ScopeParameters.lcType(landcoverIndex, 1);
    %     ScopeParameters.Rdparam(timeIndex, 1) = ScopeParameters.lcRdparam(landcoverIndex, 1);
    %     ScopeParameters.leafwidth(timeIndex, 1) = ScopeParameters.lcleafwidth(landcoverIndex, 1);
    % end
    % % Remove the intermediate variables that were required to generate the timeseries:
    % ScopeParameters = rmfield(ScopeParameters, {'lcVcmo', 'lcTparam', 'lcm', 'lcType', 'lcRdparam', 'lcleafwidth'});
    %
    % if ~isempty(Cab_file)
    %     Cabtable = load(fullfile(ModelSettings.InputPath,'dataset dataYangling', Cab_file));
    %     ScopeParameters.Cab  = interp1(Cabtable(:, 1), Cabtable(:, 2), t_);
    % else
    %     ScopeParameters.Cab = leafbio.Cab * ones(size(t_));
    % end

%% 8. 土壤含水量（暂无时间序列，保持 ScopeParameters 原值）
% SMC 由 STEMMUS 土壤模块逐步更新，此处不覆盖
% if ~isempty(SMC_file)
%     ScopeParameters.SMC = load(...);
% end

end
