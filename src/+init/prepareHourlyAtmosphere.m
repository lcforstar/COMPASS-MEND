function AtmosphereHourly = prepareHourlyAtmosphere(ForcingData, ModelSettings)
    % PREPAREHOURLVATMOSPHERE 将气象强制数据转换为小时尺度
    %
    % 语法：
    %   AtmosphereHourly = prepareHourlyAtmosphere(ForcingData, ModelSettings)
    %
    % 输入参数：
    %   ForcingData      - 从 loadForcingData 加载的强制气象数据结构体
    %   ModelSettings    - 模型设置结构体（包含 DELT、lTimeScale 等）
    %
    % 输出参数：
    %   AtmosphereHourly - 小时尺度气象参数结构体，包含以下字段：
    %      基本气温变量 (°C):
    %        Ta_hours      - 平均气温
    %        Tmax_hours    - 最高气温
    %        Tmin_hours    - 最低气温
    %        Tdew_hours    - 露点温度
    %        Ts_hours      - 土壤温度
    %        Tb_hours      - 其他温度参数
    %      相对湿度变量 (%):
    %        RH_hours      - 平均相对湿度
    %        Rhmax_hours   - 最高相对湿度
    %        Rhmin_hours   - 最低相对湿度
    %      风速变量 (m/s):
    %        U0_hours      - 纵向风速
    %        Us_hours      - 径向风速
    %        U_hours       - 合成风速
    %      辐射变量 (W/m²):
    %        Rn_hours      - 净辐射 (用户自定义转换)
    %        Rs_hours      - 短波辐射 (用户自定义转换)
    %      降雨和其他:
    %        Precip_hours  - 降雨量 (mm)
    %        N_hours       - 理论日照 (h)
    %        n_hours       - 实际日照 (h)
    %        LAI_hours     - 叶面积指数
    %        h_plant_hours - 株高 (m)
    %        h_SUR_hours   - 表面基质势 (cm H₂O)
    %        Pg_hours      - 气压 (Pa)
    %      物理参数:
    %        P_Va_hours    - 水蒸气大气压 (kPa)
    %        RHOV_A_hours  - 水蒸气密度 (g/cm³)
    %        TopPg_hours   - 地表气压 (Pa)

    DELT = ModelSettings.Time.Delt_t;               % 时间步长 (s)
    Dur_tot = ModelSettings.Time.Nmsrmn;            % 总时间步数
    Z = ModelSettings.Elevation;                    % 高程 (m)
    lTimeScale = ModelSettings.Time.lTimeScale;     % 时间尺度: 0=hourly, 1=daily

    fprintf('[Prepare Hourly Atmosphere Data]\n');
    if lTimeScale == 0
        fprintf('  Time Scale: hourly (no conversion)\n');
    elseif lTimeScale == 1
        fprintf('  Time Scale: daily (will convert to hourly)\n');
    end
    fprintf('  Time Step: %.0f sec\n', DELT);

    % =====================================================================
    % 初始化输出结构体
    % =====================================================================
    AtmosphereHourly = struct();

    % 直接赋值小时尺度数据
    if lTimeScale == 0  % 0 = hourly (无需转换)
        % 输入已经是小时尺度，无需转换
        fprintf('  [OK] Input data is hourly, no conversion needed.\n\n');

        % 温度变量 (°C)
        AtmosphereHourly.Ta_hours = ForcingData.Ta_msr;
        AtmosphereHourly.Tdew_hours = ForcingData.Tdew_msr;
        AtmosphereHourly.Ts_hours = ForcingData.Ts_msr;
        AtmosphereHourly.Tb_hours = ForcingData.Tb_msr;

        % 相对湿度变量 (%)
        AtmosphereHourly.RH_hours = ForcingData.RH_msr;

        % 风速变量 (m/s)
        AtmosphereHourly.Uu_hours = ForcingData.Uu_msr;
        AtmosphereHourly.Uv_hours = ForcingData.Uv_msr;

        % 风速处理：优先使用U_msr（如果存在且有效），否则从Uu和Uv合成
        if isfield(ForcingData, 'U_msr') && ~all(isnan(ForcingData.U_msr)) && ~all(ForcingData.U_msr == 0)
            AtmosphereHourly.U_hours = ForcingData.U_msr;  % 直接使用合成风速
        else
            AtmosphereHourly.U_hours = sqrt(ForcingData.Uu_msr.^2 + ForcingData.Uv_msr.^2);  % 从分量合成
        end

        % 降雨和日照
        AtmosphereHourly.Precip_hours = ForcingData.Precip_msr;
        AtmosphereHourly.n_hours = ForcingData.n_msr;

        % 植被和土壤参数
        AtmosphereHourly.LAI_hours = ForcingData.LAI_msr;
        AtmosphereHourly.h_plant_hours = ForcingData.h_plant_msr;
        AtmosphereHourly.h_SUR_hours = ForcingData.h_SUR_msr;

        % 气压 (Pa)
        AtmosphereHourly.Pg_hours = ForcingData.Pg_msr;

        % 辐射 (用户自定义)
        AtmosphereHourly.Rn_hours = ForcingData.Rn_msr;
        % AtmosphereHourly.Rs_hours = ... (用户需要自定义小时尺度转换)

    elseif lTimeScale == 1  % 1 = daily (需要转换)
        % 输入是天尺度，需要转换为小时尺度
        fprintf('  [Converting] Daily data to hourly...\n');

        % 数据长度和时间信息
        n_data = length(ForcingData.Ta_msr);  % 天尺度数据点数
        n_hours = Dur_tot;  % 模型需要的小时尺度时间步数

        % 初始化小时尺度数组
        Ta_hours = zeros(n_hours, 1);
        Tdew_hours = zeros(n_hours, 1);
        Ts_hours = zeros(n_hours, 1);
        Tb_hours = zeros(n_hours, 1);

        RH_hours = zeros(n_hours, 1);
       
        Uu_hours = zeros(n_hours, 1);
        Uv_hours = zeros(n_hours, 1);
        U_hours = zeros(n_hours, 1);

        Precip_hours = zeros(n_hours, 1);     
        n_hours = zeros(n_hours, 1);

        LAI_hours = zeros(n_hours, 1);
        h_plant_hours = zeros(n_hours, 1);
        
        h_SUR_hours = zeros(n_hours, 1);

        Pg_hours = zeros(n_hours, 1);

        % 天尺度到小时尺度的转换
        for kt = 1:n_hours
            % 计算当前时刻属于第几天
            day_idx = floor((kt - 1) * DELT / 86400) + 1;
            if day_idx > n_data
                day_idx = n_data;
            end

            % 计算当前时刻在当天的小时
            hour_in_day = mod((kt - 1) * DELT, 86400) / 3600;

            % ========== 温度日较差法转换 ==========
            % Ta(t) = (Tmax + Tmin)/2 + (Tmax - Tmin)/2 * cos(2*pi*(hour-15)/24)
            Tmax_val = ForcingData.Tmax_msr(day_idx);
            Tmin_val = ForcingData.Tmin_msr(day_idx);

            Ta_hours(kt) = (Tmax_val + Tmin_val) / 2 + (Tmax_val - Tmin_val) / 2 * ...
                cos(2*pi*(hour_in_day - 15) / 24);
            Tmax_hours(kt) = Tmax_val;
            Tmin_hours(kt) = Tmin_val;

            if isfield(ForcingData, 'Tdew_msr')
                Tdew_hours(kt) = ForcingData.Tdew_msr(day_idx);
            end

            % 土壤温度
            if isfield(ForcingData, 'Ts_msr')
                Ts_hours(kt) = ForcingData.Ts_msr(day_idx);
            end

            if isfield(ForcingData, 'Tb_msr')
                Tb_hours(kt) = ForcingData.Tb_msr(day_idx);
            end

            % ========== 相对湿度日较差法转换 ==========
            % RH(t) = (RHmax + RHmin)/2 + (RHmax - RHmin)/2 * cos(2*pi*(hour-6)/24)
            Rhmax_val = ForcingData.Rhmax_msr(day_idx);
            Rhmin_val = ForcingData.Rhmin_msr(day_idx);

            RH_hours(kt) = (Rhmax_val + Rhmin_val) / 2 + (Rhmax_val - Rhmin_val) / 2 * ...
                cos(2*pi*(hour_in_day - 6) / 24);
            Rhmax_hours(kt) = Rhmax_val;
            Rhmin_hours(kt) = Rhmin_val;

            % ========== 风速日较差法转换 ==========
            Uu_hours(kt) = ForcingData.Uu_msr(day_idx);
            Uv_hours(kt) = ForcingData.Uv_msr(day_idx);

            % 风速处理：优先使用U_msr（如果存在且有效），否则从Uu和Uv合成并进行日较差转换
            if isfield(ForcingData, 'U_msr') && ~isnan(ForcingData.U_msr(day_idx)) && ForcingData.U_msr(day_idx) ~= 0
                % 直接使用合成风速（假设U_msr已是小时尺度或日平均值）
                U_hours(kt) = ForcingData.U_msr(day_idx);
            else
                % 从分量合成，并应用日较差法
                Um = sqrt(U0_hours(kt)^2 + Us_hours(kt)^2) / sqrt(2);

                % 风速日较差：设最大风速在下午15点，最小风速在深夜
                Ur = 2.5;  % 默认比值，可根据风速调整
                Umax = 2 * Ur * Um / (1 + Ur);
                Umin = 2 * Um / (1 + Ur);
                U_hours(kt) = Um + (Umax - Umin) * 0.5 * cos(2*pi*(hour_in_day - 15) / 24);
            end

            % 降雨：均匀分配到该天的每个小时
            hours_per_day = 24;
            Precip_hours(kt) = ForcingData.Precip_msr(day_idx) / hours_per_day;

            % 日照
            if isfield(ForcingData, 'N_msr')
                N_hours(kt) = ForcingData.N_msr(day_idx);
            end
            if isfield(ForcingData, 'n_msr')
                n_hours_data(kt) = ForcingData.n_msr(day_idx);
            end

            % 植被参数
            if isfield(ForcingData, 'LAI_msr')
                LAI_hours(kt) = ForcingData.LAI_msr(day_idx);
            end
            if isfield(ForcingData, 'h_plant_msr')
                h_plant_hours(kt) = ForcingData.h_plant_msr(day_idx);
            end
            if isfield(ForcingData, 'h_SUR_msr')
                h_SUR_hours(kt) = ForcingData.h_SUR_msr(day_idx);
            end

            % 气压
            if isfield(ForcingData, 'Pg_msr')
                Pg_hours(kt) = ForcingData.Pg_msr(day_idx);
            end
        end

        % 赋值到输出结构体
        AtmosphereHourly.Ta_hours = Ta_hours;
        AtmosphereHourly.Tmax_hours = Tmax_hours;
        AtmosphereHourly.Tmin_hours = Tmin_hours;
        AtmosphereHourly.Tdew_hours = Tdew_hours;
        AtmosphereHourly.Ts_hours = Ts_hours;
        AtmosphereHourly.Tb_hours = Tb_hours;

        AtmosphereHourly.RH_hours = RH_hours;
        AtmosphereHourly.Rhmax_hours = Rhmax_hours;
        AtmosphereHourly.Rhmin_hours = Rhmin_hours;

        AtmosphereHourly.U0_hours = U0_hours;
        AtmosphereHourly.Us_hours = Us_hours;
        AtmosphereHourly.U_hours = U_hours;

        AtmosphereHourly.Precip_hours = Precip_hours;
        AtmosphereHourly.N_hours = N_hours;
        AtmosphereHourly.n_hours = n_hours_data;

        AtmosphereHourly.LAI_hours = LAI_hours;
        AtmosphereHourly.h_plant_hours = h_plant_hours;
        AtmosphereHourly.h_SUR_hours = h_SUR_hours;

        AtmosphereHourly.Pg_hours = Pg_hours;

        fprintf('  [OK] Conversion complete: %d days -> %d hours\n', n_data, n_hours);

    else
        error('lTimeScale must be 0 (hourly) or 1 (daily), current value: %d', lTimeScale);
    end

    % =====================================================================
    % 计算衍生物理参数（小时尺度）
    % =====================================================================
    fprintf('  Computing derived physical parameters...\n');

    % 物理常数
    Rv = 461.5;  % 水蒸气气体常数 (J/(kg·K)) = cm²/(s²·K)

    % 水蒸气大气压 (kPa)
    P_Va_hours = 0.611 * exp(17.27 * AtmosphereHourly.Ta_hours ./ (AtmosphereHourly.Ta_hours + 237.3)) .* ...
        (AtmosphereHourly.RH_hours / 100);

    % 水蒸气密度 (g/cm³)
    RHOV_A_hours = P_Va_hours * 1e4 ./ (Rv * (AtmosphereHourly.Ta_hours + 273.15));

    % 地表气压 (Pa) - 使用标准气压衰减公式
    TopPg_hours = 1013250 * exp((-1) * 9.80665 * 0.0289644 * Z ./ (8.31432 .* (AtmosphereHourly.Ts_hours + 273.15)));

    % 赋值
    AtmosphereHourly.P_Va_hours = P_Va_hours;
    AtmosphereHourly.RHOV_A_hours = RHOV_A_hours;
    AtmosphereHourly.TopPg_hours = TopPg_hours;

    % =====================================================================
    % 辐射转换（用户自定义）
    % =====================================================================
    fprintf('  NOTE: Radiation data conversion is user-defined.\n');
    fprintf('        (Rn_msr and Rs_msr are provided as daily values)\n');
    fprintf('        Please add hourly radiation conversion logic here.\n\n');

    % 示意：用户应该在这里添加辐射转换逻辑
    % AtmosphereHourly.Rn_hours = ... 用户自定义
    % AtmosphereHourly.Rs_hours = ... 用户自定义

    fprintf('[OK] Atmosphere data preparation complete.\n\n');

end
