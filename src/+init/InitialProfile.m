classdef InitialProfile
    %INITIALPROFILE 初始剖面模块 - 主协调器
    %
    % 本模块负责：
    %   1. 从 Excel 文件读取初始剖面数据
    %   2. 执行多种插值计算（温度、含水量、气压、氮浓度）
    %   3. 物理验证与自动修正
    %   4. 坐标系统转换（表层→底层 变为 底层→表层）
    %   5. 初始化模型状态变量
    %
    % 关键坐标系统说明：
    %   - Excel 数据顺序：[表层 0.0 cm → 底层 TotDepth]
    %   - 插值过程：保持表层→底层方向
    %   - 最终存储：倒置为 [底层 → 表层]（符合STEMMUS标准，索引0在底部）
    %
    % 使用示例：
    %   init_prof = init.InitialProfile(ModelSettings, VanGenuchten);
    %   state = init_prof.load_and_interpolate('InitialConditions.xlsx');

    properties
        ModelSettings   % 模型设置
        SoilProperties  % 土壤性质（包含van Genuchten参数）
        verbose = true  % 是否输出详细信息
    end

    methods
        function obj = InitialProfile(ModelSettings, SoilProperties)
            %INITIALPROFILE 构造函数
            %
            % 输入参数：
            %   ModelSettings  - 模型设置对象
            %   SoilProperties - 土壤性质对象（可选，含van Genuchten参数）

            if nargin < 1
                error('需要至少一个输入参数：ModelSettings');
            end

            obj.ModelSettings = ModelSettings;
            obj.SoilProperties = SoilProperties;
        end

        function initial_state = load_and_interpolate(obj, filename, sheetname)
            %LOAD_AND_INTERPOLATE 主流程：从Excel读取→插值→初始化
            %
            % 语法：
            %   initial_state = init_prof.load_and_interpolate(filename, sheetname)
            %
            % 输入参数：
            %   filename  - Excel文件路径
            %   sheetname - 工作表名（默认为'sheet1'）
            %
            % 输出参数：
            %   initial_state - 初始化的初始条件结构体（不会变化）

            if nargin < 3
                sheetname = 'sheet1';
            end

            fprintf('\n═══════════════════════════════════════════════════════════════\n');
            fprintf('初始剖面加载与插值\n');
            fprintf('═══════════════════════════════════════════════════════════════\n\n');

            % ===================================================================
            % 第1步：读取Excel数据
            % ===================================================================
            fprintf('[第1步] 从Excel读取初始剖面数据\n');
            fprintf('───────────────────────────────────────────────────────────────\n');

            [init_data, metadata] = io.read_excel_data(filename, sheetname);

            fprintf('  ✓ 读取成功: %s\n', filename);
            fprintf('    原始大小: %d 行 × %d 列\n', metadata.original_size(1), metadata.original_size(2));
            fprintf('    清理后:   %d 行 × %d 列\n', metadata.final_size(1), metadata.final_size(2));

            % ===================================================================
            % 第2步：提取数据列
            % ===================================================================
            fprintf('\n[第2步] 提取并验证数据列\n');
            fprintf('───────────────────────────────────────────────────────────────\n');

            % 提取必要的列
            depths_input = init_data{:, 1};  % 深度 (cm)
            fprintf('  ✓ 深度: %d 个测量点\n', length(depths_input));
            fprintf('    范围: [%.0f, %.0f] cm\n', min(depths_input), max(depths_input));

            % 提取可选的列（如果存在）
            nh4_input = obj.extract_column(init_data, 'NH4+', 2, zeros(size(depths_input)));
            no3_input = obj.extract_column(init_data, 'NO3-', 3, zeros(size(depths_input)));
            moisture_input = obj.extract_column(init_data, 'h or θ', 4, depths_input);
            temp_input = obj.extract_column(init_data, 'T', 5, depths_input);
            pressure_input= obj.extract_column(init_data, 'Pg', 6, depths_input) * 10;  % Pa → 内部单位 (g·cm⁻¹·s⁻²)，1 Pa = 10 g·cm⁻¹·s⁻²
            bd_input = obj.extract_column(init_data, 'BD', 7, depths_input);

            % ===================================================================
            % 第3步：准备目标网格
            % ===================================================================
            fprintf('\n[第3步] 准备目标网格\n');
            fprintf('───────────────────────────────────────────────────────────────\n');

            % 根据ModelSettings构建目标深度（表层→底层）
            NL = obj.ModelSettings.Grid.NL;
            %Tot_Depth = obj.ModelSettings.Grid.Tot_Depth;

            depths_target = zeros(NL + 1, 1);
            depths_target(1) = 0;
            cumulative_depth = 0;
            for j = NL:-1:1
                cumulative_depth = cumulative_depth + obj.ModelSettings.Grid.DeltZ(j);
                depths_target(NL + 2 - j) = cumulative_depth;
            end

            fprintf('  ✓ 目标网格: %d 个节点\n', length(depths_target));
            fprintf('    范围: [%.0f, %.0f] cm\n', depths_target(1), depths_target(end));

            % ===================================================================
            % 第4步：执行各变量插值
            % ===================================================================
            fprintf('\n[第4步] 执行变量插值\n');
            fprintf('───────────────────────────────────────────────────────────────\n');

            % 4.1 温度插值
            fprintf('\n  4.1 温度插值\n');
            T_target = init.calcSoilTemp(depths_input, temp_input, depths_target);

            % 4.2 含水量插值（高级方法）
            fprintf('\n  4.2 含水量插值（基质势方法）\n');
            % Restart profiles exported from COMPASS are already on the model
            % grid and column 4 is volumetric water content. In that case,
            % use theta directly so the restart conserves the soil water
            % storage exactly. The theta->h->theta path below is retained for
            % measured profiles that require interpolation across depths.
            is_theta_profile = all(isfinite(moisture_input)) && ...
                all(moisture_input >= -1e-9) && all(moisture_input <= 1 + 1e-9);
            same_depth_grid = numel(depths_input) == numel(depths_target) && ...
                max(abs(depths_input(:) - depths_target(:))) < 1e-9;
            restart_theta_direct = is_theta_profile && same_depth_grid;
            if restart_theta_direct
                fprintf('  Restart grid detected: using input theta directly to conserve storage.\n');
                theta_target = moisture_input(:);
            else
                theta_target = init.calcMoistureMatricHead(depths_input, moisture_input, ...
                    depths_target, T_target, obj.SoilProperties, obj.ModelSettings);
            end

            % 4.3 气压初始化
            fprintf('\n  4.3 气压初始化\n');
            p_gas_target = init.interpLinearWithExtrapolation(depths_input, pressure_input, depths_target);
           
            % 4.4 氮浓度插值+单位转换
            fprintf('\n  4.4 氮浓度插值与单位转换\n');

            % 处理容重：优先使用初始条件中的容重，否则使用 SoilProperties 中的值
            if all(bd_input == 0) || all(isnan(bd_input)) || isempty(bd_input)
                % 初始条件中没有容重数据，直接使用 SoilProperties 中的层级容重
                fprintf('  ⓘ 初始条件中无容重数据，使用 SoilProperties 中的值\n');

                % SoilProperties.Bulk_density 已经是倒置版本（底层→表层）
                % 但 depths_target 是表层→底层，需要倒置 SoilProperties 的容重
                bd_soil_original = flip(obj.SoilProperties.Bulk_density);  % 倒置回表层→底层

                % 构建节点级容重 (NL+1)×1，与 depths_target 对应
                % 节点i的容重 = 相邻两层的平均值（端点节点除外）
                bd_target = zeros(NL + 1, 1);
                bd_target(1) = bd_soil_original(1);        % 表层节点，取第1层
                for i = 2:NL
                    % 节点i在第(i-1)层和第i层之间，取平均值
                    bd_target(i) = (bd_soil_original(i-1) + bd_soil_original(i)) / 2;
                end
                bd_target(NL + 1) = bd_soil_original(NL);  % 底层节点，取第NL层
                fprintf('    容重范围: [%.2f, %.2f] g/cm³\n', min(bd_target), max(bd_target));
            else
                % 初始条件中有容重数据，进行插值
                fprintf('  ✓ 使用初始条件中的容重数据进行插值\n');
                bd_target = init.interpLinearWithExtrapolation(depths_input, bd_input, depths_target);
            end

            [nh4_target, no3_target] = init.calcNitrogenConcentration(...
                depths_input, nh4_input, no3_input, depths_target, bd_target, obj.ModelSettings);

            % ===================================================================
            % 第5步：物理验证与修正
            % ===================================================================
            fprintf('\n[第5步] 物理验证与修正\n');
            fprintf('───────────────────────────────────────────────────────────────\n');

            % 计算对应的基质势（按层逐个计算，使用倒置的VG参数）
            theta_s_vec = flip(obj.SoilProperties.Theta_s);
            theta_r_vec = flip(obj.SoilProperties.Theta_r);
            alpha_vec = flip(obj.SoilProperties.Alpha);
            n_vec = flip(obj.SoilProperties.n);

            h_target = zeros(NL+1, 1);
            if restart_theta_direct
                fprintf('  Restart theta mode: recalculating h from conserved theta using VG parameters.\n');
            end
            for i = 1:NL+1
                % 节点i对应VG参数 min(max(i-1, 1), NN)
                node_idx = min(max(i - 1, 1), length(theta_s_vec));
                h_target(i) = equations.van_genuchten_theta_to_h(...
                    theta_target(i), theta_s_vec(node_idx), theta_r_vec(node_idx), alpha_vec(node_idx), n_vec(node_idx));
            end

            % 执行验证
            [theta_target, h_target, corrections] = init.validatePhysicalRange(...
                theta_target, h_target, obj.SoilProperties, obj.ModelSettings);

            % ===================================================================
            % 第6步：坐标转换（表层→底层 变为 底层→表层）
            % ===================================================================
            fprintf('\n[第6步] 坐标转换\n');
            fprintf('───────────────────────────────────────────────────────────────\n');

            fprintf('  🔄 执行数据倒置 (表层→底层 变为 底层→表层)...\n');

            T_reversed = flip(T_target);
            theta_reversed = flip(theta_target);
            h_reversed = flip(h_target);
            p_gas_reversed = flip(p_gas_target);
            nh4_reversed = flip(nh4_target);
            no3_reversed = flip(no3_target);
            bd_reversed = flip(bd_target);

            fprintf('  ✓ 倒置完成\n');
            fprintf('    现在索引0代表底部（深度最大）\n');

            % ===================================================================
            % 第7步：构建状态结构体
            % ===================================================================
            fprintf('\n[第7步] 初始化状态结构体\n');
            fprintf('───────────────────────────────────────────────────────────────\n');

            % 转换含水量为 NL×2 格式（每层的两个相邻节点值）
            fprintf('\n[第7步] 含水量维度转换 (NL+1→NL×2)...\n');
            Theta_L_init = zeros(NL, 2);
            for i = 1:NL
                % 第i层对应节点i和i+1
                Theta_L_init(i, 1) = theta_reversed(i);      % 第i个节点
                Theta_L_init(i, 2) = theta_reversed(i + 1);  % 第i+1个节点
            end
            fprintf('  ✓ 含水量维度转换完成 (NL+1=%d → NL×2=%d×2)\n', NL+1, NL);

            % 创建输出初始条件结构体
            initial_state = struct();

            % 水热状态
            initial_state.theta = theta_reversed;           % 节点值 (NL+1)×1
            initial_state.Theta_L = Theta_L_init;           % 层值 NL×2
            initial_state.h = h_reversed;
            initial_state.T = T_reversed;

            % 气体状态
            initial_state.p_gas = p_gas_reversed;

            % 土壤结构参数
            initial_state.bd = bd_reversed;                 % 容重 (g/cm³)，用于后续重新计算热学参数

            % 氮状态
            initial_state.nh4 = nh4_reversed;
            initial_state.no3 = no3_reversed;

            % 元数据
            initial_state.metadata = struct();
            initial_state.metadata.depths_original = depths_input;
            initial_state.metadata.depths_target = depths_target;
            initial_state.metadata.corrections = corrections;
            initial_state.metadata.excel_file = filename;
            initial_state.metadata.excel_sheet = sheetname;

            fprintf('  ✓ 初始条件结构体初始化完成\n');
            fprintf('    节点变量: theta (NL+1)×1, h (NL+1)×1, T (NL+1)×1, p_gas (NL+1)×1\n');
            fprintf('    层变量: Theta_L (NL×2), NH4+ (NL+1)×1, NO3- (NL+1)×1\n');
            fprintf('    土壤参数: bd (NL+1)×1 - 容重用于重新计算热学参数\n');

            fprintf('\n═══════════════════════════════════════════════════════════════\n');
            fprintf('初始剖面加载完成！\n');
            fprintf('═══════════════════════════════════════════════════════════════\n\n');

        end

        function val = extract_column(obj, table_data, col_name, col_idx, default_val)
            %EXTRACT_COLUMN 从表格中安全提取列

            try
                if col_idx <= width(table_data)
                    val = table_data{:, col_idx};
                    if all(isnan(val)) || length(val) == 0
                        fprintf('  ⚠️  %s 列为空或全NaN，使用默认值\n', col_name);
                        val = default_val;
                    else
                        fprintf('  ✓ %s: %d 个值\n', col_name, length(val));
                    end
                else
                    fprintf('  ⚠️  %s 列不存在，使用默认值\n', col_name);
                    val = default_val;
                end
            catch
                fprintf('  ⚠️  %s 列提取失败，使用默认值\n', col_name);
                val = default_val;
            end
        end

    end
end
