function p_gas = calcGasPressure(depths, method, varargin)
    %CALCGASPRESSURE 土壤气体（空气）压力初始化
    %
    % 语法：
    %   p_gas = calcGasPressure(depths, method, ...)
    %
    % 输入参数：
    %   depths   - 目标深度 (cm)
    %   method   - 初始化方法：
    %              'constant'   - 常数值（所有深度相同）
    %              'linear'     - 线性变化
    %              'measured'   - 从测量数据插值
    %   varargin - 附加参数（取决于method）
    %
    % 不同方法的附加参数：
    %   method='constant': p_gas = calcGasPressure(depths, 'constant', p_value)
    %     p_value     - 气压常数值 (Pa)，默认为 101325 Pa (标准大气压)
    %
    %   method='linear': p_gas = calcGasPressure(depths, 'linear', p_top, p_bottom)
    %     p_top       - 表层气压 (Pa)
    %     p_bottom    - 底层气压 (Pa)
    %
    %   method='measured': p_gas = calcGasPressure(depths, 'measured', depths_input, p_input)
    %     depths_input - 测量深度 (cm)
    %     p_input      - 测量气压 (Pa)
    %
    % 输出参数：
    %   p_gas    - 各深度的气压 (Pa)
    %
    % 功能说明：
    %   在大多数土壤学应用中，气压被假设为恒定的大气压力。
    %   本函数提供多种初始化方式以满足不同的模型需求。
    %
    % 标准值：
    %   大气压 = 101325 Pa (标准条件下)
    %   气压随高度变化 ≈ -12 Pa/m (对流层)

    depths = depths(:);

    switch lower(method)
        case 'constant'
            % ===== 方法1：恒定气压 =====
            if nargin < 3 || isempty(varargin{1})
                p_value = 101325;  % 默认：标准大气压 (Pa)
            else
                p_value = varargin{1};
            end

            p_gas = p_value * ones(size(depths));

            fprintf('  ℹ️  气压初始化：恒定值 = %.0f Pa\n', p_value);

        case 'linear'
            % ===== 方法2：线性梯度 =====
            if nargin < 4
                error('线性方法需要2个参数：p_top, p_bottom');
            end

            p_top = varargin{1};
            p_bottom = varargin{2};

            % 线性插值
            depth_min = min(depths);
            depth_max = max(depths);
            delta_depth = depth_max - depth_min;

            if delta_depth < 1e-10
                p_gas = p_top * ones(size(depths));
            else
                p_gas = p_top + (p_bottom - p_top) * (depths - depth_min) / delta_depth;
            end

            fprintf('  ℹ️  气压初始化：线性梯度\n');
            fprintf('      表层(%.0f cm): %.0f Pa\n', depth_min, p_top);
            fprintf('      底层(%.0f cm): %.0f Pa\n', depth_max, p_bottom);

        case 'measured'
            % ===== 方法3：从测量数据插值 =====
            if nargin < 4
                error('measured方法需要2个参数：depths_input, p_input');
            end

            depths_input = varargin{1}(:);
            p_input = varargin{2}(:);

            % 线性插值+外推
            p_gas = init.interpLinearWithExtrapolation(depths_input, p_input, depths);

            fprintf('  ℹ️  气压初始化：从 %d 个测量点插值\n', length(depths_input));

        otherwise
            error('未知的方法: %s。应为 ''constant'', ''linear'' 或 ''measured''', method);
    end

    % 验证结果
    if any(isnan(p_gas)) || any(isinf(p_gas))
        warning('calcGasPressure:InvalidValue', '气压计算包含 NaN 或 Inf 值');
    end

end
