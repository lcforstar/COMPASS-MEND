function theta_interp = calcMoistureMatricHead(depths_input, theta_input, depths_target, ...
                                                 T_target, SoilProperties, ModelSettings)
    %CALCMOISTUREMATRICHEAD 含水量插值（基质势对数插值法）
    %
    % 语法：
    %   theta_interp = calcMoistureMatricHead(depths_input, theta_input, ...
    %       depths_target, T_target, SoilProperties, ModelSettings)
    %
    % 输入参数：
    %   depths_input   - 测量深度 (cm)
    %   theta_input    - 测量含水量 (体积分数)
    %   depths_target  - 目标深度 (cm)
    %   T_target       - 目标深度的温度 (°C)
    %   SoilProperties - 土壤性质对象（包含van Genuchten参数）
    %   ModelSettings  - 模型设置
    %
    % 输出参数：
    %   theta_interp   - 插值得到的含水量
    %
    % 功能说明：
    %   1. 测量点：θ → h（使用van Genuchten方程）
    %   2. 对数空间插值 h
    %   3. 模型节点：h → θ（考虑温度修正）
    %
    % 方法理由：
    %   含水量与基质势为强非线性关系。直接线性插值会产生物理不合理的中间值。
    %   通过基质势作为中介变量，在对数空间插值，确保数值稳定性。
    %
    % 参考：
    %   van Genuchten, M.T. (1980). A Closed-form Equation for Predicting the
    %   Hydraulic Conductivity of Unsaturated Soils. Soil Science Society of
    %   America Journal, 44: 892-898.

    % 确保输入为列向量
    depths_input = depths_input(:);
    theta_input = theta_input(:);
    depths_target = depths_target(:);
    T_target = T_target(:);

    % 获取分布式van Genuchten参数（节点一一对应depths_target）
    theta_s_vec = flip(SoilProperties.Theta_s);
    theta_r_vec = flip(SoilProperties.Theta_r);
    alpha_vec = flip(SoilProperties.Alpha);
    n_vec = flip(SoilProperties.n);

    % =========================================================================
    % 第1步：测量点 θ → h（按深度对应节点选择VG参数）
    % =========================================================================
    fprintf('  🔄 第1步：计算测量点的基质势 (θ → h)...\n');

    h_measured = zeros(size(theta_input));
    for i = 1:length(depths_input)
        % 根据深度在depths_target中找到对应的节点索引
        node_idx = find_nearest_node_index(depths_input(i), depths_target);

        % 使用对应节点的VG参数
        h_measured(i) = equations.van_genuchten_theta_to_h(...
            theta_input(i), theta_s_vec(node_idx), theta_r_vec(node_idx), ...
            alpha_vec(node_idx), n_vec(node_idx));
    end

    % =========================================================================
    % 第2步：对数空间插值 h
    % =========================================================================
    fprintf('  🔄 第2步：在对数空间插值基质势...\n');

    % 转换为对数空间（取绝对值避免负数）
    log_h = log10(abs(h_measured) + 1e-10);

    % 线性插值对数基质势
    log_h_interp = init.interpLinearWithExtrapolation(depths_input, log_h, depths_target);

    % 转换回物理空间
    h_interp = -10.^log_h_interp;

    % =========================================================================
    % 第3步：模型节点 h → θ（按深度对应节点选择VG参数）
    % =========================================================================
    fprintf('  🔄 第3步：计算目标深度的含水量 (h → θ)...\n');

    enable_thermal = ModelSettings.lThmrlefc;

    theta_interp = zeros(size(h_interp));
    for i = 1:length(depths_target)
        % depths_target(1)=0对应第1个VG参数；depths_target(k)对应第k-1个VG参数
        node_idx = find_nearest_node_index(depths_target(i), depths_target);

        theta_interp(i) = equations.van_genuchten_h_to_theta(...
            h_interp(i), theta_s_vec(node_idx), theta_r_vec(node_idx), ...
            alpha_vec(node_idx), n_vec(node_idx), T_target(i), enable_thermal);
    end

    fprintf('  ✓ 含水量插值完成\n');
end

function node_idx = find_nearest_node_index(depth, depths_target)
    % 根据深度找到depths_target中的节点索引
    % 优先查找精确匹配，否则找最近的相邻节点

    tolerance = 1e-10;

    % 检查是否存在精确匹配的深度值
    exact_match = abs(depths_target - depth) < tolerance;
    if any(exact_match)
        node_idx = find(exact_match, 1);
    else
        % 否则找最近的相邻节点
        [~, node_idx] = min(abs(depths_target - depth));
    end

    % 对应VG参数的索引（VG参数比depths_target少一个）
    node_idx = min(max(node_idx - 1, 1), length(depths_target) - 2);
end
