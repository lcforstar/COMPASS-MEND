function [theta_valid, h_valid, corrections] = validatePhysicalRange(theta_input, h_input, ...
                                                                    SoilProperties, ModelSettings)
    %VALIDATEPHYSICALRANGE 物理范围验证与自动修正
    %
    % 语法：
    %   [theta_valid, h_valid, corrections] = validatePhysicalRange(theta_input, h_input, ...
    %       SoilProperties, ModelSettings)
    %
    % 输入参数：
    %   theta_input   - 输入含水量 (体积分数)
    %   h_input       - 对应的基质势 (cm)
    %   SoilProperties- 土壤性质对象（含分布式VG参数）
    %   ModelSettings - 模型设置
    %
    % 输出参数：
    %   theta_valid   - 验证后的含水量
    %   h_valid       - 验证后的基质势
    %   corrections   - 修正信息结构体
    %
    % 功能说明：
    %   1. 检查含水量是否在 [θ_r, θ_s] 范围内
    %   2. 超过饱和含水量 θ_s 时，修正为田间含水量 θ_f
    %   3. 低于残余含水量 θ_r 时，修正为 θ_r + ε
    %   4. 重新计算对应的基质势
    %
    % 修正策略：
    %   θ > θ_s    → θ = θ_f (田间含水量)
    %   θ < θ_r    → θ = θ_r + 1e-6
    %   θ_r ≤ θ ≤ θ_s → 保持不变

    theta_valid = theta_input;
    h_valid = h_input;

    % 获取分布式van Genuchten参数（按节点分布）
    theta_s_vec = flip(SoilProperties.Theta_s);
    theta_r_vec = flip(SoilProperties.Theta_r);
    alpha_vec = flip(SoilProperties.Alpha);
    n_vec = flip(SoilProperties.n);

    fprintf('\n🔍 执行物理范围验证...\n');

    % 初始化修正记录
    corrections = struct();
    corrections.n_exceed_saturation = 0;
    corrections.n_below_residual = 0;
    corrections.exceed_sat_indices = [];
    corrections.below_res_indices = [];

    % 获取田间持水量（按层，NL维）
    if  ~isempty(SoilProperties.Field_mc)
        Field_mc_vec = flip(SoilProperties.Field_mc);
    else
        Field_mc_vec = [];  % 如果没有，则在循环中计算
    end

    % 对每个点进行验证
    for i = 1:length(theta_valid)
        theta = theta_valid(i);
        exceed_sat = false;
        below_res = false;

        % 获取节点i对应的VG参数
        node_idx = min(max(i - 1, 1), length(theta_s_vec));
        theta_s = theta_s_vec(node_idx);
        theta_r = theta_r_vec(node_idx);
        alpha = alpha_vec(node_idx);
        n = n_vec(node_idx);

        % 田间含水量：优先使用SoilProperties中的值，否则计算
        if ~isempty(Field_mc_vec) && node_idx <= length(Field_mc_vec)
            theta_f = Field_mc_vec(node_idx);
        else
            theta_f = 0.75 * theta_s;  % 备用值
        end

        % 检查1：超过饱和含水量
        if theta > theta_s
            theta_valid(i) = theta_f;  % 修正为田间含水量
            corrections.n_exceed_saturation = corrections.n_exceed_saturation + 1;
            corrections.exceed_sat_indices = [corrections.exceed_sat_indices, i];
            exceed_sat = true;
        end

        % 检查2：低于残余含水量
        if theta < theta_r
            theta_valid(i) = theta_r + 1e-6;
            corrections.n_below_residual = corrections.n_below_residual + 1;
            corrections.below_res_indices = [corrections.below_res_indices, i];
            below_res = true;
        end

        % 如果有修正，重新计算基质势
        if exceed_sat || below_res
            h_valid(i) = equations.van_genuchten_theta_to_h(theta_valid(i), theta_s, theta_r, alpha, n);
        end
    end

    % 打印修正摘要
    fprintf('  ✓ 验证完成\n');
    fprintf('    含水量范围: [%.4f, %.4f] cm³/cm³\n', min(theta_valid), max(theta_valid));
    fprintf('    允许范围: [%.4f, %.4f] (θ_r ~ θ_s)\n', theta_r, theta_s);

    if corrections.n_exceed_saturation > 0
        fprintf('    ⚠️  超过饱和: %d 个点，已修正为田间含水量 (%.4f)\n', ...
            corrections.n_exceed_saturation, theta_f);
    end

    if corrections.n_below_residual > 0
        fprintf('    ⚠️  低于残余: %d 个点，已修正为残余含水量 (%.4f)\n', ...
            corrections.n_below_residual, theta_r);
    end

    corrections.total_corrections = corrections.n_exceed_saturation + corrections.n_below_residual;

    if corrections.total_corrections == 0
        fprintf('    ✓ 全部数据都在物理范围内，无需修正\n');
    end

end
