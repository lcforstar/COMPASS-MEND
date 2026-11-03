function [Rl, Ztot] = initialRootBiomass(RTB, DeltZ_R, rroot, ML, landcoverClass)
% INITIALROOTBIOMASS  根系生物量初始化，计算各层根长密度
%
%   与 vegetation.hydraulics.calculateRootAreaIndex1 同出一源（Jackson 1996），
%   但计算的物理量和用途不同：
%
%     本函数             → 根长密度 Rl [cm/cm³]，供 ebal 根系吸水使用
%     calculateRootAreaIndex1 → 根面积指数 RAI [m²/m²]，供 Plantrootconductivity 土-根导水使用
%
%   同为 Jackson 1996 β 累积分布，但 β 查表方式不同：
%     本函数：按地表覆盖类型字符串（landcoverClass）查 β
%     calculateRootAreaIndex1：按 biome 编号整数（root_type 1-11）查 β
%
%   参考文献：
%     Jackson et al. (1996) A global analysis of root distributions
%       for terrestrial biomes. Oecologia 108: 389-411.
%
%   层索引约定：层1 = 地表，层ML = 最深层（top-down，与 DeltZ_R 相同方向）
%   输出 Rl 和 Ztot 均遵循此约定，Rl(1) 为表层最大值。
%
%   输入参数：
%     RTB              - 初始根系总生物量 [gDM/m²]
%     DeltZ_R          - 根系层厚度数组，从地表到底层 [cm]（top-down）
%     rroot            - 根系半径 [cm]
%     ML               - 根系层数
%     landcoverClass   - 土地覆盖类型字符串
%
%   输出参数：
%     Rl               - 各层根长密度，层1=地表 [cm/cm³]
%     Ztot             - 各层底部累积深度，层1=地表 [cm]

    % ========== 输入验证 ==========
    assert(RTB > 0, 'RTB must be positive [gDM/m²], got: %.4f', RTB);
    assert(rroot > 0, 'rroot must be positive [cm], got: %.6f', rroot);
    assert(ML >= 1, 'ML must be >= 1, got: %d', ML);
    
    DeltZ_R = reshape(DeltZ_R(1:ML), [], 1);
    assert(all(DeltZ_R > 0), 'All layer thicknesses must be positive');

    % ========== 物理常数 ==========
    root_den = 250e3;       % 根系密度 [gDM/m³]

    % ========== 单位转换 ==========
    rroot_m  = rroot / 100;         % cm → m

    % ========== Beta参数（Jackson 1996） ==========
    if startsWith(landcoverClass, 'Permanent Wetlands')
        beta = 0.993;
    elseif startsWith(landcoverClass, 'Evergreen Broadleaf')
        beta = 0.993;
    elseif startsWith(landcoverClass, 'Deciduous Broadleaf')
        beta = 0.993;
    elseif startsWith(landcoverClass, 'Mixed Forests')
        beta = 0.993;
    elseif startsWith(landcoverClass, 'Evergreen Needleleaf')
        beta = 0.993;
    elseif startsWith(landcoverClass, 'Croplands')
        beta = 0.943;
    elseif startsWith(landcoverClass, 'Open Shrublands')
        beta = 0.966;
    elseif startsWith(landcoverClass, 'Closed Shrublands')
        beta = 0.966;
    elseif startsWith(landcoverClass, 'Savannas')
        beta = 0.943;
    elseif startsWith(landcoverClass, 'Woody Savannas')
        beta = 0.943;
    elseif startsWith(landcoverClass, 'Grassland')
        beta = 0.943;
    else
        beta = 0.943;
        warning('Landcover class "%s" not recognized. Using default beta=0.943', landcoverClass);
    end

    % ========== 总根长密度（基于RTB） ==========
    % RTB [gDM/m²] / root_den [gDM/m³] = 总根体积 [m³/m²] = [m]
    % 除以单根截面积 [m²] = 总根长 [m/m²]
    Rltot = RTB / root_den / (pi * rroot_m^2);  % [m/m²]

    % ========== 累积深度（top-down，与层索引一致） ==========
    Ztot = cumsum(DeltZ_R);  % [cm]，Ztot(i) = 第i层底部深度

    % ========== Jackson 模型：各层根系比例 ==========
    % Y(d) = 1 - beta^d（d单位：cm）
    % 第i层比例 ri = Y(Ztot(i)) - Y(Ztot(i-1))
    ri = zeros(ML, 1);
    for i = 1:ML
        if i == 1
            ri(i) = 1 - beta^Ztot(i);
        else
            ri(i) = beta^Ztot(i-1) - beta^Ztot(i);
        end
    end

    % ========== 质量守恒：归一化处理 ==========
    total_fraction = sum(ri);
    truncation_pct = (1 - total_fraction) * 100;
    
    if truncation_pct > 1.0
        warning(['Root profile truncated: %.1f%% of RTB not assigned to any layer.\n' ...
                 '         Renormalizing ri to conserve mass within root domain.\n' ...
                 '         Consider increasing root_max_depth for %s (beta=%.3f).'], ...
                 truncation_pct, landcoverClass, beta);
    end
    % 归一化：确保研究域内质量守恒
    ri = ri / total_fraction;

    % ========== 各层根长密度（top-down，无需flipud） ==========
    DeltZ_R_m = DeltZ_R / 100;                      % cm → m
    Rl_m_per_m3 = ri .* Rltot ./ DeltZ_R_m;         % [m/m³]

    % 单位转换：m/m³ → cm/cm³
    % 1 m/m³ = 1e-4 cm/cm³
    Rl_top = Rl_m_per_m3 * 1e-4;                        % [cm/cm³]
    Rl=flipud(Rl_top);
  
end