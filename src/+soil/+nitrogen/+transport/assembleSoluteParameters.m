function [R1, R1O, E1, B1, F1, G1, R2, R2O, E2, B2, F2, G2, cvSnit, cvSmin, cvSden, cvSinkS, cvSinkS1, cvSvol, ...
          QNH4, QNO3, DNH4DZ, DNO3DZ, QNH4_h, QNH4_T, QNH4_A, QNH4_convection, ...
          QNH4_convection_h, QNH4_convection_T, QNH4_convection_A, QNH4_dispersion, ...
          QNH4_dispersion_h, QNH4_dispersion_T, QNH4_dispersion_A, QNO3_h, ...
          QNO3_T, QNO3_A, QNO3_convection, QNO3_convection_h, QNO3_convection_T, ...
          QNO3_convection_A, QNO3_dispersion, QNO3_dispersion_h, QNO3_dispersion_T, ...
          QNO3_dispersion_A] = assembleSoluteParameters(NL, RHO_b_Kd, Dsh, Dsh_h, Dsh_T, Dsh_A, QL_nodes, QL_nodes_h, QL_nodes_T, QL_nodes_A, ...
          Theta_Solute_L, Theta_Solute_LL, Snit, SVnit, SVden, Smin, Svol, SinkS, Sden, ...
          cvSnit, cvSvol, cvSmin, cvSden, cvSinkS, cvSinkS1, DeltZ, R1, R1O, E1, B1, F1, G1, ...
          R2, R2O, E2, B2, F2, G2, Shys, Nn, Nno3, KT, Sdis, SorgNO3)

    % =====================================================================
    % 0. 数据预处理：统一转换为列向量 (针对图片1中的行向量报错进行修复)
    % =====================================================================
    RHO_b_Kd = RHO_b_Kd(:);     % 强制转为 [NL x 1]
    % [2026-05-17 dissolution pool] backward-compat: Sdis 可选 (43rd param), 默认 0
    if nargin < 43 || isempty(Sdis)
        Sdis = zeros(NL, 2);
    end
    if nargin < 44 || isempty(SorgNO3)
        SorgNO3 = zeros(NL, 2);
    end
    DeltZ  = DeltZ(:);      % 强制转为 [NL x 1]
    Nn     = Nn(:);         % 强制转为 [NN x 1]
    Nno3   = Nno3(:);       % 强制转为 [NN x 1]
    QL_nodes   = QL_nodes(:);   % 确保 [NN x 1]
    QL_nodes_h = QL_nodes_h(:);
    QL_nodes_T = QL_nodes_T(:);
    QL_nodes_A = QL_nodes_A(:);

    % =====================================================================
    % 1. 向量化计算系数 (系数矩阵 [NL x 2])
    % =====================================================================
    % 利用广播机制：[NL x 1] ./ [NL x 2] -> 自动按行匹配
    R1O = 1 + RHO_b_Kd ./ Theta_Solute_L;
    R1  = 1 + RHO_b_Kd ./ Theta_Solute_LL;
    E1  = Dsh(:,:,1);
    % 矩阵化映射 QL_nodes：第一列取 1:NL, 第二列取 2:NL+1
    B1  = [QL_nodes(1:NL), QL_nodes(2:NL+1)]; 
    F1  = -SVnit;                              % [2026-05-16] SVnit=0 (MM 改显式), F1 闲置. 保留接口
    % [2026-05-17 dissolution pool] +Sdis: NH4 颗粒肥溶解源 (granule → dissolved)
    %   Shys: urea→NH4, Smin: SOM→NH4, +Sdis: granule→NH4, -Svol/Snit/SinkS: NH4 汇
    G1  = Shys + Smin + Sdis - Svol - Snit - SinkS(:,:,1);

    % NO3 项 (R2 = 1)
    R2  = ones(NL, 2);
    R2O = ones(NL, 2);
    E2  = Dsh(:,:,2);
    B2  = B1;
    F2  = -SVden;                  % [2026-05-16] SVden=0 (MM 改显式), F2 闲置. 保留接口
    G2  = 0.98 * Snit - Sden - SinkS(:,:,2) + SorgNO3;

    % =====================================================================
    % 2. 向量化梯度计算 (干掉循环)
    % =====================================================================
    DNH4DZ = (Nn(2:end) - Nn(1:end-1)) ./ DeltZ;
    DNO3DZ = (Nno3(2:end) - Nno3(1:end-1)) ./ DeltZ;

    % =====================================================================
    % 3. 彻底干掉 mean 函数 (消除 Profiler 中第1、3慢瓶颈)
    % =====================================================================
    % 算术平均 [NL x 1]
    E1BAR = (E1(:,1) + E1(:,2)) * 0.5;
    E2BAR = (E2(:,1) + E2(:,2)) * 0.5;
    
    % 三维张量导数项平均
    Dsh_h1BAR = (Dsh_h(:,1,1) + Dsh_h(:,2,1)) * 0.5;
    Dsh_T1BAR = (Dsh_T(:,1,1) + Dsh_T(:,2,1)) * 0.5;
    Dsh_A1BAR = (Dsh_A(:,1,1) + Dsh_A(:,2,1)) * 0.5;
    
    Dsh_h2BAR = (Dsh_h(:,1,2) + Dsh_h(:,2,2)) * 0.5;
    Dsh_T2BAR = (Dsh_T(:,1,2) + Dsh_T(:,2,2)) * 0.5;
    Dsh_A2BAR = (Dsh_A(:,1,2) + Dsh_A(:,2,2)) * 0.5;

    % =====================================================================
    % 4. 完全向量化通量计算 (核心加速点：移除原有的 140 行循环计算)
    % =====================================================================
    % 预提取左右侧节点数据
    Nn_L = Nn(1:NL);     Nn_R = Nn(2:NL+1);
    Nno3_L = Nno3(1:NL); Nno3_R = Nno3(2:NL+1);
    
    QL_L = QL_nodes(1:NL); QL_R = QL_nodes(2:NL+1);
    QLh_L = QL_nodes_h(1:NL); QLh_R = QL_nodes_h(2:NL+1);
    QLT_L = QL_nodes_T(1:NL); QLT_R = QL_nodes_T(2:NL+1);
    QLA_L = QL_nodes_A(1:NL); QLA_R = QL_nodes_A(2:NL+1);

    % --- NH4+ 通量运算 ---
    % 弥散项预计算 (左右相同，均值即本身)
    disp1 = -E1BAR .* DNH4DZ; 
    conv1_L = QL_L .* Nn_L; conv1_R = QL_R .* Nn_R;
    
    QNH4 = ( (conv1_L + disp1) + (conv1_R + disp1) ) * 0.5;
    QNH4_convection = (conv1_L + conv1_R) * 0.5;
    QNH4_dispersion = disp1;
    
    % NH4+ 导数项
    disp1_h = -Dsh_h1BAR .* DNH4DZ;
    QNH4_h = ( (QLh_L .* Nn_L + disp1_h) + (QLh_R .* Nn_R + disp1_h) ) * 0.5;
    QNH4_convection_h = (QLh_L .* Nn_L + QLh_R .* Nn_R) * 0.5;
    QNH4_dispersion_h = disp1_h;

    disp1_T = -Dsh_T1BAR .* DNH4DZ;
    QNH4_T = ( (QLT_L .* Nn_L + disp1_T) + (QLT_R .* Nn_R + disp1_T) ) * 0.5;
    QNH4_convection_T = (QLT_L .* Nn_L + QLT_R .* Nn_R) * 0.5;
    QNH4_dispersion_T = disp1_T;

    disp1_A = -Dsh_A1BAR .* DNH4DZ;
    QNH4_A = ( (QLA_L .* Nn_L + disp1_A) + (QLA_R .* Nn_R + disp1_A) ) * 0.5;
    QNH4_convection_A = (QLA_L .* Nn_L + QLA_R .* Nn_R) * 0.5;
    QNH4_dispersion_A = disp1_A;

    % --- NO3- 通量运算 ---
    disp2 = -E2BAR .* DNO3DZ;
    conv2_L = QL_L .* Nno3_L; conv2_R = QL_R .* Nno3_R;

    QNO3 = ( (conv2_L + disp2) + (conv2_R + disp2) ) * 0.5;
    QNO3_convection = (conv2_L + conv2_R) * 0.5;
    QNO3_dispersion = disp2;

    % NO3- 导数项
    disp2_h = -Dsh_h2BAR .* DNO3DZ;
    QNO3_h = ( (QLh_L .* Nno3_L + disp2_h) + (QLh_R .* Nno3_R + disp2_h) ) * 0.5;
    QNO3_convection_h = (QLh_L .* Nno3_L + QLh_R .* Nno3_R) * 0.5;
    QNO3_dispersion_h = disp2_h;

    disp2_T = -Dsh_T2BAR .* DNO3DZ;
    QNO3_T = ( (QLT_L .* Nno3_L + disp2_T) + (QLT_R .* Nno3_R + disp2_T) ) * 0.5;
    QNO3_convection_T = (QLT_L .* Nno3_L + QLT_R .* Nno3_R) * 0.5;
    QNO3_dispersion_T = disp2_T;

    disp2_A = -Dsh_A2BAR .* DNO3DZ;
    QNO3_A = ( (QLA_L .* Nno3_L + disp2_A) + (QLA_R .* Nno3_R + disp2_A) ) * 0.5;
    QNO3_convection_A = (QLA_L .* Nno3_L + QLA_R .* Nno3_R) * 0.5;
    QNO3_dispersion_A = disp2_A;

    % =====================================================================
    % 5. 向量化积分计算 (干掉最后的累加循环)
    % =====================================================================
    % 预计算每个单元层面的平均源项
    Svol_avg   = (Svol(:,1) + Svol(:,2)) * 0.5;
    Snit_avg   = (Snit(:,1) + Snit(:,2)) * 0.5;
    Smin_avg   = (Smin(:,1) + Smin(:,2)) * 0.5;
    Sden_avg   = (Sden(:,1) + Sden(:,2)) * 0.5;
    SinkS1_avg = (SinkS(:,1,1) + SinkS(:,2,1)) * 0.5;
    SinkS2_avg = (SinkS(:,1,2) + SinkS(:,2,2)) * 0.5;

    % 利用矩阵转置乘法 DeltZ' * S 直接求积分，结果存入当前 KT 时刻
    cvSvol(KT)   = DeltZ' * Svol_avg;
    cvSnit(KT)   = DeltZ' * Snit_avg;
    cvSmin(KT)   = DeltZ' * Smin_avg;
    cvSden(KT)   = DeltZ' * Sden_avg;
    cvSinkS(KT)  = DeltZ' * SinkS1_avg;
    cvSinkS1(KT) = DeltZ' * SinkS2_avg;

end
