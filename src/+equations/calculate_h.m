function Inith = calculate_h(Theta_s, Theta_r, InitX0, m, n, Alpha)
    % 初始化输出向量
    Inith = zeros(size(InitX0));
    
    % 计算有效饱和度相关项，添加 max 处理防止 InitX0 < Theta_r
    Se_inv = (Theta_s - Theta_r) ./ (max(InitX0, Theta_r + eps) - Theta_r + eps);
    
    % 计算 term3：((theta_s-theta_r)/(theta-theta_r))^(1/m) - 1
    term3 = Se_inv.^(1./m) - 1;
    
    % 关键：防止由于浮点数误差产生的微小负值导致复数结果
    term3 = max(term3, 0); 
    
    % 计算水头
    Inith = -(term3.^(1./n)) ./ Alpha;
    
    % 处理边界情况：饱和状态
    Inith(InitX0 >= Theta_s) = -1e-6; 
    
    % 处理边界情况：极度干旱状态（低于残余含水量）
    Inith(InitX0 <= Theta_r) = -1e7;
end