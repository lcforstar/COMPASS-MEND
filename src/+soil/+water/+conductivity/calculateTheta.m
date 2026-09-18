function theta = calculateTheta(subRoutine, theta_m, heat_term, gamma_hh, theta_s, theta_r, lamda, phi_s, alpha, n, m)
% calculateTheta - 按子路径计算 VG / Clapp-Hornberger 含水量
%                  Compute water content via VG or Clapp-Hornberger subroutine
%
% 函数定义 / Function Definition:
%   theta = calculateTheta(subRoutine, theta_m, heat_term, gamma_hh, theta_s, theta_r, lamda, phi_s, alpha, n, m)
%
% 输入 / Inputs:
%   subRoutine - 计算路径索引（0=lFreezing+解析, 1=Clapp, 2=VG无lFreezing）
%   theta_m    - 毛管含水量上限 / Capillary water content upper bound
%   heat_term  - 有效压力水头（可为 hh 或 hh+hh_frez）
%   gamma_hh   - 冻结分区因子 / Freezing partition factor
%   theta_s/r  - 饱和/残余含水量 / Saturated / residual water content
%   lamda      - Clapp-Hornberger 指数 / Clapp-Hornberger exponent
%   phi_s      - 进气压头 / Air-entry pressure head
%   alpha, n, m - van Genuchten 参数 / VG parameters
%
% 输出 / Outputs:
%   theta      - 体积含水量 / Volumetric water content
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-10

    switch subRoutine
        case 0
            % lFreezing=1 + 解析差分路径 / lFreezing=1 analytical path
            % heat_term = hh or hh + hh_frez
            SoilConstants = io.getSoilConstants();
            hd = SoilConstants.hd;
            theta = gamma_hh * theta_r + (theta_m - gamma_hh * theta_r) / (1 + abs(alpha * heat_term)^n)^m;
            if heat_term <= hd
                theta = 0;
            end
        case 1
            % Clapp-Hornberger 路径 / Clapp-Hornberger path
            % heat_term = hh or hh + hh_frez
            theta = theta_s * (heat_term / phi_s)^(-1 * lamda);
            if heat_term <= -1e7
                theta = theta_r;
            elseif heat_term >= phi_s
                theta = theta_s;
            end
        case 2
            % VG 无 lFreezing 路径 / VG without lFreezing path
            % heat_term = hh
            theta = theta_r + (theta_s - theta_r) / (1 + abs(alpha * heat_term)^n)^m;
    end
end
