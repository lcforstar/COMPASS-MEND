function [Gamma_h, Gamma_hh] = updateGammaH(i, SoilVariables)

% 获取土壤滞后模型常数
% io.getSoilConstants() 返回原 STEMMUS-SCOPE io 包中的土壤物理常数结构体，
% 若迁移至本项目请替换为直接赋值或从 SoilProperties 中提取：
%   hd: 极干（气干）阈值水头 (cm)，|hh| >= |hd| 时滞后系数 Gamma = 0
%       典型值约为 -1e7 cm（接近气干状态）
%   hm: 中间阈值水头 (cm)，|hm| < |hh| < |hd| 时 Gamma 线性插值 (0~1)
%       典型值约为 -9899 cm（接近凋萎点）
SoilConstants = io.getSoilConstants();

hd = SoilConstants.hd;
hm = SoilConstants.hm;

Gamma_h  = SoilVariables.Gamma_h;
Gamma_hh = SoilVariables.Gamma_hh;

hh = SoilVariables.hh;

if abs(hh(i)) >= abs(hd)
    Gamma_h(i)  = 0;
    Gamma_hh(i) = Gamma_h(i);
elseif abs(hh(i)) >= abs(hm)
    % Gamma_h(i)=1-log(abs(hh(i)))/log(abs(hm));
    Gamma_h(i)  = log(abs(hd) / abs(hh(i))) / log(abs(hd) / abs(hm));
    Gamma_hh(i) = Gamma_h(i);
else
    Gamma_h(i)  = 1;
    Gamma_hh(i) = Gamma_h(i);
end

end
