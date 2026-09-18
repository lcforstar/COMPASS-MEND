function L = Latent(SoilVariables, cfg)
% Latent - 液态水气化潜热 / Latent heat of vaporization of water (kJ/g)
%
% 输入 / Inputs:
%   SoilVariables - 包含 TT 的结构体
%   cfg           - ModelSettings，含 Grid.NN

NN = cfg.Grid.NN;
TT = SoilVariables.TT;

L = zeros(1, NN);
for MN = 1:NN
    L(MN) = (2.501e6 - 2369.2 * TT(MN)) / 1000;
end
end
