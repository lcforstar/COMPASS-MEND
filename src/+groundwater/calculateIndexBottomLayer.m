function [indxBotmLayer, indxBotmLayer_R] = calculateIndexBottomLayer(soilThick, gw_Dep, GridSettings)
% calculateIndexBottomLayer - Find soil layer index containing the water table
%
%   [indxBotmLayer, indxBotmLayer_R] = calculateIndexBottomLayer(soilThick, gw_Dep, GridSettings)
%
% STEMMUS indexes layers bottom-to-top (1=deepest), opposite of MODFLOW.
%
% Inputs:
%   soilThick    - cumulative layer thickness top-to-bottom [cm], NN×1
%   gw_Dep       - depth to water table [cm]
%   GridSettings - struct with fields .NL and .NN (use cfg.Grid)
%
% Outputs:
%   indxBotmLayer   - layer index (bottom-to-top, STEMMUS convention)
%   indxBotmLayer_R - layer index (top-to-bottom, MODFLOW convention)
%
% Reference: Mostafa et al. (2022) GMD preprint, MODFLOW-STEMMUS BMI coupling

indxBotmLayer_R = [];

for i = 1:GridSettings.NL
    midThick = (soilThick(i) + soilThick(i + 1)) / 2;
    if gw_Dep >= soilThick(i) && gw_Dep < soilThick(i + 1)
        if gw_Dep < midThick
            indxBotmLayer_R = i;
        else
            indxBotmLayer_R = i + 1;
        end
        break
    elseif gw_Dep >= soilThick(i + 1)
        continue
    end
end

indxBotmLayer = GridSettings.NN - indxBotmLayer_R + 1;
end
