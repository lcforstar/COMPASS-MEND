function groundWaterDepth = calculateGroundWaterDepth(topLevel, headBotmLayer, Tot_Depth)
% calculateGroundWaterDepth - Compute depth from land surface to water table
%
%   groundWaterDepth = calculateGroundWaterDepth(topLevel, headBotmLayer, Tot_Depth)
%
% Inputs:
%   topLevel      - elevation of top aquifer surface [cm], from MODFLOW via BMI
%   headBotmLayer - groundwater head at bottom layer [cm], from MODFLOW via BMI
%   Tot_Depth     - total soil column depth [cm]
%
% Output:
%   groundWaterDepth - depth from land surface to water table [cm]
%
% Reference: Mostafa et al. (2022) GMD preprint, MODFLOW-STEMMUS BMI coupling

groundWaterDepth = topLevel - headBotmLayer;

if groundWaterDepth <= 0
    warning('The soil is fully saturated up to the land surface level!');
    groundWaterDepth = 1.0; % minimum value to avoid model crash
end
end
