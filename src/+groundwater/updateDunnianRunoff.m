function [runoffDunnian, update_Precip_msr] = updateDunnianRunoff(Precip_msr, groundWaterDepth)
% updateDunnianRunoff - Compute Dunnian runoff when soil is fully saturated
%
%   [runoffDunnian, update_Precip_msr] = updateDunnianRunoff(Precip_msr, groundWaterDepth)
%
% When gw_Dep <= 1 cm (soil fully saturated), all precipitation becomes
% Dunnian runoff and net infiltration is set to zero.
%
% Inputs:
%   Precip_msr      - precipitation time series [mm], full array
%   groundWaterDepth - depth to water table [cm]
%
% Outputs:
%   runoffDunnian    - Dunnian runoff [mm], same size as Precip_msr
%   update_Precip_msr - modified precipitation after subtracting Dunnian runoff
%
% Reference: Mostafa et al. (2022) GMD preprint, MODFLOW-STEMMUS BMI coupling

runoffDunnian     = zeros(size(Precip_msr));
update_Precip_msr = Precip_msr;

if groundWaterDepth <= 1.0
    runoffDunnian     = Precip_msr;
    update_Precip_msr = zeros(size(Precip_msr));
end
end
