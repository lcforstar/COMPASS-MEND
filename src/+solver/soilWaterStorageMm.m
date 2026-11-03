function mm = soilWaterStorageMm(Theta_layers, DeltZ)
%SOILWATERSTORAGEMM  Column soil-water storage [mm] from the [NL x 2] layer-endpoint theta (WBC-3).
%
%   mm = solver.soilWaterStorageMm(Theta_layers, DeltZ) builds the bottom-up NN-node theta vector the
%   legacy storage integrator expects -- [Theta(1:NL,1); Theta(NL,2)] (COMPASS.m:2130) -- and calls
%   solver.calculateSoilWaterStorageMm. Used for the WBC-3 storage-change term on the accepted
%   Y.Theta_L (liquid water layer state; under water-only no-freeze this is the total-water storage).
%
%   Inputs:
%     Theta_layers [in]  [m3/m3, NL x 2] bottom-up layer endpoints (col1 = lower node, col2 = upper).
%     DeltZ        [in]  [cm, NL x 1] layer thickness, bottom-up.
%
%   Returns:
%     mm  [mm] depth-integrated column water storage.
%
%   Reference: legacy node construction COMPASS.m:2130; integrator solver.calculateSoilWaterStorageMm.
%   Editor: Jizixing Zhu; Date: 2026-07-08.
    NL = size(Theta_layers, 1);
    node = [Theta_layers(1:NL, 1); Theta_layers(NL, 2)];   % [m3/m3, NN] bottom-up node theta
    mm = solver.calculateSoilWaterStorageMm(node, DeltZ);
end
