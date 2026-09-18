function storage_mm = calculateSoilWaterStorageMm(theta_nodes, DeltZ)
% calculateSoilWaterStorageMm Integrate nodal volumetric water content.
%
% theta_nodes is bottom-to-top with NN nodes. DeltZ is bottom-to-top with
% NL elements in cm. The result is equivalent water depth in mm.

theta_nodes = theta_nodes(:);
DeltZ = DeltZ(:);

n_elem = min(numel(DeltZ), numel(theta_nodes) - 1);
if n_elem <= 0
    storage_mm = NaN;
    return
end

theta_left = theta_nodes(1:n_elem);
theta_right = theta_nodes(2:n_elem + 1);
theta_elem = 0.5 * (theta_left + theta_right);
storage_cm = sum(theta_elem .* DeltZ(1:n_elem), 'omitnan');
storage_mm = storage_cm * 10;
end
