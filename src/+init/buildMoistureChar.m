function soil = buildMoistureChar(soil)
%BUILDMOISTURECHAR  Stage B2: characteristic water contents from the VG curve.
%
%   soil = init.buildMoistureChar(soil) augments the per-layer soil parameter
%   struct (from init.buildSoilParameters) with the characteristic water
%   contents that the legacy SoilProperties.calculate_characteristic_
%   moisture_contents derives from the van Genuchten retention curve at three
%   fixed matric heads. Computed per layer directly from the layer VG params,
%   which is identical to the legacy "per-type then distribute" because the VG
%   params are already the distributed, bottom-up per-layer values.
%
%   Inputs:
%     soil   struct  built.parameters.soil with [NL x 1] theta_s/theta_r/alpha/n
%
%   Returns soil with these added [NL x 1] fields (bottom-up, [m3/m3]):
%     field_capacity  theta at -0.01 MPa (h = -101.9716 cm)
%     wilting_point   theta at -1.5  MPa (h = -30591 cm)
%     air_dry         theta at -10   MPa (h = -101970 cm)
%     xk              theta_r + 0.02  (air-water interface water content)
%
%   No thermal correction (T = 0, enable_thermal = 0), matching the legacy.
%   Reference: legacy SoilProperties.calculate_characteristic_moisture_contents
%   (matric heads -0.01 / -1.5 / -10 MPa; XK = theta_r + 0.02).
% Editor: Jizixing Zhu; Date: 2026-06-27.

    h_field = -101.9716;    % [cm]  -0.01 MPa  field capacity
    h_wilt  = -30591;       % [cm]  -1.5  MPa  wilting point
    h_air   = -101970;      % [cm]  -10   MPa  air-dry

    theta_s = soil.theta_s(:);
    theta_r = soil.theta_r(:);
    alpha   = soil.alpha(:);
    n       = soil.n(:);
    NL      = numel(theta_s);

    field_capacity = zeros(NL, 1);
    wilting_point  = zeros(NL, 1);
    air_dry        = zeros(NL, 1);
    for i = 1:NL
        field_capacity(i) = equations.van_genuchten_h_to_theta(h_field, ...
            theta_s(i), theta_r(i), alpha(i), n(i), 0, 0);
        wilting_point(i)  = equations.van_genuchten_h_to_theta(h_wilt, ...
            theta_s(i), theta_r(i), alpha(i), n(i), 0, 0);
        air_dry(i)        = equations.van_genuchten_h_to_theta(h_air, ...
            theta_s(i), theta_r(i), alpha(i), n(i), 0, 0);
    end

    soil.field_capacity = field_capacity;
    soil.wilting_point  = wilting_point;
    soil.air_dry        = air_dry;
    soil.xk             = theta_r + 0.02;
end
