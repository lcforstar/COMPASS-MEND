function [U_nodal, U_T, U_h, U_P, V] = assembleNodalInternalEnergy(U_endpoint, CTh, CTT, CTa, dz)
%ASSEMBLENODALINTERNALENERGY  Aggregate endpoint internal energy to FEM nodes.
%
%   [U_nodal, U_T, U_h, U_P, V] = soil.thermal.assembleNodalInternalEnergy(...
%       U_endpoint, CTh, CTT, CTa, dz)
%   averages the element-endpoint internal-energy state and its derivatives
%   onto FEM nodes using control-volume weights (§5.9.2, Eqs. 5.29a–5.29c).
%
%   COMPASS node ordering is bottom-up: node 1 = deepest (bottom boundary),
%   node NN = surface.  DeltZ(1) = deepest layer, DeltZ(NL) = top layer.
%   Element e spans nodes (e, e+1); its left endpoint (col 1) is the lower
%   node, right endpoint (col 2) is the upper node.
%
%   The control-volume weight for node i is:
%     V(1)       = dz(1)/2
%     V(i)       = dz(i-1)/2 + dz(i)/2   for interior nodes
%     V(NN)      = dz(NL)/2
%
%   Inputs:
%     U_endpoint  [J/cm³]       endpoint internal energy (NL×2)
%     CTh         [J·cm⁻³·cm⁻¹] ∂U/∂h at endpoints (NL×2)
%     CTT         [J·cm⁻³·K⁻¹]  ∂U/∂T at endpoints (NL×2)
%     CTa         [J·cm⁻³·uP⁻¹] ∂U/∂P_g at endpoints (NL×2)
%     dz          [cm]          element thicknesses, bottom-up (NL×1)
%
%   Returns:
%     U_nodal  [J/cm³]       control-volume nodal internal energy (NN×1)
%     U_T      [J·cm⁻³·K⁻¹]  ∂U_nodal/∂T  (NN×1)
%     U_h      [J·cm⁻³·cm⁻¹] ∂U_nodal/∂h  (NN×1)
%     U_P      [J·cm⁻³·uP⁻¹] ∂U_nodal/∂P_g (NN×1)
%     V        [cm]          control-volume weights (NN×1)
%
%   Reference: docs/compass_soil_water_heat_physics.md §5.9.2 Eqs. 5.29a–5.29c.
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%
%   See also: soil.thermal.evaluateInternalEnergyEndpoints

    arguments
        U_endpoint (:,:) double
        CTh        (:,:) double
        CTT        (:,:) double
        CTa        (:,:) double
        dz         (:,1) double
    end

    NL = size(U_endpoint, 1);
    NN = NL + 1;
    assert(size(U_endpoint, 2) == 2, ...
        'assembleNodalInternalEnergy:badShape', ...
        'U_endpoint must be NL×2.');
    assert(isequal(size(CTh), [NL,2]) && isequal(size(CTT), [NL,2]) ...
        && isequal(size(CTa), [NL,2]), ...
        'assembleNodalInternalEnergy:badShape', ...
        'CTh, CTT, CTa must all be NL×2.');
    assert(all(isfinite(U_endpoint(:))) && all(isfinite(CTh(:))) ...
        && all(isfinite(CTT(:))) && all(isfinite(CTa(:))), ...
        'assembleNodalInternalEnergy:nonFinite', ...
        'U_endpoint, CTh, CTT, CTa must all be finite.');
    assert(numel(dz) == NL && all(isfinite(dz)) && all(dz > 0), ...
        'assembleNodalInternalEnergy:badDz', ...
        'dz must be NL×1 with finite positive values.');

    halfDz = dz / 2;

    % Control-volume weights (bottom-up).
    V = zeros(NN, 1);
    V(1)       = halfDz(1);                    % bottom node
    V(2:NL)    = halfDz(1:NL-1) + halfDz(2:NL); % interior
    V(NN)      = halfDz(NL);                    % surface node

    invV = 1 ./ V;

    U_nodal = zeros(NN, 1);
    U_T     = zeros(NN, 1);
    U_h     = zeros(NN, 1);
    U_P     = zeros(NN, 1);

    % Bottom node (i=1): only element 1, left endpoint (col 1, lower node).
    U_nodal(1) = halfDz(1) * U_endpoint(1,1) * invV(1);
    U_T(1)     = halfDz(1) * CTT(1,1)        * invV(1);
    U_h(1)     = halfDz(1) * CTh(1,1)        * invV(1);
    U_P(1)     = halfDz(1) * CTa(1,1)        * invV(1);

    % Interior nodes (i=2..NL): right endpoint of element i-1 (col 2, upper
    % node of element below) + left endpoint of element i (col 1, lower node
    % of element above).
    for i = 2:NL
        U_nodal(i) = (halfDz(i-1) * U_endpoint(i-1,2) ...
                    + halfDz(i)   * U_endpoint(i,1))   * invV(i);
        U_T(i)     = (halfDz(i-1) * CTT(i-1,2) ...
                    + halfDz(i)   * CTT(i,1))          * invV(i);
        U_h(i)     = (halfDz(i-1) * CTh(i-1,2) ...
                    + halfDz(i)   * CTh(i,1))          * invV(i);
        U_P(i)     = (halfDz(i-1) * CTa(i-1,2) ...
                    + halfDz(i)   * CTa(i,1))          * invV(i);
    end

    % Surface node (i=NN): only element NL, right endpoint (col 2, upper node).
    U_nodal(NN) = halfDz(NL) * U_endpoint(NL,2) * invV(NN);
    U_T(NN)     = halfDz(NL) * CTT(NL,2)        * invV(NN);
    U_h(NN)     = halfDz(NL) * CTh(NL,2)        * invV(NN);
    U_P(NN)     = halfDz(NL) * CTa(NL,2)        * invV(NN);
end
