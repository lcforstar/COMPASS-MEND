function node = mapLayerToNodeMean(layer, NN, NL)
% Map layer-centered values to nodes by legacy simple averaging:
% node 1 = layer 1, interior node = adjacent-layer mean, node NN = layer NL.
% (Distinct from equations.interpolateElementToNode, which is thickness-weighted.)
% Editor: Jizixing Zhu; Date: 2026-06-27.
    node = zeros(NN, 1);
    node(1) = layer(1);
    for i = 2:NL
        node(i) = (layer(i-1) + layer(i)) / 2;
    end
    node(NN) = layer(NL);
end
