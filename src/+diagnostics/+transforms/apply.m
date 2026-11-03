function out = apply(value, transform_name, NL, NN)
%APPLY  Apply a pre-write shape/orientation transform to a diagnostic value.
%
% COMPASS internal representations are mixed:
%   - some fields are NN-node-vector (hh, TT, P_g)
%   - some are NL×nD element matrix (Theta_LL, K_h)
%   - some come in reversed depth order (bottom-up)
% NetCDF output wants a uniform layout per variable: usually (tS, NN) with
% node 1 = surface. This function does the conversion at write-time.
%
% Inputs:
%   value          [in]   value as held in COMPASS memory (any shape)
%   transform_name [in]   'none' | 'edge_to_node' | 'edge_to_node_reversed_layer'
%                         | 'flip_node' | 'scalar'
%   NL             [in]   layer count
%   NN             [in]   node count
%
% Returns:
%   out            transformed value ready for buffer assignment

    switch transform_name
        case 'none'
            out = value;
        case 'scalar'
            out = double(value(1));
        case 'edge_to_node'
            % NL×2 element matrix → NN-vector (node value = average of two endpoints)
            % Endpoint convention: matrix(ML, 1) = top of element, (ML, 2) = bottom.
            % Node 1 = top of element 1; node ML = bottom of (ML-1) ≡ top of ML; node NN = bottom of NL.
            if isequal(size(value), [NL, 2])
                out = zeros(NN, 1);
                out(1) = value(1, 1);
                out(NN) = value(NL, 2);
                for ML = 2:NL
                    out(ML) = 0.5 * (value(ML-1, 2) + value(ML, 1));
                end
            elseif isvector(value) && numel(value) == NN
                out = value(:);
            else
                % unknown shape, pass through
                out = value;
            end
        case 'edge_to_node_reversed_layer'
            % Convert bottom-up element endpoints to bottom-up nodes first,
            % then reverse the node vector to the public surface-to-bottom
            % output order. Reversing rows before endpoint assembly is wrong
            % because it does not also exchange the endpoint columns.
            if isequal(size(value), [NL, 2])
                nodeBottomUp = zeros(NN, 1);
                nodeBottomUp(1) = value(1, 1);
                nodeBottomUp(NN) = value(NL, 2);
                for ML = 2:NL
                    nodeBottomUp(ML) = ...
                        0.5 * (value(ML-1, 2) + value(ML, 1));
                end
                out = flip(nodeBottomUp);
            elseif isvector(value) && numel(value) == NN
                out = flip(value(:));
            else
                out = value;
            end
        case 'flip_node'
            % some COMPASS fields are stored bottom→top (e.g. flux_vars from QL),
            % flip so node 1 = surface
            out = flip(value(:));
        otherwise
            error('diagnostics:transforms:Unknown', 'unknown transform: %s', transform_name);
    end
end
