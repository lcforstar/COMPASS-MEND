function y_interp = interpLinearWithExtrapolation(x_data, y_data, x_nodes)
%INTERPLINEARWITHEXTRAPOLATION  Linear interpolation with left/right extrapolation.
%
%   y_interp = interpLinearWithExtrapolation(x_data, y_data, x_nodes)
%
%   Inputs:
%     x_data   original x coordinates (row or column vector)
%     y_data   original y values (same size as x_data)
%     x_nodes  target x coordinates (row or column vector)
%
%   Returns:
%     y_interp interpolated values (same size as x_nodes)
%
%   For each point in x_nodes: linear interpolation inside the data range;
%   linear extrapolation using the boundary slope outside it.
%
%   Editor: Jizixing Zhu; Date: 2026-02-27.
%   Change: 2026-06-27 Jizixing Zhu -- doc comments translated to English
%     (kept public init.X: tests/test_buildState calls it for gate-A, so it
%     cannot move to +init/private).

    x_data = x_data(:);
    y_data = y_data(:);
    x_nodes = x_nodes(:);

    n = length(x_data);
    if n > 1
        slopes = (y_data(2:n) - y_data(1:n-1)) ./ (x_data(2:n) - x_data(1:n-1));
    else
        slopes = [];
    end

    y_interp = zeros(size(x_nodes));

    for k = 1:length(x_nodes)
        x_val = x_nodes(k);

        if x_val >= x_data(1) && x_val <= x_data(end)
            y_interp(k) = linear_interp_point(x_data, y_data, slopes, x_val);
        elseif x_val < x_data(1)
            if n > 1
                y_interp(k) = y_data(1) + slopes(1) * (x_val - x_data(1));
            else
                y_interp(k) = y_data(1);
            end
        else
            if n > 1
                y_interp(k) = y_data(n) + slopes(n-1) * (x_val - x_data(n));
            else
                y_interp(k) = y_data(n);
            end
        end
    end
end

function y_val = linear_interp_point(x_data, y_data, slopes, x_val)
% Editor: Jizixing Zhu; Date: 2026-06-27.
    n = length(x_data);

    if abs(x_val - x_data(n)) < 1e-10
        y_val = y_data(n);
        return;
    end

    for j = 1:n - 1
        if x_val >= x_data(j) && x_val < x_data(j + 1)
            y_val = y_data(j) + slopes(j) * (x_val - x_data(j));
            return;
        end
    end

    y_val = y_data(n);
end
