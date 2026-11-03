function [RAI, rootProfile] = calculateRootDistribution(cfg, plantHydro, mode)
% calculateRootDistribution - Root distribution profile for plant hydraulics
%
% Inputs:
%   cfg        - Model settings with Grid.NL, Grid.DeltZ, Grid.Tot_Depth
%   plantHydro - Plant hydraulic struct containing baseline RAI and Lmax
%   mode       - 'jackson'  : use the baseline Jackson-style profile
%                'measured' : read measured root length density from CSV
%
% Output:
%   RAI        - NLx1 layer RAI profile
%   rootProfile - diagnostic/source metadata. For measured mode this includes
%                 layer root length density [m/m3].
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-16 Jizixing Zhu — read the soil-column and rooting depths
%   from plantHydro in the measured-profile branch.

rootProfile = struct('mode', char(mode), ...
    'root_length_density_m_m3', [], ...
    'layer_top_cm', [], ...
    'layer_bottom_cm', [], ...
    'layer_center_depth_cm', [], ...
    'measured_file', '', ...
    'total_root_length_m_m2', NaN);

NL = cfg.Grid.NL;
DeltZ = cfg.Grid.DeltZ;
RL = plantHydro.RL;
Lmax = plantHydro.Lmax;

if nargin < 3 || isempty(mode)
    mode = 'jackson';
end

switch mode
    case 'jackson'
        % Use the baseline profile prepared in Init_Plant.
        RAI = plantHydro.RAI(:);

    case 'measured'
        measured_file = '';
        if isprop(cfg.Veg, 'root_measured_file')
            measured_file = strtrim(char(cfg.Veg.root_measured_file));
        end
        if isempty(measured_file)
            error('calculateRootDistribution:MissingMeasuredFile', ...
                'root_dist_mode=1 requires cfg.Veg.root_measured_file.');
        end
        if ~isfile(measured_file)
            measured_file = fullfile(pwd, measured_file);
        end
        if ~isfile(measured_file)
            error('calculateRootDistribution:MeasuredFileNotFound', ...
                'Measured root file not found: %s', measured_file);
        end
        rootProfile.measured_file = measured_file;

        T = readtable(measured_file, 'VariableNamingRule', 'preserve');
        names = matlab.lang.makeValidName(T.Properties.VariableNames);
        T.Properties.VariableNames = names;

        depth_col = find(strcmpi(names, 'depth_cm'), 1);
        if isempty(depth_col)
            error('calculateRootDistribution:MissingDepthColumn', ...
                'Measured root CSV must contain a depth_cm column.');
        end

        rld_col = find(strcmpi(names, 'root_length_density_m_m3'), 1);
        rld_unit = 'm_m3';
        if isempty(rld_col)
            rld_col = find(strcmpi(names, 'root_length_density_cm_cm3'), 1);
            rld_unit = 'cm_cm3';
        end
        if isempty(rld_col)
            error('calculateRootDistribution:MissingRLDColumn', ...
                ['Measured root CSV must contain root_length_density_m_m3 ', ...
                 'or root_length_density_cm_cm3.']);
        end

        depth_cm = T{:, depth_col};
        rld = T{:, rld_col};
        valid = isfinite(depth_cm) & isfinite(rld) & depth_cm >= 0 & rld >= 0;
        depth_cm = depth_cm(valid);
        rld = rld(valid);
        if isempty(depth_cm)
            error('calculateRootDistribution:EmptyMeasuredRootData', ...
                'Measured root CSV contains no valid non-negative rows.');
        end
        if strcmp(rld_unit, 'cm_cm3')
            rld = rld * 10000;  % cm/cm3 -> m/m3
        end

        [depth_cm, order] = sort(depth_cm(:));
        rld = rld(order);
        [depth_cm, ia] = unique(depth_cm, 'stable');
        rld = rld(ia);

        measured_top_cm = [0; depth_cm(1:end-1)];
        measured_bottom_cm = depth_cm;

        layer_top_cm = zeros(NL, 1);
        layer_bottom_cm = zeros(NL, 1);
        layer_center_depth_cm = zeros(NL, 1);
        depth_from_bottom_cm = 0;
        for ML = 1:NL
            layer_top_cm(ML) = RL - (depth_from_bottom_cm + DeltZ(ML));
            layer_bottom_cm(ML) = RL - depth_from_bottom_cm;
            layer_center_depth_cm(ML) = 0.5 * (layer_top_cm(ML) + layer_bottom_cm(ML));
            depth_from_bottom_cm = depth_from_bottom_cm + DeltZ(ML);
        end

        rld_layer = zeros(NL, 1);
        for ML = 1:NL
            layer_top = max(layer_top_cm(ML), 0);
            layer_bottom = min(layer_bottom_cm(ML), Lmax);
            if layer_bottom <= layer_top
                continue;
            end

            overlap_sum = 0;
            weighted_sum = 0;
            for k = 1:numel(rld)
                overlap = min(layer_bottom, measured_bottom_cm(k)) - ...
                          max(layer_top, measured_top_cm(k));
                if overlap > 0
                    weighted_sum = weighted_sum + rld(k) * overlap;
                    overlap_sum = overlap_sum + overlap;
                end
            end
            if overlap_sum > 0
                rld_layer(ML) = weighted_sum / overlap_sum;
            else
                rld_layer(ML) = interp1(depth_cm, rld, layer_center_depth_cm(ML), 'linear', 0);
            end
        end
        rld_layer = max(rld_layer, 0);

        r_root_m = 5e-4;
        if isprop(cfg.Veg, 'root_rroot_mm') && cfg.Veg.root_rroot_mm > 0
            r_root_m = cfg.Veg.root_rroot_mm * 1e-3;
        end
        DeltZ_m = DeltZ(:) / 100;
        RAI = 2 * pi * r_root_m * rld_layer .* DeltZ_m;
        if sum(RAI) <= 0
            error('calculateRootDistribution:ZeroMeasuredRAI', ...
                'Measured root profile produced zero RAI on the model grid.');
        end
        rootProfile.root_length_density_m_m3 = rld_layer(:);
        rootProfile.layer_top_cm = layer_top_cm(:);
        rootProfile.layer_bottom_cm = layer_bottom_cm(:);
        rootProfile.layer_center_depth_cm = layer_center_depth_cm(:);
        rootProfile.total_root_length_m_m2 = sum(rld_layer(:) .* DeltZ_m(:));

    otherwise
        error('calculateRootDistribution:UnknownMode', ...
            'Unknown mode "%s". Use ''jackson'' or ''measured''.', mode);
end

RAI = RAI(:);
end
