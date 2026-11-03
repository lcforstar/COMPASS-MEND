function domain = buildDomain(gridCfg)
%BUILDDOMAIN  Stage A: construct the 1-D vertical soil grid from [grid].
%
%   domain = init.buildDomain(gridCfg) builds the vertical discretization,
%   reproducing the legacy non-uniform spacing table through
%   init.buildVerticalGrid. It owns only domain geometry; soil horizon mapping,
%   root distribution and initial-condition interpolation happen in later
%   builders.
%
%   Inputs:
%     gridCfg  struct  the normalized grid config (built.grid from init.buildGrid):
%                      Tot_Depth [cm], NL/requested layer count [-],
%                      lEqlspace (0/1), max_elements [-].
%
%   Conventions (MUST match the legacy main loop):
%     - depth is positive downward; surface = 0, column bottom = tot_depth.
%     - NODE ordering is BOTTOM-UP: node 1 = column bottom (deepest), node NN =
%       surface. (Legacy state arrays use this after InitialProfile's flip.)
%     - dz is LAYER thickness, also bottom-up: dz(1) = deepest layer (coarse),
%       dz(NL) = top layer (1 cm in the current default depth schedule).
%
%   Returns domain with fields:
%     dz        [cm,  NL x 1]  layer thickness, bottom-up
%     NL        [-]            actual number of layers (non-uniform grid may
%                              differ from the requested gridCfg.NL)
%     NN        [-]            number of nodes = NL + 1
%     z_node    [cm,  NN x 1]  node depth below surface, bottom-up
%                              (z_node(1) = tot_depth, z_node(NN) = 0)
%     sumdz     [cm,  NL x 1]  cumsum(flip(dz)), matches legacy GridSettings.SUMDELTZ
%     tot_depth [cm]           total column depth
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    tot_depth = gridCfg.Tot_Depth;

    if isfield(gridCfg, 'lEqlspace') && gridCfg.lEqlspace
        % uniform grid: NL equal layers of tot_depth/NL
        NL = gridCfg.NL;
        dz = (tot_depth / NL) * ones(NL, 1);   % [cm, NL x 1]
    else
        % non-uniform grid (dense near surface): legacy spacing table.
        maxElements = 5000;
        if isfield(gridCfg, 'max_elements')
            maxElements = gridCfg.max_elements;
        end
        [dz, ~, NL, ~] = init.buildVerticalGrid(gridCfg.NL, tot_depth, maxElements);   % [cm, NL x 1] bottom-up
    end

    NN = NL + 1;

    % node depths below surface, bottom-up (node 1 = bottom = tot_depth).
    % Top-down cumulative first (0 at surface), then flip to bottom-up so it
    % matches the legacy state-array ordering. dz(NL) is the top (thinnest) layer.
    depth_top_down = zeros(NN, 1);             % [cm], index 1 = surface (0)
    acc = 0;
    for k = 1:NL
        acc = acc + dz(NL + 1 - k);            % add top layer first, go downward
        depth_top_down(k + 1) = acc;
    end
    z_node = flip(depth_top_down);             % bottom-up: z_node(1)=tot_depth, z_node(NN)=0

    domain = struct();
    domain.dz        = dz(:);
    domain.NL        = NL;
    domain.NN        = NN;
    domain.z_node    = z_node;
    domain.sumdz     = cumsum(flip(dz(:)));    % matches legacy GridSettings.SUMDELTZ
    domain.tot_depth = tot_depth;
end
