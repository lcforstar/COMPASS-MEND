function grid = buildGrid(rawConfig)
%BUILDGRID  Normalize [grid] scalars (the geometry itself is init.buildDomain).
%
%   grid = init.buildGrid(rawConfig) maps the friendly site.toml grid names onto
%   the legacy grid scalars (Tot_Depth, NL, NS, ...). This produces only the
%   grid *configuration*; the discretized node geometry is built separately by
%   init.buildDomain(grid). Pure config translation -- no physics.
%
%   Inputs:
%     rawConfig  struct  raw.config (uses .grid)
%
%   Returns grid  struct  legacy grid scalars:
%     Tot_Depth [cm], NL (requested layer count), NS/nS/nD, max_elements,
%     lEqlspace (1 = uniform spacing).
%
%   Reference: legacy GridSettings (field names).
% Editor: Jizixing Zhu; Date: 2026-06-27.

    g = rawConfig.grid;
    grid = struct();
    grid.Tot_Depth    = g.tot_depth;
    grid.NL           = g.NL;            % requested layer count (actual NL is built.domain.NL)
    grid.NS           = g.NS;
    grid.nS           = g.nS;
    grid.nD           = g.nD;
    grid.max_elements = g.max_elements;
    grid.lEqlspace    = double(logical(g.uniform_grid));
end
