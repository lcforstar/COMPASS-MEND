function matrices = assembleEquationMatrices(equations, dz, rootSink)
%ASSEMBLEEQUATIONMATRICES  Dispatch process-specific soil FEM matrix assembly.
%
%   matrices = soil.assembleEquationMatrices(equations, dz, rootSink)
%   converts the physical coefficient banks from soil.buildEquationCoefficients
%   into the legacy-compatible water, heat, and dry-air matrix banks. Each
%   enabled equation retains its original one-pass element assembly. Empty
%   heat or dry-air banks remain empty. No RHS, boundary condition, linear
%   solve, or state mutation occurs here.
%
%   Inputs:
%     equations [struct] coefficient container with .water/.heat/.dry_air
%     dz        [cm, NL x 1] element thickness, bottom-up
%     rootSink  [1/s, NL x 2] root-water sink for the water load vector
%
%   Returns:
%     matrices.water   [struct] C1/C2/C4/C5/C6/C7/C9 matrix bank
%     matrices.heat    [struct] C1-C7 heat matrix bank, or empty
%     matrices.dry_air [struct] C1-C7 dry-air matrix bank, or empty
%
%   Reference:
%     soil.water.assembleWaterEquationMatrices
%     soil.thermal.assembleThermalEquationMatrices
%     soil.gas.assembleDryAirEquationMatrices
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- restore the process-specific one-pass
%           assemblers and preserve the heat pressure-column switch.

    NL = numel(dz);
    if nargin < 3 || isempty(rootSink)
        rootSink = zeros(NL, 2);
    end
    assert(isstruct(equations) && isfield(equations, 'water') ...
        && isfield(equations, 'heat') && isfield(equations, 'dry_air'), ...
        'assembleEquationMatrices:badEquations', ...
        'equations must contain water, heat, and dry_air coefficient banks.');

    airEnabled = ~isempty(fieldnames(equations.dry_air));
    cfg = struct();
    cfg.Grid = struct('NL', NL, 'NN', NL + 1, 'DeltZ', dz(:));
    cfg.lSoilairefc = airEnabled;

    matrices = struct();
    matrices.water = soil.water.assembleWaterEquationMatrices( ...
        equations.water, rootSink, cfg);
    matrices.heat = struct();
    if ~isempty(fieldnames(equations.heat))
        matrices.heat = soil.thermal.assembleThermalEquationMatrices( ...
            equations.heat, cfg);
    end
    matrices.dry_air = struct();
    if airEnabled
        matrices.dry_air = soil.gas.assembleDryAirEquationMatrices( ...
            equations.dry_air, cfg);
    end
end
