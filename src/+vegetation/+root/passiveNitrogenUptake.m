function [SinkS, SPUptake, uptakeDiag] = passiveNitrogenUptake( ...
        Srt, DeltZ, C_NH4, C_NO3)
%PASSIVENITROGENUPTAKE  Water-linked passive NH4/NO3 mass-flow uptake.
%
%   [SinkS,SPUptake,uptakeDiag] = vegetation.root.passiveNitrogenUptake(...)
%   applies the original passive formulation at every layer endpoint:
%
%     SinkS_species = max(Srt,0) * max(C_species,0)
%
%   Hydraulic redistribution does not return nitrogen to soil because only
%   positive root-water extraction contributes. There is no plant-demand or
%   empirical concentration cap in this fallback.
%
%   Inputs:
%     Srt       [in] [1/s, NLx2]       root-water uptake at endpoints
%     DeltZ     [in] [cm, NLx1]        layer thickness
%     C_NH4     [in] [ug/cm3, NNx1]    nodal ammonium concentration
%     C_NO3     [in] [ug/cm3, NNx1]    nodal nitrate concentration
%
%   Returns:
%     SinkS     [out] [ug/cm3/s, NLx2x2] endpoint NH4/NO3 sinks
%     SPUptake  [out] [ug/cm2/s, 1x2]    depth-integrated uptake
%     uptakeDiag [out] struct             passive-mode audit values
%
%   Reference: pre-M1 applyRootUptake passive branch and
%   docs/root_water_and_n_uptake_method.md, Part 5.
% Editor: Jizixing Zhu; Date: 2026-07-11.

    NL = numel(DeltZ);
    assert(isequal(size(Srt), [NL, 2]), ...
        'passiveNitrogenUptake:badSrtShape', 'Srt must be NLx2.');
    assert(numel(C_NH4) >= NL + 1 && numel(C_NO3) >= NL + 1, ...
        'passiveNitrogenUptake:badConcentrationShape', ...
        'NH4 and NO3 concentrations must contain at least NL+1 nodes.');

    waterSink = max(Srt, 0);
    nh4Nodes = max(C_NH4(:), 0);
    no3Nodes = max(C_NO3(:), 0);
    SinkS = zeros(NL, 2, 2);
    SinkS(:, 1, 1) = waterSink(:, 1) .* nh4Nodes(1:NL);
    SinkS(:, 2, 1) = waterSink(:, 2) .* nh4Nodes(2:NL+1);
    SinkS(:, 1, 2) = waterSink(:, 1) .* no3Nodes(1:NL);
    SinkS(:, 2, 2) = waterSink(:, 2) .* no3Nodes(2:NL+1);

    dz = DeltZ(:);
    SPUptake = [ ...
        sum(0.5 * (SinkS(:, 1, 1) + SinkS(:, 2, 1)) .* dz), ...
        sum(0.5 * (SinkS(:, 1, 2) + SinkS(:, 2, 2)) .* dz)];
    uptakeDiag = struct('scheme', 'passive_mass_flow', ...
        'passive_nh4', SPUptake(1), 'passive_no3', SPUptake(2), ...
        'active_total', 0);
end
