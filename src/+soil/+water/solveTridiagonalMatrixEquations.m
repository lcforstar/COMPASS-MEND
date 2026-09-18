function [SoilVariables, CHKWATER] = solveTridiagonalMatrixEquations(WaterMatrix, RHS, SoilVariables, DSTMAX, hcrit, cfg, isRaining) %#ok<INUSD>
% solveTridiagonalMatrixEquations  Solve Richards equation tridiagonal system via Thomas algorithm
%   (Eq. 4.25, STEMMUS Technical Notes p.41), with zero-pivot guard and
%   groundwater coupling support.
%
% DSTMAX and isRaining remain in the compatibility signature; DSTMAX is unused.
% Editor: Jizixing Zhu; Date: 2026-07-12.
% Change: 2026-07-12 Jizixing Zhu -- retain positive heads when the transactional
% variably saturated solver sets cfg.AllowPositiveHead=true.
% Change: 2026-07-13 Jizixing Zhu -- delegate new transactional solves to
% solver.solveTridiagonal; retain the legacy compatibility sweep otherwise.
% Change: 2026-07-14 Jizixing Zhu -- share the absolute/correction head-update
% safeguards through soil.water.applyHeadUpdateSafeguards.

if nargin < 7 || isempty(isRaining)
    isRaining = false;
end

C4   = WaterMatrix.C4;
C4_a = WaterMatrix.C4_a;
hh   = SoilVariables.hh(:);
hh_prev = hh;
NN   = length(hh);

indxBotm  = cfg.Grid.indxBotmLayer;
GW_active = indxBotm > 1;
allowPositiveHead = isfield(cfg, 'AllowPositiveHead') && cfg.AllowPositiveHead;
my_eps = 1e-20;

if isfield(cfg, 'UseCommonThomas') && cfg.UseCommonThomas
    [RHS, C4(:, 1)] = solver.solveTridiagonal( ...
        C4_a(1:NN-1), C4(:, 1), C4(:, 2), RHS, [indxBotm, NN]);
else

% ── Forward elimination (start from indxBotm, skip GW-saturated layers) ─────
val = C4(indxBotm, 1);
s   = sign(val); if s == 0, s = 1; end
C4(indxBotm, 1) = max(abs(val), my_eps) * s;
RHS(indxBotm)   = RHS(indxBotm) / C4(indxBotm, 1);

for ML = indxBotm+1:NN
    val_ML = C4(ML, 1) - C4_a(ML-1) * C4(ML-1, 2) / C4(ML-1, 1);
    s_ML   = sign(val_ML); if s_ML == 0, s_ML = 1; end
    C4(ML, 1) = max(abs(val_ML), my_eps) * s_ML;
    RHS(ML)   = (RHS(ML) - C4_a(ML-1) * RHS(ML-1)) / C4(ML, 1);
end

% ── Back substitution ────────────────────────────────────────────────────
for ML = (NN-1):-1:indxBotm
    RHS(ML) = RHS(ML) - C4(ML, 2) * RHS(ML+1) / C4(ML, 1);
end
end

% ── GW coupling: fill saturated layers below indxBotm with hydrostatic head ──
if GW_active
    for ML = indxBotm:-1:2
        RHS(ML-1) = RHS(ML) + cfg.Grid.DeltZ(ML-1);
    end
end

[hh, CHKWATER] = soil.water.applyHeadUpdateSafeguards( ...
    hh_prev, RHS(1:NN), hcrit, isRaining, allowPositiveHead, indxBotm);

% The MODFLOW-coupled saturated segment is not a Richards/Picard unknown.
% applyHeadUpdateSafeguards may damp a large interface jump in every entry,
% including the hydrostatic continuation formed above. Reconstruct that
% continuation *after* damping so the externally imposed interface head and
% the gravity relation h(i-1)=h(i)+dz(i-1) remain exact.
if GW_active
    for ML = indxBotm:-1:2
        hh(ML-1) = hh(ML) + cfg.Grid.DeltZ(ML-1);
    end
end

SoilVariables.hh = hh;
end
