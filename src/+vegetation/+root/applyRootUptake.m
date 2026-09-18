function [Srt, SinkS, SPUptake, uptakeDiag] = applyRootUptake( ...
        Srt_in, NL, DeltZ, RL, Lmax, Nefc, C_NH4, C_NO3, SoilState, RLD, ...
        daisyParams, nacryParams)
%APPLYROOTUPTAKE  Plant→soil N uptake dispatch (Phase C M1 coupled DAISY + Nacry).
%
% Integration point between plant water uptake (Srt) and soil solute transport
% (NH4/NO3 sinks). Normalizes Srt shape, then delegates N uptake computation
% to soil.nitrogen.root_uptake.nitrogenUptake which implements the coupled
% DAISY single-root + Nacry membrane M-M physics.
%
% Inputs:
%   Srt_in       [in]  [1/s]    Root water uptake — three accepted shapes:
%                                NN×1 node-based: maps Srt(ML,:) = [RWU(ML), RWU(ML+1)]
%                                NL×1 layer-based: duplicated to both nodes
%                                NL×2 layer-node form: used directly
%   NL           [in]  [-]      Number of soil layers.
%   DeltZ        [in]  [cm, NL×1]  Layer thickness.
%   RL           [in]  [cm]     Reserved (legacy interface compatibility).
%   Lmax         [in]  [cm]     Reserved (legacy interface compatibility).
%   Nefc         [in]  [0/1]    N cycling switch; 0 → returns zero sinks.
%   C_NH4        [in]  [μg/cm³, NN×1]  NH4 concentration at nodes.
%   C_NO3        [in]  [μg/cm³, NN×1]  NO3 concentration at nodes.
%   SoilState    [in]  struct   Required by nitrogenUptake (see its doc).
%                                Fields: Theta_LL, Theta_s, RHO_bulk, RHOKG,
%                                        DIF1111, TT.
%   RLD          [in]  [m-root/m³-soil, NL×1]  Root length density per layer.
%   daisyParams  [in]  struct   r_root_m, SRL.
%   nacryParams  [in]  struct   Vmax/Km/Cmin per species + solver tol/max_iter.
%
% Returns:
%   Srt          [NL×2, 1/s]   Pass-through layer-node RWU profile.
%   SinkS        [NL×2×2]      N sink terms (NH4/NO3); zero if Nefc=0.
%   SPUptake     [1×2]         Depth-integrated uptake [NH4, NO3].
%   uptakeDiag   struct        Per-layer DAISY+Nacry diagnostics, see nitrogenUptake.
%
% Reference:
%   cn_coupling_path_and_literature_audit.md §12.8 (Phase C M1 implementation plan)
%
% Phase C M1 (2026-05-18): signature simplified — removed lActRSU, wt, N_demand
% (M1 always uses coupled DAISY+Nacry; M2 will modulate C_min(t) inside nacryParams).

% ---- 1. Normalize Srt shape ----
if isvector(Srt_in)
    sv = Srt_in(:);
    if numel(sv) == NL + 1
        Srt = [sv(1:NL), sv(2:NL+1)];   % NN×1 node form
    elseif numel(sv) == NL
        Srt = [sv, sv];                  % NL×1 layer form duplicated
    else
        error('applyRootUptake:BadSrtSize', ...
            'Vector Srt_in must have length NL or NL+1.');
    end
else
    if size(Srt_in, 1) < NL || size(Srt_in, 2) < 2
        error('applyRootUptake:BadSrtShape', ...
            'Matrix Srt_in must be at least NLx2.');
    end
    Srt = Srt_in(1:NL, 1:2);
end

% RL/Lmax kept in signature for legacy callers; not used here.
unused_interface = RL + Lmax; %#ok<NASGU>

% ---- 2. N cycling off → return zero sinks ----
SinkS    = zeros(NL, 2, 2);
SPUptake = zeros(1, 2);
uptakeDiag = struct();

if Nefc ~= 1
    return;
end

% ---- 3. Delegate to coupled DAISY+Nacry physics ----
[SinkS, SPUptake, uptakeDiag] = soil.nitrogen.root_uptake.nitrogenUptake( ...
    Srt, NL, DeltZ, C_NH4, C_NO3, SoilState, RLD, daisyParams, nacryParams);

end
