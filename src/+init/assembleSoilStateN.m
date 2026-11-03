function SoilStateN = assembleSoilStateN(built)
%ASSEMBLESOILSTATEN  Reconstruct the legacy SoilStateN slice from the new split containers.
%
%   SoilStateN = init.assembleSoilStateN(built) gathers the soil-state fields that
%   vegetation.root.applyRootUptake / soil.nitrogen.root_uptake.nitrogenUptake consume
%   (legacy COMPASS.m:1204-1209), from the new-chain containers. This is a discrete
%   intermediate for the PR-D root-N-uptake wiring (NOT the NitrogenVariables
%   assemble/scatter bridge); it reads only, mutates nothing.
%
%   The one non-trivial field is TT (soil node temperature): legacy sets
%   SoilStateN.TT = SoilVariables.TT (COMPASS.m:1209), an NN x 1 bottom-up node array
%   (defineSoilVariables.m:16). Its new-chain origin is the prognostic Y.T
%   (buildInitialState.m -> bottom-up, node 1 = deepest). Y.T and SoilVariables.TT
%   share the bottom-up node frame, so TT = built.Y.T(:) with NO flip. This is pinned
%   by tests/test_assembleSoilStateN.m (KT=1 parity vs legacy SoilVariables.TT).
%
%   Inputs:
%     built  struct  init container: .Y.T, .workspace.derived.Theta_LL,
%                    .parameters.soil.theta_s, .parameters.soil_n.{RHO_bulk,RHOKG,DIF1111}
%
%   Returns SoilStateN struct:
%     .Theta_LL  [NL x 2, -]        liquid water content (solver working frame)
%     .Theta_s   [NL x 1, -]        saturated water content (parameters.soil, NOT soil_n)
%     .RHO_bulk  [NL x 1, g/cm3]    dry bulk density
%     .RHOKG     [NL x 1, -]        NH4 adsorption buffer (Kd * bulk density)
%     .DIF1111   [1 x NL x nS]      solute molecular diffusion (NH4=:,:,1 / NO3=:,:,2)
%     .TT        [NN x 1, degC]     soil node temperature (bottom-up = legacy SoilVariables.TT)
%
%   Reference: COMPASS.m:1204-1209 (SoilStateN assembly); defineSoilVariables.m:16 (TT frame).
%   LEGACY ADAPTER LAYER (Julia-port): a transitional boundary adapter feeding the reused verbatim
%   root-uptake physics its old SoilStateN shape; deletes when that physics ports to (dY,Y,p,aux,t).
% Editor: Jizixing Zhu; Date: 2026-07-04.
% Change: 2026-07-04 Jizixing Zhu -- doc: mark as the transitional legacy adapter layer (Julia-port).
    SoilStateN = struct();
    SoilStateN.Theta_LL = built.workspace.derived.Theta_LL;   % NL x 2 [-] working frame
    SoilStateN.Theta_s  = built.parameters.soil.theta_s;      % NL x 1 [-] (parameters.soil, NOT soil_n -- trap)
    SoilStateN.RHO_bulk = built.parameters.soil_n.RHO_bulk;   % NL x 1 [g/cm3]
    SoilStateN.RHOKG    = built.parameters.soil_n.RHOKG;      % NL x 1 [-]
    SoilStateN.DIF1111  = built.parameters.soil_n.DIF1111;    % 1 x NL x nS (singleton lead dim)
    SoilStateN.TT       = built.Y.T(:);                       % NN x 1 [degC] bottom-up node temp (= legacy SoilVariables.TT)
end
