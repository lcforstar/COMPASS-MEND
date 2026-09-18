function localMeteo = buildDemandMeteoView(built, kt)
%BUILDDEMANDMETEOVIEW  2-slot air-temperature view so runPlantNDemand reads real-kt forcing under kt2=2.
%
%   localMeteo = vegetation.integrated.buildDemandMeteoView(built, kt) returns a minimal
%   MeteoData-shaped struct whose only field, Ta, is a 2-element column [Ta; Ta] with
%   Ta = built.forcing.hourly.Ta(kt), the REAL step-kt air temperature.
%
%   Why: the PR-D plant-C/N bridge calls the reused legacy functions with a LOCAL kt2=2
%   (so their pool reads at idx_prev=max(1,kt2-1)=1 hit the accepted scalar). But
%   runPlantNDemand ALSO indexes real forcing with the SAME argument -- `MeteoData.Ta(KT)`
%   (runPlantNDemand.m:135) is its one in-function forcing-index read. Calling it with the
%   real MeteoData under kt2=2 would read Ta of forcing step 2, not step kt. This 2-slot
%   view neutralises that: with kt2=2, `localMeteo.Ta(2)` returns the real step-kt Ta.
%   runPlantNDemand reads no other MeteoData field and calls no MeteoData method, so a plain
%   struct with a single Ta field suffices (kept minimal; the reused function stays verbatim).
%
%   Frame note: built.forcing.hourly.Ta(kt) == built.workspace.scope_step.meteo.Ta
%   (prepareScopeStep.m:258 direct slice, no transform) == legacy MeteoData.Ta(kt). The raw
%   forcing is used here (not the scope_step slice) so the source is the same field legacy
%   runPlantNDemand read, independent of whether prepareCanopyStep has run.
%
%   Inputs:
%     built  struct  init container; reads built.forcing.hourly.Ta (tS x 1) [degC]
%     kt     [-]     real outer-step index (1-based)
%
%   Returns:
%     localMeteo  struct  .Ta = [Ta_kt; Ta_kt]  [2 x 1, degC]  (slot 2 = the kt2=2 read)
%
%   Variable class (per compass-coding-standard §1):
%     built.forcing.hourly.Ta -> forcing (in);  localMeteo -> a per-step forcing view (out).
%
%   Reference: COMPASS.m:1257-1259 (glue call to runPlantNDemand); runPlantNDemand.m:135
%     (MeteoData.Ta(KT)); prepareScopeStep.m:258 (Ta slice); docs/pr_d_plant_cn_spec.md §4.1.
%   LEGACY ADAPTER LAYER (Julia-port): a transitional boundary adapter giving the reused verbatim
%   runPlantNDemand its old MeteoData.Ta(KT) shape under the 2-slot remap; deletes when that physics
%   ports to (dY,Y,p,aux,t) and reads forcing at the real kt directly.
% Editor: Jizixing Zhu; Date: 2026-07-04.
% Change: 2026-07-04 Jizixing Zhu -- doc: mark as the transitional legacy adapter layer (Julia-port).
    Ta = built.forcing.hourly.Ta(kt);      % [degC] real step-kt air temperature (== legacy MeteoData.Ta(kt))
    localMeteo = struct('Ta', [Ta; Ta]);   % [2 x 1] 2-slot view; runPlantNDemand(kt2=2) reads Ta(2)
end
