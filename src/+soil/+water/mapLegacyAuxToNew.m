function aux = mapLegacyAuxToNew(SoilVariables)
%MAPLEGACYAUXTONEW  Pull the 13 Stage-B refresh-aux fields from a legacy SoilVariables
%   struct into the built.workspace.derived naming/orientation, for the Stage-B oracle
%   (legacy UpdateSoilWaterContent output vs new updateWaterContentAuxiliary output).
%
%   Inputs:
%     SoilVariables  [in]  legacy soil-variables struct (UpdateSoilWaterContent output)
%
%   Returns:
%     aux  struct of the 13 compared aux fields named as built.workspace.derived:
%          Theta_LL/Theta_UU/Theta_II/KL_h/KfL_h/KfL_T/DTheta_LLh/DTheta_UUh/Ratio_ice
%          (endpoint [NL x nD], same name), Gamma_hh (node [NN x 1], same name),
%          COR/CORh reshaped from legacy row [1 x NN] to new column [NN x 1], and Se
%          SEMANTICALLY REMAPPED (see below), not pulled from SoilVariables.Se.
%
%   Method: 10 fields are same-name pulls; COR/CORh change orientation; Se is a semantic
%   remap. The new frame defines Se = Theta_LL/POR (air/pore-space saturation, for the gas
%   and vapor auxiliaries; see computeEndpointHydraulics.m calculateSeArray). The legacy
%   main path instead writes SoilVariables.Se = the VG effective saturation Se_eff returned
%   by calculateTheta_LL (Zone 4: (Theta_LL-theta_r)/(theta_s-theta_r)) -- a legacy field-
%   definition inconsistency: UpdateSoilWaterContent.m's own doc says Se is "for the gas
%   equation" and its resync branch (:288) writes Theta_LL/POR, yet the main path stores the
%   VG Se_eff. So SoilVariables.Se is NOT the physical oracle for the new Se_air; we recompute
%   the air-basis Se from the legacy Theta_LL/POR (both bit-identical to the new frame in the
%   unsaturated regime). The K path is unaffected: KL_h/KfL_h use their OWN internal VG
%   SeK=(Theta_LL-theta_r)/(theta_s-theta_r) on both sides and are compared directly.
%   Written as an explicit field list so a future rename on either side is caught.
%
%   Reference: docs/superpowers/specs/2026-07-07-stage4-water-only-migration-oracle-design.md
%   sec 4 (Stage B); field set confirmed against src/+soil/+water/updateWaterContentAuxiliary.m
%   and src/UpdateSoilWaterContent.m.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-07 Jizixing Zhu -- Se is now a semantic remap to the new air-basis
%           Se_air (Theta_LL/POR), not a same-name pull of the legacy VG Se_eff.
    same = {'Theta_LL','Theta_UU','Theta_II','KL_h','KfL_h','KfL_T', ...
            'DTheta_LLh','DTheta_UUh','Ratio_ice','Gamma_hh'};
    aux = struct();
    for i = 1:numel(same)
        aux.(same{i}) = SoilVariables.(same{i});
    end
    aux.COR  = SoilVariables.COR(:);    % legacy [1 x NN] row -> new [NN x 1] column
    aux.CORh = SoilVariables.CORh(:);
    % Se: remap legacy Theta_LL/POR to the new-frame air-basis Se_air (POR [NL x 1]
    % broadcasts across the [NL x nD] endpoints), clamped [0,1] like calculateSeArray.
    aux.Se = min(1, max(0, SoilVariables.Theta_LL ./ SoilVariables.POR(:)));
end
