function [W_new, psi_v2_MPa, diag] = Plant_hydraulics_v2_wrapper(varargin)
%PLANT_HYDRAULICS_V2_WRAPPER  Deprecated alias for updateOrganStorageFromPsi.
%
%   Prefer vegetation.hydraulics.updateOrganStorageFromPsi.
%   Editor: Jizixing Zhu; Date: 2026-07-17.
    [W_new, psi_v2_MPa, diag] = vegetation.hydraulics.updateOrganStorageFromPsi( ...
        varargin{:});
end
