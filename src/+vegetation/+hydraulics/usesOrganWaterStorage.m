function tf = usesOrganWaterStorage(scheme)
%USESORGANWATERSTORAGE  True when organ water storage W is active.
%
%   tf = vegetation.hydraulics.usesOrganWaterStorage(scheme) is true for
%   storage_from_psi and storage_dynamic (and their legacy aliases after
%   normalizeHydraulicScheme).
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.
    s = vegetation.hydraulics.normalizeHydraulicScheme(scheme);
    tf = any(strcmp(s, {'storage_from_psi', 'storage_dynamic'}));
end
