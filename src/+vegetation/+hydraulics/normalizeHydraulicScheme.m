function scheme = normalizeHydraulicScheme(raw)
%NORMALIZEHYDRAULICSCHEME  Map config aliases to a canonical hydraulic scheme.
%
%   scheme = vegetation.hydraulics.normalizeHydraulicScheme(raw) returns one of:
%     'steady_spac'      - v1 steady SPAC (calSPAC); no organ storage W
%     'storage_from_psi' - diagnose W from psi via retention (VOD path)
%     'storage_dynamic'  - prognostic dW/dt storage (shares V2 wrapper today)
%
%   Legacy aliases (still accepted):
%     'equilibrium' -> 'steady_spac'
%     'dynamic'     -> 'storage_from_psi'
%
%   Empty or missing input defaults to 'steady_spac'.
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.
    if nargin < 1 || isempty(raw)
        scheme = 'steady_spac';
        return;
    end
    key = lower(strtrim(char(raw)));
    switch key
        case {'steady_spac', 'equilibrium'}
            scheme = 'steady_spac';
        case {'storage_from_psi', 'dynamic'}
            scheme = 'storage_from_psi';
        case 'storage_dynamic'
            scheme = 'storage_dynamic';
        otherwise
            error('PlantHydraulic:UnknownHydraulicScheme', ...
                ['Unknown hydraulic_scheme "%s". Use steady_spac, ', ...
                 'storage_from_psi, or storage_dynamic ', ...
                 '(legacy: equilibrium, dynamic).'], char(raw));
    end
end
