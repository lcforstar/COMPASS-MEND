function s = getSwitchString(cfg, name, dflt)
% Read a string-valued switch from cfg.switches, returning dflt if absent.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    s = dflt;
    if isfield(cfg, 'switches') && isfield(cfg.switches, name)
        s = char(cfg.switches.(name));
    end
end
