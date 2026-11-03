function tf = getSwitch(cfg, name, dflt)
% Read a logical on/off switch from cfg.switches, returning dflt if absent.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    tf = dflt;
    if isfield(cfg, 'switches') && isfield(cfg.switches, name)
        tf = logical(cfg.switches.(name));
    end
end
