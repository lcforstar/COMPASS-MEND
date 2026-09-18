function reportDryAirDebug(tag, kt, kit, deltT, SoilVariables, TransportCoefficient, AirVars, cfg, errWater, errEnergy, PBeforeAir)
%REPORTDRYAIRDEBUG Print compact diagnostics for dry-air Picard failures.
% Editor: Jizixing Zhu; Date: 2026-07-13.
% Change: 2026-07-13 Jizixing Zhu -- report gas pressure mobility as Kg.

if ~(cfg.lSoilairefc == 1)
    return;
end

fprintf('[DRY-AIR-DEBUG] %s KT=%d KIT=%d dt=%.3g Water=% .3e Energy=% .3e Freezing=%d SFCC=%d\n', ...
    tag, kt, kit, deltT, errWater, errEnergy, cfg.lFreezing, cfg.lSFCC);

printRange('hh', SoilVariables.hh);
printRange('TT', SoilVariables.TT);
printRange('P_g', SoilVariables.P_g);
printRange('P_gg', SoilVariables.P_gg);
if nargin >= 11 && ~isempty(PBeforeAir)
    dP = SoilVariables.P_gg(:) - PBeforeAir(:);
    printRange('dP_air', dP);
    [~, idx] = max(abs(dP));
    fprintf('    max|dP_air| node=%d value=% .6e Pa\n', idx, dP(idx));
end

printRange('Theta_LL', SoilVariables.Theta_LL);
printRange('Theta_g', SoilVariables.Theta_g);
printRange('DTheta_LLh', SoilVariables.DTheta_LLh);
printRange('KfL_h', SoilVariables.KfL_h);
printRange('KfL_T', SoilVariables.KfL_T);
if hasMember(SoilVariables, 'Ratio_ice')
    printRange('Ratio_ice', SoilVariables.Ratio_ice);
end

if ~isempty(AirVars)
    printRange('Caa', AirVars.Caa);
    printRange('Kaa', AirVars.Kaa);
    printRange('Vaa', AirVars.Vaa);
    printRange('Cah', AirVars.Cah);
    printRange('CaT', AirVars.CaT);
    printRange('Kah', AirVars.Kah);
    printRange('KaT', AirVars.KaT);
    if hasMember(AirVars, 'QL'),    printRange('QL', AirVars.QL);       end
    if hasMember(AirVars, 'QL_A'),  printRange('QL_A', AirVars.QL_A);   end
    if hasMember(AirVars, 'DPgDZ'), printRange('DPgDZ', AirVars.DPgDZ); end
end

printRange('Kg', TransportCoefficient.Kg);
printRange('D_Vg', TransportCoefficient.D_Vg);
printRange('V_A', TransportCoefficient.V_A);
end

function tf = hasMember(obj, name)
if isstruct(obj)
    tf = isfield(obj, name);
elseif isobject(obj)
    tf = isprop(obj, name);
else
    tf = false;
end
end

function printRange(name, value)
v = value(:);
finiteMask = isfinite(v);
nBad = nnz(~finiteMask);
if any(finiteMask)
    vf = v(finiteMask);
    fprintf('    %-12s min=% .6e max=% .6e mean=% .6e bad=%d\n', ...
        name, min(vf), max(vf), mean(vf), nBad);
else
    fprintf('    %-12s all non-finite bad=%d\n', name, nBad);
end
end
