function reportAcceptedSoilDiagnostics(kt, currentDate, SoilVariables, TransportCoefficient, AirVars, errWater, errEnergy, substepCount)
%REPORTACCEPTEDSOILDIAGNOSTICS Print compact diagnostics for accepted states.
% Editor: Jizixing Zhu; Date: 2026-07-13.
% Change: 2026-07-13 Jizixing Zhu -- report gas pressure mobility as Kg.

if nargin < 8 || isempty(substepCount)
    substepCount = NaN;
end

if isempty(currentDate) || (isdatetime(currentDate) && isnat(currentDate))
    dateText = 'NaT';
else
    dateText = char(currentDate);
end

fprintf('[SOIL-DIAG] KT=%d Date=%s AcceptedSubsteps=%d Water=% .3e Energy=% .3e\n', ...
    kt, dateText, substepCount, errWater, errEnergy);

printRange('hh', SoilVariables.hh);
printRange('TT', SoilVariables.TT);
printRange('Theta_LL', SoilVariables.Theta_LL);
printRange('Theta_II', SoilVariables.Theta_II);
printRange('Theta_g', SoilVariables.Theta_g);
printRange('DTheta_LLh', SoilVariables.DTheta_LLh);
printRange('DTheta_UUh', SoilVariables.DTheta_UUh);
printRange('KL_h', SoilVariables.KL_h);
printRange('KfL_h', SoilVariables.KfL_h);
printRange('KfL_T', SoilVariables.KfL_T);
printRange('Ratio_ice', SoilVariables.Ratio_ice);

if ~isempty(AirVars)
    if hasMember(AirVars, 'QL'),    printRange('QL', AirVars.QL);       end
    if hasMember(AirVars, 'QL_A'),  printRange('QL_A', AirVars.QL_A);   end
    if hasMember(AirVars, 'DPgDZ'), printRange('DPgDZ', AirVars.DPgDZ); end
    if hasMember(AirVars, 'Cah'),   printRange('Cah', AirVars.Cah);     end
    if hasMember(AirVars, 'CaT'),   printRange('CaT', AirVars.CaT);     end
end

if ~isempty(TransportCoefficient)
    if hasMember(TransportCoefficient, 'Kg'), printRange('Kg', TransportCoefficient.Kg); end
    if hasMember(TransportCoefficient, 'D_Vg'),   printRange('D_Vg', TransportCoefficient.D_Vg);     end
    if hasMember(TransportCoefficient, 'V_A'),    printRange('V_A', TransportCoefficient.V_A);       end
end
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
if isempty(value)
    fprintf('    %-12s empty\n', name);
    return;
end

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
