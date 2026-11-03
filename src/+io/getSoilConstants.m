function SoilConstants = getSoilConstants()
%GETSOILCONSTANTS  Legacy IO wrapper for soil-water model constants.
%
%   SoilConstants = io.getSoilConstants() is kept for legacy callers that
%   still import constants through the io package. New soil-water code should
%   read built.parameters.soil.constants or call soil.water.constants().
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    SoilConstants = soil.water.constants();

    % Legacy heterogeneity initialization still reads these two fields through
    % io.getSoilConstants. They are intentionally not part of the canonical
    % built.parameters.soil.constants container.
    SoilConstants.CHST = 0;
    SoilConstants.Phi_S = [-17.9 -17 -17 -19 -10 -10];
end
