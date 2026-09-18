function [SoilVariables, VanGenuchten] = applySoilHeteroEffect(SoilProperties, SoilVariables, VanGenuchten)

 
    SoilVariables.Phi_s = []; % see issue 139, item 5
    SoilVariables.Lamda = []; % see issue 139, item 5

    % get model settings
    ModelSettings = io.getModelSettings();

    % Get soil constants for StartInit
    SoilConstants = io.getSoilConstants();

    for i = 1:6
        if ModelSettings.SWCC == 0
            if SoilConstants.CHST == 0
                SoilVariables.Phi_s(i) = SoilConstants.Phi_S(i);
                SoilVariables.Lamda(i) = SoilProperties.Coef_Lamda(i);
            end
        end
    end

   
   