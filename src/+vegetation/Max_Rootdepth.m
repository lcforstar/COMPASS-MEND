%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction - Root - Depth         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
function [bbx] = Max_Rootdepth(cfg, root_max_depth)

    % get model settings
    NL        = cfg.Grid.NL;
    NN        = cfg.Grid.NN;
    DeltZ     = cfg.Grid.DeltZ;
   Tot_Depth = cfg.Grid.Tot_Depth;

    % BR = 10:1:650; %% [gC /m^2 PFT]
    % rroot =  0.5*1e-3 ; % 3.3*1e-4 ;%% [0.5-6 *10^-3] [m] root radius
    %%% OUTPUTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Root lenth distribution %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_depth = root_max_depth;
    Tot_Depth = Tot_Depth;
    DeltZ = DeltZ;
    Elmn_Lnth = 0;
    for ML = 1:NL
        Elmn_Lnth = Elmn_Lnth + DeltZ(ML);
        if Elmn_Lnth < Tot_Depth - R_depth
            bbx(ML) = 0;
        elseif Elmn_Lnth >= Tot_Depth - R_depth && Elmn_Lnth <= Tot_Depth - 5
            bbx(ML) = 1;
        else
            bbx(ML) = 0;
        end
    end
end
