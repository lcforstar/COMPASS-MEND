function Fout = meanleaf(canopy, F, choice, Ps)

    nl    = canopy.nlayers;
    nli   = canopy.nlincl;
    nlazi = canopy.nlazi;
    lidf  = canopy.lidf;

    % author: Dr. ir. Christiaan van der Tol (tol@itc.nl)
    % date:     7   December 2007
    % update:   11  February 2008 made modular (Joris Timmermans)

    % update:   25 Feb 2013 Wout Verhoef : Propose name change, remove globals
    %                                      and use canopy-structure for input
    %
    % function [F2,F3] = F1tot(F,choice,Ps)
    % calculates the layer average and the canopy average of leaf properties
    % per layer, per leaf angle and per leaf azimuth (36)
    %
    % Input:
    %   F       input matrix (3D)   [nli, nlazi,nl]
    %   choice  integration method  'angles'            : integration over leaf angles
    %                               'angles_and_layers' : integration over leaf layers and leaf angles
    %   Ps      fraction sunlit per layer [nl]
    %
    % Output:
    %   Fout    in case of choice = 'angles': [nl]
    %           in case of choice = 'angles_and_layers': [1]
    Fout = zeros(nli, nlazi, nl);
    switch choice
        %% integration over leaf angles
        case 'angles'

            for j = 1:nli
                Fout(j, :, :)    = F(j, :, :) * lidf(j);         % [nli, nlazi,nl]
            end
            Fout                = sum(sum(Fout)) / nlazi;     % [1,1,nl]
            Fout                = permute(Fout, [3 1 2]);    % [nl]

            %% integration over layers only
        case 'layers'
            % not implemented
            %% integration over both leaf angles and layers
        case 'angles_and_layers'
            for j = 1:nli
                Fout(j, :, :)      = F(j, :, :) * lidf(j);
            end

            for j = 1:nl
                Fout(:, :, j)    = Fout(:, :, j) * Ps(j);
            end
            Fout                = sum(sum(sum(Fout))) / nlazi / nl;

        case 'angles_and_layers_weighted'
            % Like angles_and_layers but replaces the uniform 1/nl layer
            % weight with canopy.layer_weight. Requires canopy.layer_weight
            % [nl×1] and Ps [nl×1] as the fourth argument.
            % Returns the Ps-weighted, layer_weight-scaled canopy mean:
            %   Fout = sum( layer_weight .* Ps .* mean_angles(F) )
            % Note: this is the sunlit branch only. To get a full two-branch
            % canopy integral, compute the shaded branch explicitly:
            %   Ftot = sum(LAI_layer .* ((1-Ps) .* Fsh + Ps .* meanleaf(F,'angles')))
            for j = 1:nli
                Fout(j, :, :) = F(j, :, :) * lidf(j);
            end
            F1d  = permute(sum(sum(Fout), 1), [3 1 2]) / nlazi;  % [nl×1]
            Fout = sum(canopy.layer_weight .* Ps(1:nl) .* F1d);
    end
