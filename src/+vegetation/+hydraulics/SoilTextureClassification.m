function [tx] =  SoilTextureClassification(FOC, FOS, FOSL, NS)
%SoilTextureClassification based on particle composition.
%   This function implements the Clapp and Hornberger (1978) soil texture
%   classification scheme. It determines the soil texture class based on
%   the percentage of clay, sand, and silt in the soil.
%
%   Input Parameters:
%   - FOC: Fraction of Clay in soil (%, between 0-100)
%   - FOS: Fraction of Sand in soil (%, between 0-100)
%   - FOSL: Fraction of Silt in soil (%, between 0-100)
%   - NS: Number of soil texture
%
%   Output Parameter:
%   - tx: Soil texture class index (integer from 1 to 11, or NaN if invalid)
%         See comments within the code for the meaning of each index value.
%
%   Example:
%   [textureClass] = vegetation.hydraulics.SoilTextureClassification(20, 30, 50); % Classifies soil with 20% clay, 30% sand, 50% silt
%
%   Reference:
%   Clapp, R. B., & Hornberger, G. M. (1978). Empirical equations for some soil hydraulic properties. 
%   Water Resources Research, 14(4), 601-604.
%
%   Author: [Zhu jizixing]
%   Date: 2025-09-09
    sand = FOS.*100;
    silt = FOSL.*100;
    clay = FOC.*100;
    tx = zeros(1,NS);
    for i = 1 : NS
        if silt(i) + 1.5 * clay(i) < 15
            tx(i) = 1; % sand(i)
        elseif silt(i) + 1.5 * clay(i) >= 15 && silt(i) + 2 * clay(i) < 30
            tx(i) = 2; % loamy sand(i)
        elseif (7 <= clay(i) && clay(i) < 20 && sand(i) > 52 && silt(i) + 2 * clay(i) >= 30) || (clay(i) < 7 && silt(i) < 50 && silt(i) + 2 * clay(i) >= 30)
            tx(i) = 3; % sandy loam
        elseif (silt(i) >= 50 && 12 <= clay(i) && clay(i) < 27) || (50 <= silt(i) && silt(i) < 80 && clay(i) < 12)
            tx(i) = 4; % silt(i) loam or silt(i)
        elseif (7 <= clay(i) && clay(i) < 27 && 28 <= silt(i) && silt(i) < 50 && sand(i) <= 52)
            tx(i) = 5; % loam
        elseif (20 <= clay(i) && clay(i) < 35 && silt(i) < 28 && sand(i) > 45)
            tx(i) = 6; % sandy clay(i) loam
        elseif (27 <= clay(i) && clay(i) < 40 && sand(i) <= 20)
            tx(i) = 7; % silt(i) clay(i) loam
        elseif (27 <= clay(i) && clay(i) < 40 && 20 < sand(i) && sand(i) <= 45)
            tx(i) = 8; % clay(i) loam
        elseif (clay(i) >= 35 && sand(i) > 45)
            tx(i) = 9; % sandy clay(i)
        elseif (clay(i) >= 40 && silt(i) >= 40)
            tx(i) = 10; % silt(i) clay(i)
        elseif (clay(i) >= 40 && sand(i) <= 45 && silt(i) < 40)
            tx(i) = 11; % clay(i)
        else
            tx(i) = NaN;
        end
    end
end

   
    
    
    

