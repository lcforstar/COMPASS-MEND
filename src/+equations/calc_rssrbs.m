function [rss, rbs] = calc_rssrbs(SoilProperties, LAI, rbs)

    % an alias

    aa = 3.8;
    rss = exp((aa + 4.1) - aa * (SoilProperties.SMC - SoilProperties.Theta_r_type(1)) / (SoilProperties.Field_mc_type(1) - SoilProperties.Theta_r_type(1)));
    rbs            = rbs * LAI / 4.3;
