function [Zr1, Zr2, RAI1, RAI2] = calculateRootAreaIndex(root_type, root_depth)
% calculateRootDistribution 计算土壤根系分布参数
% 该函数根据输入的根系类型和深度，计算表层和下层根系深度及其面积指数
% from : A global analysis of root distributions for terrestrial biomes  
% A global budget for fine root biomass, surface area and nutrient contents
% 
%
% 输入参数:
%   root_type:    根系类型索引（整数，应在RAIlist和betalist的有效索引范围内）
%   root_depth:   根系总深度 (cm)

%
% 输出参数:
%   Zr1:  表层根系深度 (cm)
%   Zr2:  深层层根系深度 (cm)
%   RAI1: 表层根系面积指数
%   RAI2: 深层层根系面积指数

    % 定义根系参数列表（这些值应根据实际研究数据进行调整）
    RAIlist = [4.6, 10, 5.5, 11.6, 11.0, 9.8, 79.1, 6.3, 7.4, 42.5, 5.2];  % leaf fine root area index(m2.m-2)
    betalist = [0.943, 0.961, 0.970, 0.950, 0.980, 0.967, 0.943, 0.982, 0.972, 0.972, 0.909]; % beta
    
    % 计算根系深度分布
    Zr1 = min(30, root_depth);           % 表层根系深度（不超过30cm）
    Zr2 = max(0, root_depth - 30);       % 下层根系深度（总深度减30cm，且不小于0）
    
    % (RAI)
    if Zr2 > 0
       cumulative_rootfraction = betalist(root_type)^30;
    else
       cumulative_rootfraction = 0;
    end
    
    RAI1 = RAIlist(root_type) * (1 - cumulative_rootfraction);  % 表层RAI
    RAI2 = RAIlist(root_type) - RAI1;                           % 下层RAI
    
 
end