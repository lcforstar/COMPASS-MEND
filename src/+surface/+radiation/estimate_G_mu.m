function G_mu = estimate_G_mu(mu, leaf_angle_distribution)
%ESTIMATE_G_MU 估算叶片投影函数G(μ)
%   输入:
%     mu: 太阳天顶角余弦
%     leaf_angle_distribution: 叶角分布类型
%       'spherical' - 球型分布（默认）
%       'horizontal' - 水平分布
%       'vertical' - 垂直分布
%
%   输出:
%     G_mu: 叶片在太阳方向的投影函数

    if nargin < 2
        leaf_angle_distribution = 'spherical';
    end
    
    switch leaf_angle_distribution
        case 'spherical'
            % 球型叶角分布：G(μ) = 0.5
            G_mu = 0.5;
        case 'horizontal'
            % 水平叶角分布：G(μ) = μ
            G_mu = mu;
        case 'vertical'
            % 垂直叶角分布：G(μ) = 2/π * sqrt(1 - μ^2)
            G_mu = (2/pi) * sqrt(1 - mu^2);
        otherwise
            error('不支持的叶角分布类型');
    end
end