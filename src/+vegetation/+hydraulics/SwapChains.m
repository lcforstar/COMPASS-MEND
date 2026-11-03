function [theta, logp1, swapflag] = SwapChains(temps, theta, logp1)
    n = length(temps);
    swapflag = zeros(size(temps));
    swap = randperm(n, 2);
    logA = (1 / temps(swap(2)) - 1 / temps(swap(1))) * (logp1(swap(1)) - logp1(swap(2)));
    if log(rand) < logA
        logp_tmp = logp1(swap(1));
        logp1(swap(1)) = logp1(swap(2));
        logp1(swap(2)) = logp_tmp;
        theta_tmp = theta(swap(1), :);
        theta(swap(1), :) = theta(swap(2), :);
        theta(swap(2), :) = theta_tmp;
        swapflag(swap) = 1;
    end
end
