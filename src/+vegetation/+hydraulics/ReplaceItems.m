function theta = ReplaceItems(theta, theta_star, accept)
    for i = 1:length(accept)
        if accept(i)
            theta(i) = theta_star(i);
        end
    end
end
