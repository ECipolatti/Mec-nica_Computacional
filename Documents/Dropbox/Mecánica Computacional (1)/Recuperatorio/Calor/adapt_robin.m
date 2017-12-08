function [ROB] = adapt_robin(X, h, phi_inf)
    N = length(X);
    ROB = [X, h*ones(N,1), phi_inf*ones(N,1)];
end

