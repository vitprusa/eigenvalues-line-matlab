function y = q_paine(x)
% Paine's potential q(x) = exp(x), used on the interval [0, pi].
%
% Pass as a function handle to any of the eig_* solvers, e.g.
%   lambdas = eig_dst(1000, 0, pi, @q_paine);
arguments
    x double
end

y = exp(x);

end
