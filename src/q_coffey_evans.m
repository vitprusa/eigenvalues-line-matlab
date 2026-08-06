function y = q_coffey_evans(x, beta)
% Coffey--Evans potential q(x) = -2*beta*cos(2x) + beta^2*sin(2x)^2, used on
% the interval [-pi/2, pi/2].  Default beta = 30.
%
% Pass as a function handle to any of the eig_* solvers, e.g.
%   lambdas = eig_dst(1001, -pi/2, pi/2, @q_coffey_evans);
% For a value of beta other than the default, wrap it:
%   lambdas = eig_dst(1001, -pi/2, pi/2, @(x) q_coffey_evans(x, 50));
arguments
    x double
    beta (1, 1) {mustBeReal} = 30
end

y = -2*beta * cos(2*x) + beta^2 * sin(2*x).^2;

end
