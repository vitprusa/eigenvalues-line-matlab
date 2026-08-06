function lambdas = eig_fd(N, a, b, q)
% Second order finite difference method.
%
% Returns the N eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet boundary
% conditions y(a) = y(b) = 0, discretised by the standard three point stencil
% on the uniform interior grid of N points.
%
%   lambdas = eig_fd(10000, 0, pi, @q_paine)
arguments
    % matrix size / number of eigenvalues
    N (1, 1) {mustBeInteger, mustBePositive}
    % endpoints of the interval
    a (1, 1) {mustBeReal}
    b (1, 1) {mustBeReal, mustBeGreaterThan(b, a)}
    % potential q
    q function_handle
end

h = (b - a) / (N + 1);
x = a + h * (1:N)';

D2 = spdiags(ones(N, 1) * [1 -2 1], -1:1, N, N) / h^2;

A = -full(D2) + diag(q(x));

lambdas = sort(real(eig(A)));

end
