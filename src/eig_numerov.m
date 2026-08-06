function lambdas = eig_numerov(N, a, b, q)
% Numerov method.
%
% Returns the N eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet boundary
% conditions y(a) = y(b) = 0, discretised by the Numerov scheme on the uniform
% interior grid of N points.  The scheme
%
%   (y_{j-1} - 2y_j + y_{j+1})/h^2 = (f_{j-1} + 10 f_j + f_{j+1})/12,
%   f = (q - lambda) y,
%
% gives the generalised eigenvalue problem A y = lambda B y solved below.
%
%   lambdas = eig_numerov(1000, 0, pi, @q_paine)
%
% See eig_numerov_corrected for the variant with the asymptotic correction.
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

qv = q(x);

% spdiags takes the sub- and superdiagonals from the same column vector with
% the appropriate shift, so both off-diagonal entries pick up q at their own
% grid point, as the Numerov stencil requires.
main = 24 + 10 * h^2 * qv;
off = -(12 - h^2 * qv);

A = spdiags([off, main, off], -1:1, N, N);
B = h^2 * spdiags(ones(N, 1) * [1 10 1], -1:1, N, N);

lambdas = sort(real(eig(full(A), full(B))));

end
