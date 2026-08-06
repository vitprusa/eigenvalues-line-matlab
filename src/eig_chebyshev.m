function lambdas = eig_chebyshev(N, a, b, q)
% Chebyshev differentiation matrix method (CDM).
%
% Returns the N eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet boundary
% conditions y(a) = y(b) = 0.  The operator is discretised on N + 2 Chebyshev
% points of the second kind; the two boundary rows and columns are removed to
% impose the boundary conditions, leaving an N by N matrix.
%
%   lambdas = eig_chebyshev(1000, 0, pi, @q_paine)
%
% Requires Chebfun (chebpts, diffmat).
arguments
    % matrix size / number of eigenvalues
    N (1, 1) {mustBeInteger, mustBePositive}
    % endpoints of the interval
    a (1, 1) {mustBeReal}
    b (1, 1) {mustBeReal, mustBeGreaterThan(b, a)}
    % potential q
    q function_handle
end

M = N + 2;                          % grid including both endpoints

x = chebpts(M, [a b]);
D2 = diffmat(M, 2, [a b]);

% Drop the boundary points: Dirichlet conditions
x_int = x(2:end-1);
D2_int = D2(2:end-1, 2:end-1);

A = -D2_int + diag(q(x_int));

lambdas = sort(real(eig(A)));

end
