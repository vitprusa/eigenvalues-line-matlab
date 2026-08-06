function lambdas = eig_chebyshev_mapped(N, a, b, q, alpha)
% Mapped barycentric Chebyshev differentiation matrix method (MBCDM).
%
% Returns the N eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet boundary
% conditions y(a) = y(b) = 0.  The Chebyshev points are pushed towards a
% uniform distribution by the Kosloff--Tal-Ezer map
%
%   g(s) = asin(alpha*s) / asin(alpha),   x = a + (b-a)*(g(s)+1)/2,
%
% which relaxes the O(N^-2) clustering of the grid at the endpoints.  The
% parameter alpha in (0, 1] controls the strength of the map: alpha -> 0
% recovers the plain Chebyshev grid, alpha = 1 gives a uniform grid.
%
%   lambdas = eig_chebyshev_mapped(1000, 0, pi, @q_paine)
%   lambdas = eig_chebyshev_mapped(1001, -pi/2, pi/2, @q_coffey_evans, 1)
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
    % strength of the Kosloff--Tal-Ezer map
    alpha (1, 1) {mustBeReal, mustBePositive, mustBeLessThanOrEqual(alpha, 1)} = 0.999
end

M = N + 2;                          % grid including both endpoints
L = b - a;

s = chebpts(M);                     % Chebyshev points on [-1, 1]

% Map and its first two derivatives
asin_alpha = asin(alpha);
g   = asin(alpha * s) / asin_alpha;
gp  = alpha ./ sqrt(1 - (alpha * s).^2) / asin_alpha;
gpp = (alpha^3 * s) ./ (1 - alpha^2 * s.^2).^(3/2) / asin_alpha;

x       = a + (L/2) * (g + 1);
dx_ds   = (L/2) * gp;
d2x_ds2 = (L/2) * gpp;

% Chain rule: d2/dx2 = (dx/ds)^-2 * d2/ds2 - (d2x/ds2)/(dx/ds)^3 * d/ds.
% For alpha = 1 the map is singular at s = +-1, which pollutes the first and
% last row only; those rows are removed with the boundary points below.
D_s  = diffmat(M, 1);
D_s2 = diffmat(M, 2);
D2 = diag(1 ./ dx_ds.^2) * D_s2 - diag(d2x_ds2 ./ dx_ds.^3) * D_s;

% Drop the boundary points: Dirichlet conditions
x_int = x(2:end-1);
D2_int = D2(2:end-1, 2:end-1);

A = -D2_int + diag(q(x_int));

lambdas = sort(real(eig(A)));

end
