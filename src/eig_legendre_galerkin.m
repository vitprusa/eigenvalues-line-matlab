function lambdas = eig_legendre_galerkin(N, a, b, q)
% Legendre--Galerkin method with Chebyshev-computed coefficients (LGCC).
%
% Returns the N eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet boundary
% conditions y(a) = y(b) = 0.  The problem is pulled back to [-1, 1] by
% x = a + (b-a)*(t+1)/2, under which
%
%   -y_tt + ((b-a)/2)^2 q y = mu y,   lambda = mu / ((b-a)/2)^2,
%
% and discretised in the compact Legendre basis phi_k = P_k - P_{k+2}, which
% satisfies the boundary conditions.  The stiffness matrix S is diagonal and
% the mass matrix T pentadiagonal; the potential matrix is assembled from the
% Legendre coefficients of q*phi_j, computed via Chebfun.
%
%   lambdas = eig_legendre_galerkin(500, 0, pi, @q_paine)
%
% Requires Chebfun (chebfun, legpoly, legcoeffs).  The assembly is O(N^2)
% chebfun operations and is by far the slowest of the methods here.
arguments
    % number of basis functions / number of eigenvalues
    N (1, 1) {mustBeInteger, mustBePositive}
    % endpoints of the interval
    a (1, 1) {mustBeReal}
    b (1, 1) {mustBeReal, mustBeGreaterThan(b, a)}
    % potential q
    q function_handle
end

L = b - a;
ks = (0:N-1)';

% Stiffness: integral of phi_i' phi_j' over [-1, 1]
S = diag(4*ks + 6);

% Mass: integral of phi_i phi_j over [-1, 1]
T = zeros(N);
for i = 1:N
    k = ks(i);
    T(i, i) = 2/(2*k+1) + 2/(2*k+5);
    if i + 2 <= N
        T(i, i+2) = -2/(2*k+5);
        T(i+2, i) = -2/(2*k+5);
    end
end

% Potential pulled back to [-1, 1]
q_hat = chebfun(@(t) (L/2)^2 * q(a + (L/2) * (t + 1)), [-1 1]);

phi = cell(N, 1);
for i = 1:N
    phi{i} = legpoly(ks(i)) - legpoly(ks(i) + 2);
end

% Potential matrix: integral of q_hat phi_i phi_j, obtained from the Legendre
% coefficients c of q_hat*phi_j using orthogonality of the P_k.
n_coeffs = N + 3;
M = zeros(N);
for j = 1:N
    c = legcoeffs(chebfun(q_hat .* phi{j}, n_coeffs));

    for i = 1:N
        k = ks(i);
        ck = 0;
        if k + 1 <= numel(c)
            ck = c(k+1);
        end
        ck2 = 0;
        if k + 3 <= numel(c)
            ck2 = c(k+3);
        end
        M(i, j) = ck * (2/(2*k+1)) - ck2 * (2/(2*k+5));
    end
end

% Enforce symmetry, lost to the truncation of the coefficient expansion
M = triu(M) + triu(M, 1)';

mus = sort(real(eig(S + M, T)));

lambdas = mus / (L/2)^2;

end
