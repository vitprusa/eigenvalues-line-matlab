function lambdas = eig_numerov_corrected(N, a, b, q)
% Numerov method with the asymptotic correction.
%
% Returns the N eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet boundary
% conditions y(a) = y(b) = 0, computed by eig_numerov and then corrected in
% the manner of Paine, de Hoog and Anderssen: the k-th Numerov eigenvalue is
% shifted by the error the scheme commits on the same problem with q = 0,
%
%   lambda_k <- lambda_k + (k*pi/(b-a))^2 - 12*s^2/(h^2*(3 - s^2)),
%   s = sin(k*pi/(2(N+1))),
%
% which removes the leading discretisation error for the higher modes.
%
%   lambdas = eig_numerov_corrected(1000, 0, pi, @q_paine)
arguments
    % matrix size / number of eigenvalues
    N (1, 1) {mustBeInteger, mustBePositive}
    % endpoints of the interval
    a (1, 1) {mustBeReal}
    b (1, 1) {mustBeReal, mustBeGreaterThan(b, a)}
    % potential q
    q function_handle
end

L = b - a;
h = L / (N + 1);

k = (1:N)';

% Eigenvalues of the q = 0 problem, exact and as seen by the Numerov scheme
s2 = sin(k * pi / (2*(N + 1))).^2;
lambda_exact = (k * pi / L).^2;
lambda_discrete = 12 * s2 ./ (h^2 * (3 - s2));

lambdas = eig_numerov(N, a, b, q) + lambda_exact - lambda_discrete;

end
