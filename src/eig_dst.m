function lambdas = eig_dst(N, a, b, q)
% Discrete sine transform method.
%
% Returns the N smallest eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet
% boundary conditions y(a) = y(b) = 0, discretised on the uniform interior
% grid of N points using a DST of size N.
%
%   lambdas = eig_dst(1000, 0, pi, @q_paine)
%
% Requires dst/idst (Signal Processing Toolbox).
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
x = a + h * (1:N)';

% Eigenvalues of -d^2/dx^2 on the sine basis sin(k*pi*(x-a)/L): differentiation
% is diagonal in DST space, so -y'' is idst(kvec .* dst(y)).
kvec = (pi * (1:N)' / L).^2;

% Apply the operator to the identity to form the matrix column by column
% (dst acts column-wise), then add the multiplication operator by q.
A = idst(kvec .* dst(eye(N))) + diag(q(x));

lambdas = sort(real(eig(A)));

end
