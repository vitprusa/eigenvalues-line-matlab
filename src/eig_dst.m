function lambdas = eig_dst(N, a, b, q)
% Discrete sine transform method.
%
% Returns the N smallest eigenvalues of -y'' + q(x)y on [a, b] with Dirichlet
% boundary conditions y(a) = y(b) = 0, discretised on the uniform interior
% grid of N points using a DST of size N.
%
%   lambdas = eig_dst(1000, 0, pi, @q_paine)
%
% Requires dst/idst.
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

% The matrix is symmetric in exact arithmetic: the DST-I matrix is itself
% symmetric and is its own inverse up to a scalar, so the second-derivative part
% is S*Lambda*S, and the multiplication by q is diagonal.  The round-off of the
% FFT leaves it a few tens of ulp short of that -- 1.2e-9 on a norm of 2.5e5 at
% N = 500, a relative 5e-15 -- and issymmetric asks for bit-for-bit equality, so
% eig would take the general path.  Symmetrising restores the property the
% operator has and lets eig take the symmetric path, for the same spectrum to
% round-off: a relative 1e-14 at N = 1000.  The eig call itself is then about
% eight times faster; end to end this function gains rather less, two to three
% times, the assembly above dominating once the solve is cheap.
A = (A + A.') / 2;

lambdas = sort(eig(A));

end
