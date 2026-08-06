function lambdas = dst_eig(mm, q)
% Returns mm eigenvalues of linear operator -y'' + q(x)y on interval [0, pi]
% with Dirichlet boundary conditions, operator is discretised using DST
% based technique with matrix size mm.
arguments
    % size of matrix/number of eigenvalues
    mm (1, 1) {mustBeInteger, mustBePositive}
    % function q
    q function_handle
end

% Second derivative out of vector of grid values x, operator form
% Grid values are arranged in a column vector
% Symbol kvec denotes the vector of negative squares, differentiation in DST space
kvec = @(x) -((1:length(x)).^2);
D2op = @(x) idst(kvec(x)' .* dst(x));

% Multiplication of grid values vector x with grid values q(x) function
% Grid values are arranged in a column vector

xgrid = @(x) pi/(length(x)+1) * (1:length(x));             % pro [0,pi]   
%xgrid = @(x) -pi/2 + (pi/(length(x)+1)) * (1:length(x))';    % pro [-pi/2,pi/2]
Qop = @(x) (q(xgrid(x)') .* x);

% Sturm--Liouville operator -y'' + q(x)y
% Grid values are arranged in a column vector -- we need this for eigs
% solver in the matrix--free form
Lop = @(x) -D2op(x) + Qop(x);

% Solve the eigenvalue problem using explicitely formed matrix, return
% eigenvalues
L = Lop(eye(mm));
[~, D] = eig(L);
lambdas = sort(diag(D));

end