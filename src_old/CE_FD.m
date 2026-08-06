clear all; clc

N = 10001;           
beta = 30;

a = -pi/2;
b =  pi/2;

h = (b - a) / (N + 1);
x = a + h * (1:N)';

start_cpu = cputime;

D2 = (1/h^2) * spdiags([ 1  -2  1 ], -1:1, N, N);

q = -2*beta * cos(2*x) + beta^2 * sin(2*x).^2;
Q = spdiags(q, 0, N, N);

A = -full(D2) + full(Q);

lambda_FD = sort(eig(A));

end_cpu = cputime;
cpu_s = end_cpu - start_cpu;

% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('CoffeyEvans_FD_10001.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== FINITE DIFFERENCE METHOD (Coffey-Evans) ===\n');
% fprintf(fid, 'N                    = %d\n', N);
% fprintf(fid, 'beta                 = %d\n', beta);
% fprintf(fid, 'Interval             = [-π/2, π/2]\n');
% fprintf(fid, 'q(x)                 = -2β cos(2x) + β² sin²(2x)\n');
% fprintf(fid, 'Boundary conditions  = y(-π/2) = y(π/2) = 0\n');
% fprintf(fid, 'Number of eigenvalues = %d\n\n', N);
% 
% % Column title + separator line
% fprintf(fid, '%6s     %20s\n', 'k', 'lambda');
% fprintf(fid, '%6s     %20s\n', '------', '--------------------');
% 
% % Data
% for i = 1:N
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambda_FD(i));
% end
% 
% fclose(fid);


