clear all; clc

N = 10000;

h = pi / (N + 1);
x = h * (1:N)';

start_cpu_FD = cputime;

D2_FD = (1 / h^2) * spdiags([1 -2 1],-1:1,N,N);
Exp = diag(exp(x));

A_FD = - D2_FD + Exp;

lambda_FD = eig(A_FD);

end_cpu_FD = cputime;
cpu_s_FD = end_cpu_FD - start_cpu_FD;

% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
%
% lambda_FD = sort(real(lambda_FD));
% 
% k = (1:N)';
% 
% fid = fopen('SLP_exp_FD_10000.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== FINITE DIFFERENCE METHOD ===\n');
% fprintf(fid, 'N                    = %d\n', N);
% fprintf(fid, 'Interval             = [0, π]\n');
% fprintf(fid, 'q(x)                 = exp(x)\n');
% fprintf(fid, 'Boundary conditions  = y(0) = y(π) = 0\n');
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

