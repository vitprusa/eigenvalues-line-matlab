clear all; clc
N = 1000;

start_cpu_CDM = cputime;                 

M = N + 2;
t = chebpts(M, 2);                 

D2_full = diffmat(M, 2);

D2_Cheb = D2_full(2:end-1, 2:end-1);

t_inter = t(2:end-1);
x_inter = (pi / 2) * (t_inter + 1);

scale = (2 / pi)^2;
D2_Cheb_scaled = scale * D2_Cheb;

Exp = diag(exp(x_inter));

A_Cheb = -D2_Cheb_scaled + Exp;

lambda_Cheb = sort(eig(A_Cheb));

end_cpu_CDM = cputime;
cpu_s_CDM = end_cpu_CDM - start_cpu_CDM;


% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('SLP_exp_CDM_1000.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== CHEBYSHEV DIFFERENTIATION MATRIX METHOD ===\n');
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
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambda_Cheb(i));
% end
% 
% fclose(fid);