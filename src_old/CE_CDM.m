clear all; clc

N = 1002;  
beta = 30;            

start_cpu_CDM = cputime;

M = N + 1;
t = chebpts(M, 2);

D2_full = diffmat(M, 2);
D2_Cheb = D2_full(2:end-1, 2:end-1);

t_inter = t(2:end-1);
x_inter = (pi/2) * t_inter;
            
q = -2*beta * cos(2*x_inter) + beta^2 * sin(2*x_inter).^2;
Q = diag(q);
scale = (2 / pi)^2;
D2_Cheb_scaled = scale * D2_Cheb;
A_Cheb = -D2_Cheb_scaled + Q;
 
lambda_Cheb = sort(real(eig(A_Cheb)));

end_cpu_CDM = cputime;
cpu_s_CDM = end_cpu_CDM - start_cpu_CDM;


% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('CoffeyEvans_CDM_1001.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== CHEBYSHEV DIFFERENTIATION MATRIX METHOD (Coffey-Evans) ===\n');
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
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambda_Cheb(i));
% end
% 
% fclose(fid);

