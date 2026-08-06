clear all; clc

% Specify function q(x)
q = @(x) exp(x);
% Specify number of eigenvalues/matrix size
mm = 1000;

% Call dst_eig and get the eigenvalues
start_cpu_DST = cputime;

lambdas = dst_eig(mm, q);

end_cpu_DST = cputime;
cpu_s_DST = end_cpu_DST - start_cpu_DST;
fprintf('CPU time pro DST: %.4f sekund\n', cpu_s_DST)

% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:mm)';
% 
% fid = fopen('SLP_exp_DST_500.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== DISCRETE SINE TRANSFORM METHOD ===\n');
% fprintf(fid, 'N                    = %d\n', mm);
% fprintf(fid, 'Interval             = [0, π]\n');
% fprintf(fid, 'q(x)                 = exp(x)\n');
% fprintf(fid, 'Boundary conditions  = y(0) = y(π) = 0\n');
% fprintf(fid, 'Number of eigenvalues = %d\n\n', mm);
% 
% % Column title + separator line
% fprintf(fid, '%6s     %20s\n', 'k', 'lambda');
% fprintf(fid, '%6s     %20s\n', '------', '--------------------');
% 
% % Data
% for i = 1:mm
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambdas(i));
% end
% 
% fclose(fid);