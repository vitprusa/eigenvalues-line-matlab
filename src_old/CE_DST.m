clear all; clc

% Specify function q(x)
beta = 30;
q = @(x) -2*beta * cos(2*x) + beta^2 * sin(2*x).^2;
% Specify number of eigenvalues/matrix size
mm = 1001;

% Call dst_eig and get the eigenvalues
start_cpu_DST = cputime;

lambdas = dst_eig(mm, q);

end_cpu_DST = cputime;
cpu_s_DST = end_cpu_DST - start_cpu_DST;

% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:mm)';
% 
% fid = fopen('CoffeyEvans_DST_1001.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== FINITE DIFFERENCE METHOD (Coffey-Evans) ===\n');
% fprintf(fid, 'N                    = %d\n', mm);
% fprintf(fid, 'beta                 = %d\n', beta);
% fprintf(fid, 'Interval             = [-π/2, π/2]\n');
% fprintf(fid, 'q(x)                 = -2β cos(2x) + β² sin²(2x)\n');
% fprintf(fid, 'Boundary conditions  = y(-π/2) = y(π/2) = 0\n');
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