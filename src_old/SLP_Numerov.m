clear all; clc;

N = 1000;
h = pi / (N + 1);
x = h * (1:N)';

start_cpu = cputime;

q = exp(x);

main_A  = 24 + 10 * h^2 * q;
off = -(12 - h^2 * q(1:end));

A = spdiags([off, main_A, off], [-1 0 1], N, N);

B = h^2 * spdiags([1 10 1], [-1 0 1], N, N);

lambda_numerov = sort(real(eig(full(A), full(B))));

% asymptotic correction
k = (1:N)';
s2 = sin(k*h/2).^2;
eps2 = k.^2 - 12 * s2 ./ (h^2 * (3 - s2));

lambda_corrected = lambda_numerov + eps2;

end_cpu = cputime;
cpu_s = end_cpu - start_cpu;


% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('SLP_exp_Numerov_correction_1000.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== NUMEROV METHOD (WITH CORRECTION) ===\n');
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
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambda_corrected(i));
% end
% 
% fclose(fid);