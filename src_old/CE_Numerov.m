clear all; clc;

N = 501;
beta = 30;

a = -pi/2;
b =  pi/2;

h = (b - a) / (N + 1);
x = a + h * (1:N)';

start_cpu = cputime;

q = beta^2 * sin(2*x).^2 - 2*beta * cos(2*x);

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
% fid = fopen('CoffeyEvans_Numerov_correction_501.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== NUMEROV METHOD (Coffey-Evans) (WITH CORRECTION) ===\n');
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
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambda_corrected(i));
% end
% 
% fclose(fid);