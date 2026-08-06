function eigenvalues = mbcdm_slp_example2_chebfun(N, alpha)

p = chebfunpref();
p.cheb2Prefs.chebfuneps = 1e-15;
p.cheb2Prefs.plotting = 'off';
p.bvpTol = 1e-12;

s = chebpts(N+2);

asin_alpha = asin(alpha);
g = asin(alpha * s) / asin_alpha;
x_mapped = (g + 1) * (pi / 2);

D_s = diffmat(N+2, 1);

gp = alpha ./ sqrt(1 - (alpha * s).^2) / asin_alpha;

gpp = (alpha^3 * s) ./ (1 - alpha^2 * s.^2).^(3/2) / asin_alpha;

scale = pi / 2;
dx_ds = scale * gp;
d2x_ds2 = scale * gpp;

diag_inv_dx_ds2 = diag(1 ./ (dx_ds.^2));
diag_d2x_over_dx3 = diag(d2x_ds2 ./ (dx_ds.^3));
D_s2 = diffmat(N+2, 2);
D2 = diag_inv_dx_ds2 * D_s2 - diag_d2x_over_dx3 * D_s;

Q = diag(exp(x_mapped));
A = -D2 + Q;
A_int = A(2:N+1, 2:N+1);

lambda = eig(A_int);

eigenvalues = sort(real(lambda));

end

start_cpu = cputime;

eigs = mbcdm_slp_example2_chebfun(1000, 0.999);

end_cpu = cputime;
cpu_s = end_cpu - start_cpu;

% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('SLP_exp_MBCDM_1000.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== MAPED BARYCENTRIC CHEBYSHEV DIFFERENTIATION MATRIX METHOD ===\n');
% fprintf(fid, 'N                    = %d\n', N);
% fprintf(fid, 'alpha                = %d\n', 0.999);
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
%     fprintf(fid, '%6d     %20.10f\n', k(i), eigs(i));
% end
% 
% fclose(fid);