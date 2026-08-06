function eigenvalues = mbcdm_coffey_evans(N, beta, alpha)

p = chebfunpref();
p.cheb2Prefs.chebfuneps = 1e-14;
p.cheb2Prefs.plotting   = 'off';
p.bvpTol = 1e-12;

s = chebpts(N+1);

asin_alpha = asin(alpha);
g  = asin(alpha * s) / asin_alpha;
x  = (pi/2) * g;

gp  = alpha ./ sqrt(1 - (alpha * s).^2) / asin_alpha;
gpp = (alpha^3 * s) ./ (1 - alpha^2 * s.^2).^(3/2) / asin_alpha;

dx_ds  = (pi/2) * gp;
d2x_ds2 = (pi/2) * gpp;

D_s  = diffmat(N+1, 1);
D_s2 = diffmat(N+1, 2);

inv_dx_ds2       = 1 ./ (dx_ds .^ 2);
d2x_over_dx3     = d2x_ds2 ./ (dx_ds .^ 3);

D2 = diag(inv_dx_ds2) * D_s2  -  diag(d2x_over_dx3) * D_s;

q_vals = -2*beta * cos(2*x) + beta^2 * sin(2*x).^2;
Q = diag(q_vals);

A = -D2 + Q;

A_int = A(2:end-1, 2:end-1);

lambda = eig(A_int);
eigenvalues = sort(real(lambda));

end

start_cpu = cputime;

N = 1001;
eigenval = mbcdm_coffey_evans(N,30,1);

end_cpu = cputime;
cpu_s = end_cpu - start_cpu;


% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('CoffeyEvans_MBCDM_1001.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== MAPED BARYCENTRIC CHEBYSHEV DIFFERENTIATION MATRIX METHOD (Coffey-Evans) ===\n');
% fprintf(fid, 'N                    = %d\n', N);
% fprintf(fid, 'alpha                 = %d\n', 1);
% fprintf(fid, 'beta                 = %d\n', 30);
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
%     fprintf(fid, '%6d     %20.10f\n', k(i), eigenval(i));
% end
% 
% fclose(fid);