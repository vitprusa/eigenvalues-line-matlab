clear all; clc

function lambdas = LGCC_CoffeyEvans(N, beta)
    if nargin < 2, beta = 30; end

    m = N - 1;
    ks = (0:m-1)';

    S = diag(4*ks + 6);

    T = zeros(m);
    for i = 1:m
        k = ks(i);
        T(i,i) = 2/(2*k+1) + 2/(2*k+5);
        if i+2 <= m
            val = -2/(2*k+5);
            T(i, i+2) = val;
            T(i+2, i) = val;
        end
    end

    % transformed potential to interval [-1,1]
    q_new = chebfun(@(t) (pi^2/4) * (-2*beta*cos(pi*t) + beta^2*sin(pi*t).^2), [-1 1]);

    phi = cell(m,1);
    for i = 1:m
        k = ks(i);
        phi{i} = legpoly(k) - legpoly(k+2);
    end

    M = zeros(m);
    for j = 1:m
        fj = q_new .* phi{j};
        I_fj = chebfun(fj, N+1);
        c = legcoeffs(I_fj);

        for i = 1:m
            k = ks(i);
            ck  = c(k+1);
            ck2 = (k+3 <= length(c)) * c(k+3);
            M(i,j) = ck * (2/(2*k+1)) - ck2 * (2/(2*k+5));
        end
    end

    hatM = triu(M) + triu(M,1)';

    A = S + hatM;
    [~, D] = eig(A, T);
    mus = sort(real(diag(D)));

    lambdas = 4 * mus / pi^2;
end

beta = 30;
N = 1002;                              

start_cpu = cputime;

lambdas_lgcc = LGCC_CoffeyEvans(N, beta);

end_cpu = cputime;
cpu_s = end_cpu - start_cpu;


% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('CoffeyEvans_LGCC_1001.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== LGCC METHOD (Coffey-Evans) ===\n');
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
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambdas_lgcc(i));
% end
% 
% fclose(fid);