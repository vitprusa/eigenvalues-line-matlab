clear all; clc

function lambdas = LGCC_Chebfun_Example31(N)
    m = N - 1;
    ks = (0:m-1)';

    S = diag(4*ks + 6);

    T = zeros(m);
    for i = 1:m
        k = ks(i);
        T(i,i) = 2/(2*k+1) + 2/(2*k+5);
        if i+2 <= m
            T(i,i+2) = -2/(2*k+5);
            T(i+2,i) = -2/(2*k+5);
        end
    end

    % transformed potential to interval [-1,1]
    q = chebfun(@(t) (pi^2/4) .* exp((pi/2).*(t + 1)), [-1 1]);

    phi = cell(m,1);
    for i = 1:m
        k = ks(i);
        phi{i} = legpoly(k) - legpoly(k+2);
    end

    M = zeros(m);
    for j = 1:m
        fj = q .* phi{j};
        I_fj = chebfun(fj, N+1);
        c = legcoeffs(I_fj);

        for i = 1:m
            k = ks(i);
            ck  = c(k+1);
            ck2 = (k+2 <= N) * c(k+3);
            M(i,j) = ck * (2/(2*k+1)) - ck2 * (2/(2*k+5));
        end
    end

    hatM = triu(M) + triu(M,1)';

    A = S + hatM;
    mus = sort(real(eig(A, T)));

    lambdas = 4 * mus / pi^2;
end

N = 501;

start_cpu = cputime;

lambdas = LGCC_Chebfun_Example31(N);

end_cpu = cputime;
cpu_s = end_cpu - start_cpu;

% % ==============================================================
% % === SAVE TO A TEXT FILE ===
% % ==============================================================
% 
% k = (1:N)';
% 
% fid = fopen('SLP_exp_LGCC_500.txt', 'w');
% 
% % Title with parameters
% fprintf(fid, '=== LGCC METHOD ===\n');
% fprintf(fid, 'N                    = %d\n', N);
% fprintf(fid, 'Interval             = [0, π]\n');
% fprintf(fid, 'q(x)                 = exp(x)\n');
% fprintf(fid, 'Boundary conditions  = y(0) = y(π) = 0\n');
% fprintf(fid, 'Number of eigenvalues = %d\n\n', N-1);
% 
% % Column title + separator line
% fprintf(fid, '%6s     %20s\n', 'k', 'lambda');
% fprintf(fid, '%6s     %20s\n', '------', '--------------------');
% 
% % Data
% for i = 1:N
%     fprintf(fid, '%6d     %20.10f\n', k(i), lambdas(i));
% end
% 
% fclose(fid);