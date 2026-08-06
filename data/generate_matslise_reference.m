function generate_matslise_reference(n_eig, tol)
% Compute reference eigenvalues with MATSLISE and write them to CSV.
%
% Produces  paine_matslise.csv  and  coffey_evans_matslise.csv  in the
% directory holding this file.  Both problems use the same potentials as the
% solvers in ../src, so the reference and the methods under test cannot drift
% apart.
%
%   generate_matslise_reference          % 500 eigenvalues, tol = 1e-12
%   generate_matslise_reference(500, 1e-10)
arguments
    % number of eigenvalues
    n_eig (1, 1) {mustBeInteger, mustBePositive} = 500
    % input tolerance passed to computeEigenvalues
    tol (1, 1) {mustBeReal, mustBePositive} = 1e-12
end

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
matslise_dir = fullfile(fileparts(root), 'MATSLISE2');

addpath(genpath(fullfile(matslise_dir, 'source')));
addpath(fullfile(root, 'src'));

one = @(x) ones(size(x));       % p and w for -y'' + q(x)y = lambda y

% Throwaway call so that the first timed case does not also pay for JIT
% compilation and class loading.
computeEigenvalues(slp(one, @q_paine, one, 0, pi, 1, 0, 1, 0), 0, 2, 1e-10, true);

% MATSLISE indexes eigenvalues from zero, the article from one, so the last
% index requested is n_eig - 1.
kmax = n_eig - 1;

run_case(fullfile(here, 'paine_matslise.csv'), ...
    'Paine problem 1', ...
    'q(x) = exp(x)', ...
    0, pi, ...
    'slp(one, @q_paine, one, 0, pi, 1, 0, 1, 0)', ...
    slp(one, @q_paine, one, 0, pi, 1, 0, 1, 0), ...
    n_eig, kmax, tol, matslise_dir);

run_case(fullfile(here, 'coffey_evans_matslise.csv'), ...
    'Coffey--Evans problem, beta = 30 (default)', ...
    'q(x) = -2*beta*cos(2x) + beta^2*sin(2x)^2,  beta = 30', ...
    -pi/2, pi/2, ...
    'slp(one, @q_coffey_evans, one, -pi/2, pi/2, 1, 0, 1, 0)', ...
    slp(one, @q_coffey_evans, one, -pi/2, pi/2, 1, 0, 1, 0), ...
    n_eig, kmax, tol, matslise_dir);

end


function run_case(csv_path, title, q_text, a, b, ctor_text, o, n_eig, kmax, tol, matslise_dir)

call_text = sprintf('computeEigenvalues(o, 0, %d, %g, true)', kmax, tol);

fprintf('%s: computing %d eigenvalues ...\n', title, n_eig);

t = tic;
E = computeEigenvalues(o, 0, kmax, tol, true);
elapsed = toc(t);

if ~E.success
    error('MATSLISE failed for %s: %s', title, E.msg);
end

fprintf('  done in %.3f s, %d eigenvalues, max status flag %d\n', ...
        elapsed, numel(E.eigenvalues), max(E.status));

n = numel(E.eigenvalues);

fid = fopen(csv_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# %s\n', title);
fprintf(fid, '# Reference eigenvalues of the Sturm-Liouville problem\n');
fprintf(fid, '#     -y'''' + q(x) y = lambda y,   y(a) = y(b) = 0\n');
fprintf(fid, '#     %s\n', q_text);
fprintf(fid, '#     [a, b] = [%.15g, %.15g]\n', a, b);
fprintf(fid, '#\n');
fprintf(fid, '# METHOD\n');
fprintf(fid, '#   MATSLISE 2, CP methods: piecewise constant reference potential plus\n');
fprintf(fid, '#   perturbation corrections, on an automatically constructed mesh.  The\n');
fprintf(fid, '#   requested input tolerance was %g.\n', tol);
fprintf(fid, '#   MATSLISE directory: %s\n', matslise_dir);
fprintf(fid, '#\n');
fprintf(fid, '# SOFTWARE REFERENCE (Chicago style)\n');
fprintf(fid, '#   Ledoux, Veerle, and Marnix Van Daele. "Matslise 2.0: A Matlab Toolbox\n');
fprintf(fid, '#     for Sturm-Liouville Computations." ACM Transactions on Mathematical\n');
fprintf(fid, '#     Software 42, no. 4 (2016): 29:1-29:18.\n');
fprintf(fid, '#     https://doi.org/10.1145/2839299.\n');
fprintf(fid, '#   Ledoux, Veerle, Marnix Van Daele, and Guido Vanden Berghe. "MATSLISE:\n');
fprintf(fid, '#     A MATLAB Package for the Numerical Solution of Sturm-Liouville and\n');
fprintf(fid, '#     Schrodinger Equations." ACM Transactions on Mathematical Software 31,\n');
fprintf(fid, '#     no. 4 (2005): 532-554. https://doi.org/10.1145/1114268.1114273.\n');
fprintf(fid, '#   Ledoux, Veerle, and Marnix Van Daele. MATSLISE. Version 2. MATLAB\n');
fprintf(fid, '#     package. SourceForge, 2015. https://sourceforge.net/projects/matslise/.\n');
fprintf(fid, '#\n');
fprintf(fid, '# COMMAND\n');
fprintf(fid, '#   The potentials are the function handles used by the solvers in ../src,\n');
fprintf(fid, '#   so the reference matches the methods under test exactly.\n');
fprintf(fid, '#     addpath(genpath(fullfile(matslise_dir, ''source'')));\n');
fprintf(fid, '#     addpath(fullfile(root, ''src''));\n');
fprintf(fid, '#     one = @(x) ones(size(x));\n');
fprintf(fid, '#     o   = %s;\n', ctor_text);
fprintf(fid, '#     E   = %s;\n', call_text);
fprintf(fid, '#   Run in full, from this directory, by:\n');
fprintf(fid, '#     matlab -batch "generate_matslise_reference(%d, %g)"\n', n_eig, tol);
fprintf(fid, '#\n');
fprintf(fid, '# TIMING\n');
fprintf(fid, '#   WALL_CLOCK = %.6f s\n', elapsed);
fprintf(fid, '#   WALL_CLOCK_PER_EIGENVALUE = %.6g s\n', elapsed / n);
fprintf(fid, '#   Wall clock for the computeEigenvalues call above, mesh construction\n');
fprintf(fid, '#   included; one call returning all %d eigenvalues.  A throwaway call\n', n);
fprintf(fid, '#   precedes the timed ones, so this excludes JIT compilation and class\n');
fprintf(fid, '#   loading (warm start).\n');
fprintf(fid, '#\n');
fprintf(fid, '# ENVIRONMENT\n');
fprintf(fid, '#   MATLAB %s on %s\n', version, computer);
fprintf(fid, '#\n');
fprintf(fid, '# COLUMNS\n');
fprintf(fid, '#   k                index counting from one (convention of the article)\n');
fprintf(fid, '#   matslise_index   index counting from zero, as returned by MATSLISE\n');
fprintf(fid, '#   eigenvalue       lambda_k\n');
fprintf(fid, '#   estimated_error  MATSLISE error estimate for that eigenvalue\n');
fprintf(fid, '#   status           MATSLISE status flag, 0 = correct calculation\n');
fprintf(fid, '#\n');

fprintf(fid, 'k,matslise_index,eigenvalue,estimated_error,status\n');
for i = 1:n
    fprintf(fid, '%d,%d,%.15g,%.6e,%d\n', ...
            E.indices(i) + 1, E.indices(i), E.eigenvalues(i), E.errors(i), E.status(i));
end

fprintf('  written %s\n', csv_path);

end
