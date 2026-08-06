function make_eigenvalues_head_paper_tables(name)
%MAKE_EIGENVALUES_HEAD_PAPER_TABLES Paper eigenvalue-comparison tables (self-computed).
%
%   Line (one-dimensional Sturm-Liouville) counterpart of the plane project's
%   experiments_paper/eigenvalues_head_paper. It GENERATES ITS OWN DATA: for
%   each test problem it recomputes the spectrum of
%
%       -y'' + q(x) y = lambda y,   y(a) = y(b) = 0
%
%   with every method in src/ -- DST, finite differences, Chebyshev
%   differentiation matrix, mapped barycentric Chebyshev, Numerov, Numerov with
%   the asymptotic correction and Legendre-Galerkin -- at N = 500 and N = 1000,
%   finite differences additionally at N = 10000, timing every run with tic/toc.
%
%   For each problem it writes
%     - one CSV per (method, N) into results_paper/eigenvalues_head/cache/, each
%       carrying the DOF count and the measured time in its metadata header.
%       These are the cached spectra: a run whose CSV is present is read back
%       rather than recomputed, so a layout-only regeneration costs no
%       computation at all. Delete a CSV to compute that run again.
%     - one LaTeX snippet <problem>_eigenvalues_head_transposed.tex into
%       results_paper/eigenvalues_head/, holding a booktabs table in the
%       transposed layout: one column per run (the MATSLISE reference plus every
%       method and size), one row per eigenvalue index, so that a single
%       eigenvalue reads across all discretisations. The DOF and timing are the
%       two leading body rows. It carries the leading eight eigenvalues and then
%       lambda_50, lambda_100, lambda_400 and lambda_500, each on an empty row
%       of its own.
%   The table is \scriptsize; the digits are set per problem by FMT_MODE and
%   FMT_N (see PROBLEM_CONFIGS). Wider than the text width, it is scaled down
%   with \resizebox.
%
%   The reference is the MATSLISE spectrum in data/, read and not recomputed;
%   see data/README.md. It holds 500 eigenvalues per problem, so every index
%   shown here has a reference value. DOF is the size of the matrix actually
%   solved, which is N for every method (the Chebyshev grids carry N + 2 points,
%   two of which the Dirichlet conditions remove).
%
%   Called with no argument (or an empty one) it does both problems; pass a
%   problem name (paine, coffey_evans) to regenerate just that one.
%
%   Requires dst/idst and Chebfun.

    NEIG  = 8;                       % leading eigenvalues shown per table
    EXTRA = [50 100 400 500];        % deeper eigenvalues of the extended table

    here         = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(here));
    addpath(fullfile(project_root, 'src'));

    out_dir   = fullfile(project_root, 'results_paper', 'eigenvalues_head');
    cache_dir = fullfile(out_dir, 'cache');
    if ~exist(cache_dir, 'dir')
        mkdir(cache_dir);
    end

    cfgs = problem_configs(project_root);
    if nargin >= 1 && ~isempty(name)
        cfgs = cfgs(strcmp({cfgs.name}, name));
        if isempty(cfgs)
            error('head_paper:unknownProblem', ...
                  'Unknown problem "%s" (known: paine, coffey_evans).', name);
        end
    end

    for i = 1:numel(cfgs)
        cfg = cfgs(i);
        fprintf('=== %s ===\n', cfg.name);
        rows = compute_problem(cfg, cache_dir);

        tex = fullfile(out_dir, sprintf('%s_eigenvalues_head_transposed.tex', cfg.name));
        write_latex_transposed(tex, cfg, rows, NEIG, EXTRA);
        fprintf('Wrote %s\n', tex);
    end
end


function cfgs = problem_configs(project_root)
%PROBLEM_CONFIGS Per-problem interval, potential, method schedule and reference.
%
%   The method schedule is shared by both problems: every method at N = 500 and
%   N = 1000, finite differences additionally at N = 10000. LABEL is what the
%   tables print, KEY what the file names use.
%
%   FMT_MODE and FMT_N set how many digits the eigenvalue cells carry, per
%   problem: 'decimals' prints FMT_N digits after the point, 'significant'
%   prints FMT_N significant digits by varying the decimals with the magnitude.
    methods = struct( ...
        'label', {'DST', 'FD', 'CDM', 'MBCDM', 'Numerov', 'Numerov+AC', 'LGCC'}, ...
        'key',   {'dst', 'fd', 'chebyshev', 'chebyshev_mapped', ...
                  'numerov', 'numerov_corrected', 'legendre_galerkin'}, ...
        'fun',   {@eig_dst, @eig_fd, @eig_chebyshev, @eig_chebyshev_mapped, ...
                  @eig_numerov, @eig_numerov_corrected, @eig_legendre_galerkin}, ...
        'N',     {[500 1000], [500 1000 10000], [500 1000], [500 1000], ...
                  [500 1000], [500 1000], [500 1000]});

    cfgs = struct('name', {}, 'pretty', {}, 'q_text', {}, 'a', {}, 'b', {}, ...
                  'q', {}, 'methods', {}, 'reference_csv', {}, ...
                  'fmt_mode', {}, 'fmt_n', {});

    cfgs(end+1) = struct( ...
        'name', 'paine', 'pretty', 'Paine', ...
        'q_text', '$q(x) = \mathrm{e}^{x}$', ...
        'a', 0, 'b', pi, 'q', @q_paine, 'methods', methods, ...
        'reference_csv', fullfile(project_root, 'data', 'paine_matslise.csv'), ...
        'fmt_mode', 'significant', 'fmt_n', 6);

    cfgs(end+1) = struct( ...
        'name', 'coffey_evans', 'pretty', 'Coffey--Evans', ...
        'q_text', ['$q(x) = -2\beta\cos 2x + \beta^{2}\sin^{2} 2x$, ', ...
                   '$\beta = 30$'], ...
        'a', -pi/2, 'b', pi/2, 'q', @q_coffey_evans, 'methods', methods, ...
        'reference_csv', fullfile(project_root, 'data', 'coffey_evans_matslise.csv'), ...
        'fmt_mode', 'significant', 'fmt_n', 6);
end


function rows = compute_problem(cfg, cache_dir)
%COMPUTE_PROBLEM Compute every method/run for one problem; cache CSVs; collect rows.
%
%   Each row keeps the whole spectrum, not just the leading values: the extended
%   table reads deeper indices out of the same rows, and the writers mark
%   anything past the end of a row with "---".
    rows = struct('method_label', {}, 'group', {}, 'dof', {}, 'time', {}, ...
                  'eigs', {}, 'bold', {});

    [ref_evals, ref_time] = read_matslise_csv(cfg.reference_csv);
    rows(end+1) = mk('\texttt{MATSLISE}', 'truth', '---', ref_time, ref_evals, true);

    for mi = 1:numel(cfg.methods)
        m = cfg.methods(mi);
        for k = 1:numel(m.N)
            N   = m.N(k);
            csv = fullfile(cache_dir, sprintf('%s_%s_N%d-eigenvalues.csv', ...
                           cfg.name, m.key, N));
            if exist(csv, 'file')
                % Cached: read the eigenvalues, DOF count and measured time
                % back, no recompute.
                [evals, dofs, tsec] = read_head_csv(csv);
                fprintf('  %-11s N = %-5d  (cached)\n', m.label, N);
            else
                try
                    [evals, dofs, tsec] = run_one(m, N, cfg);
                catch ME
                    warning('head_paper:run', '%s %s N = %d failed: %s', ...
                            cfg.name, m.label, N, ME.message);
                    continue;
                end
                write_head_csv(csv, cfg, m.label, N, dofs, tsec, evals);
                fprintf('  %-11s N = %-5d  dofs = %-6d  time = %8.2f s  ->  %s\n', ...
                        m.label, N, dofs, tsec, csv);
            end
            rows(end+1) = mk(m.label, m.label, num2str(dofs), ...
                             sprintf('%.2f', tsec), evals, false); %#ok<AGROW>
        end
    end
end


function [evals, dofs, tsec] = run_one(m, N, cfg)
%RUN_ONE One timed spectrum computation for method M at size N.
%
%   The timing covers the solver call only, and never a cold call: warm-up is
%   spent on discarded runs by ENSURE_WARM instead. Left in, it lands on
%   whichever run happens to go first and makes a coarse grid look slower than a
%   finer one.
    ensure_warm(m, cfg);
    t = tic;
    evals = m.fun(N, cfg.a, cfg.b, cfg.q);
    tsec = toc(t);
    evals = sort(real(evals(:)), 'ascend');
    dofs = numel(evals);
end


function ensure_warm(m, cfg)
%ENSURE_WARM Discarded runs of method M, before its first timing on a problem.
%
%   Warm-up covers MATLAB's JIT compilation and, for the Chebyshev and
%   Legendre-Galerkin methods, Chebfun's own load. It is done at a small fixed
%   size rather than at the first scheduled one: N = 500 of the Legendre-Galerkin
%   method costs minutes, and a warm-up is by definition thrown away.
%
%   Lazy on purpose: a run whose CSV is cached never reaches RUN_ONE, so a
%   layout-only regeneration still costs no computation at all.
    persistent warmed
    if isempty(warmed)
        warmed = struct();
    end
    key = sprintf('%s_%s', m.key, cfg.name);
    if isfield(warmed, key)
        return;
    end
    warmed.(key) = true;   % set first, so that a failed warm-up is not retried
    fprintf('  warm-up: %s on %s (discarded)\n', m.label, cfg.name);
    try
        for i = 1:2
            m.fun(40, cfg.a, cfg.b, cfg.q);
        end
    catch
        % Ignored: the real run reports its own failure through the caller.
    end
end


function s = mk(method_label, group, dof, time, eigs, bold)
    s = struct('method_label', method_label, 'group', group, 'dof', dof, ...
               'time', time, 'eigs', eigs(:)', 'bold', bold);
end


function [evals, tsec] = read_matslise_csv(path)
%READ_MATSLISE_CSV Reference eigenvalues and wall clock from a data/ CSV.
%
%   The column holding the eigenvalues is located by name from the header row,
%   the timing by the WALL_CLOCK key of the commented header.
    fid = fopen(path, 'r');
    if fid == -1
        error('head_paper:reference', ...
              ['Could not open reference CSV %s. Generate it with ', ...
               'data/generate_matslise_reference.m.'], path);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

    evals = [];
    tsec  = '---';
    col   = NaN;
    while true
        line = fgetl(fid);
        if ~ischar(line); break; end
        s = strtrim(line);
        if isempty(s); continue; end
        if s(1) == '#'
            tok = regexp(s, 'WALL_CLOCK\s*=\s*([0-9.]+)', 'tokens', 'once');
            if ~isempty(tok); tsec = sprintf('%.2f', str2double(tok{1})); end
            continue;
        end
        parts = strtrim(strsplit(s, ','));
        if isnan(col)
            col = find(strcmp(parts, 'eigenvalue'), 1);
            if isempty(col)
                error('head_paper:reference', ...
                      'No "eigenvalue" column in %s.', path);
            end
            continue;
        end
        if numel(parts) >= col
            v = str2double(parts{col});
            if ~isnan(v); evals(end+1, 1) = v; end %#ok<AGROW>
        end
    end
end


function [evals, dofs, tsec] = read_head_csv(path)
%READ_HEAD_CSV Read a paper-head CSV back: eigenvalues, DOF count, measured time.
    fid = fopen(path, 'r');
    if fid == -1
        error('head_paper:csvread', 'Could not open %s.', path);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    dofs = NaN; tsec = NaN; evals = [];
    while true
        line = fgetl(fid);
        if ~ischar(line); break; end
        s = strtrim(line);
        if isempty(s); continue; end
        if s(1) == '#'
            m = regexp(s, 'dofs\s*=\s*(\d+)', 'tokens', 'once');
            if ~isempty(m); dofs = str2double(m{1}); end
            m = regexp(s, 'Computation time:\s*([0-9.]+)', 'tokens', 'once');
            if ~isempty(m); tsec = str2double(m{1}); end
            continue;
        end
        if strncmpi(s, 'n,', 2); continue; end
        parts = strsplit(s, ',');
        if numel(parts) >= 2
            v = str2double(parts{2});
            if ~isnan(v); evals(end+1, 1) = v; end %#ok<AGROW>
        end
    end
end


function write_head_csv(csv, cfg, method, N, dofs, tsec, evals)
%WRITE_HEAD_CSV One spectrum CSV with DOF/time metadata, then n,lambda_n rows.
    fid = fopen(csv, 'w');
    if fid == -1
        error('head_paper:csv', 'Could not open %s for writing.', csv);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, ['# Problem: %s (%s, paper eigenvalues-head run, ', ...
                  'dense eig full spectrum)\n'], cfg.name, method);
    fprintf(fid, '#   -y'''' + q(x) y = lambda y,  y(a) = y(b) = 0\n');
    fprintf(fid, '#   [a, b] = [%.15g, %.15g]\n', cfg.a, cfg.b);
    fprintf(fid, '# Resolution N = %d, dofs = %d\n', N, dofs);
    fprintf(fid, '# Computation time: %.4f s\n', tsec);
    fprintf(fid, 'n,lambda_n\n');
    for i = 1:numel(evals)
        fprintf(fid, '%d,%.12g\n', i, evals(i));
    end
end


function t = caption_bits(cfg)
%CAPTION_BITS Caption and header-comment wording.
    labels = {cfg.methods.label};
    tt = cellfun(@(m) sprintf('\\texttt{%s}', m), labels, 'UniformOutput', false);
    t.method_list = sprintf('%s and %s', strjoin(tt(1:end-1), ', '), tt{end});
    t.method_list_plain = strjoin(labels, ', ');

    ns = unique([cfg.methods.N]);
    t.lev_note = sprintf('N = %s', strjoin(arrayfun(@(n) num2str(n), ns, ...
                         'UniformOutput', false), ', '));

    if abs(cfg.a) < 1e-12
        t.interval = sprintf('0, \\pi');
    else
        t.interval = '-\pi/2, \pi/2';
    end
    t.interval = sprintf('$%s$', t.interval);
    t.interval = t.interval(2:end-1);   % the caption supplies its own $[...]$
end


function write_latex_transposed(tex, cfg, rows, NEIG, extras)
%WRITE_LATEX_TRANSPOSED Transposed booktabs table: one row per eigenvalue.
%
%   One column per run (the reference plus every method and size), one row per
%   eigenvalue index, so that a single eigenvalue can be read across all
%   discretisations. The method name spans its runs as a \multicolumn group head,
%   the DOF and timing become the two leading body rows, and the reference column
%   stays bold.
%
%   EXTRAS are eigenvalue indices past the leading NEIG, each written on a row of
%   its own after an empty row, so that the jump in the index is visible. A
%   non-empty EXTRAS also wraps the tabular in \resizebox, the deeper eigenvalues
%   running into six figures before the decimal point and setting every column
%   width. Pass [] for the leading indices alone.
    fid = fopen(tex, 'w');
    if fid == -1
        error('head_paper:tex', 'Could not open %s for writing.', tex);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

    t = caption_bits(cfg);
    ncol = numel(rows);
    colspec = ['l', repmat('r', 1, ncol)];

    % Contiguous column groups: the single reference column, then one per method.
    starts = 1;
    for i = 2:ncol
        if ~strcmp(rows(i).group, rows(i - 1).group)
            starts(end + 1) = i; %#ok<AGROW>
        end
    end
    stops = [starts(2:end) - 1, ncol];

    fprintf(fid, ['%% Eigenvalue comparison for the %s problem \\texttt{%s}: ', ...
        '%s,\n'], cfg.pretty, latex_name(cfg.name), t.method_list_plain);
    fprintf(fid, ['%% %s, first %d eigenvalues, self-computed with timing.\n', ...
        '%% Transposed layout: one column per run, one row per eigenvalue.\n'], ...
        t.lev_note, NEIG);
    if ~isempty(extras)
        fprintf(fid, ['%% Extended with the eigenvalues %s;\n', ...
            '%% the tabular is scaled down if it exceeds the text width.\n'], ...
            index_list(extras));
    end
    fprintf(fid, '%%\n%% Requires in the preamble:\n');
    fprintf(fid, ['%%   \\usepackage{booktabs}\n', ...
                  '%%   \\usepackage{amsmath}\n%%   \\usepackage{listings}\n']);
    if ~isempty(extras)
        fprintf(fid, '%%   \\usepackage{graphicx}          %% \\resizebox\n');
    end
    fprintf(fid, '\n');

    fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
    % The size and \tabcolsep changes scope to the tabular (grouped) so that the
    % caption keeps the normal body size.
    fprintf(fid, '  {\\scriptsize\n  \\setlength{\\tabcolsep}{3pt}\n');
    if ~isempty(extras)
        % Scaled DOWN only: \width is the natural width of the tabular, so a
        % table that already fits is passed through at its own size rather than
        % blown up to fill the text width. Trailing %% so that the line break
        % adds no space before the tabular.
        fprintf(fid, ['  \\resizebox{\\ifdim\\width>\\textwidth\\textwidth', ...
                      '\\else\\width\\fi}{!}{%%\n']);
    end
    fprintf(fid, '  \\begin{tabular}{%s}\n    \\toprule\n', colspec);

    % Group head: method name spanning its runs, underlined by \cmidrule.
    fprintf(fid, '   ');
    for g = 1:numel(starts)
        lab = rows(starts(g)).method_label;
        if ~strcmp(rows(starts(g)).group, 'truth')
            lab = sprintf('\\texttt{%s}', lab);
        end
        fprintf(fid, ' & \\multicolumn{%d}{c}{%s}', stops(g) - starts(g) + 1, lab);
    end
    fprintf(fid, ' \\\\\n   ');
    for g = 1:numel(starts)
        fprintf(fid, ' \\cmidrule(lr){%d-%d}', starts(g) + 1, stops(g) + 1);
    end
    fprintf(fid, '\n');

    % Run metadata, then the eigenvalues.
    fprintf(fid, '    DOF');
    fprintf(fid, ' & %s', rows.dof);
    fprintf(fid, ' \\\\\n    Time (s)');
    fprintf(fid, ' & %s', rows.time);
    fprintf(fid, ' \\\\\n    \\midrule\n');

    for i = 1:NEIG
        write_eig_row(fid, rows, i, cfg);
    end
    % The deeper eigenvalues, each after an empty row: the index jumps from one
    % to the next, and a rule would read as a new block of the same sequence.
    for i = extras(:)'
        fprintf(fid, '   %s \\\\\n', repmat(' &', 1, ncol));
        write_eig_row(fid, rows, i, cfg);
    end

    if isempty(extras)
        fprintf(fid, '    \\bottomrule\n  \\end{tabular}}\n');
    else
        fprintf(fid, '    \\bottomrule\n  \\end{tabular}}}\n');   % tabular, \resizebox, size group
    end

    timing_clause = [' Timing shows the time for the matrix assembly and ', ...
        'full spectrum computation using \texttt{MATLAB}''s ', ...
        '\lstinline{eig} function with the default settings.'];
    fprintf(fid, ['  \\caption{Eigenvalues of the Sturm--Liouville operator ', ...
        '$-y'''' + q(x)y$ on $[%s]$ with %s, the %s problem \\texttt{%s}. ', ...
        'Numerical eigenvalues computed by each discretisation %s with varying ', ...
        'degrees of freedom (DOF); compared with the reference eigenvalues from ', ...
        '\\texttt{MATSLISE}.%s}\n'], ...
        t.interval, cfg.q_text, cfg.pretty, latex_name(cfg.name), ...
        t.method_list, timing_clause);
    fprintf(fid, '  \\label{tab:eigenvalues_head_%s_transposed}\n', cfg.name);
    fprintf(fid, '\\end{table}\n');
end


function write_eig_row(fid, rows, i, cfg)
%WRITE_EIG_ROW One eigenvalue row of the transposed table: index, then each run.
%
%   The notation is decided once for the whole row, so that a row is not part
%   fixed and part scientific.
    use_sci = row_needs_scientific(rows, i, cfg);
    fprintf(fid, '    $\\lambda_{%d}$', i);
    for j = 1:numel(rows)
        fprintf(fid, ' & %s', cell_str(rows(j), i, cfg, use_sci));
    end
    fprintf(fid, ' \\\\\n');
end


function tf = row_needs_scientific(rows, i, cfg)
%ROW_NEEDS_SCIENTIFIC Whether eigenvalue row i is to be written in scientific notation.
%
%   Fixed notation is the default, but a value near zero needs as many decimals
%   as it has leading zeros before it carries FMT_N significant digits, and one
%   such cell sets the width of its whole column. The Coffey--Evans ground state
%   is the case in point: it is zero, the methods return it as 1e-11 or so, and
%   six significant digits of that is seventeen decimals, which alone doubles
%   the natural width of the table. Such a row is written in scientific notation
%   instead, every cell of it, so that the row keeps one notation throughout.
    MAX_DECIMALS = 8;
    tf = false;
    if ~strcmp(cfg.fmt_mode, 'significant')
        return;
    end
    for j = 1:numel(rows)
        if i <= numel(rows(j).eigs)
            v = rows(j).eigs(i);
            if v ~= 0 && cfg.fmt_n - 1 - floor(log10(abs(v))) > MAX_DECIMALS
                tf = true;
                return;
            end
        end
    end
end


function s = cell_str(row, i, cfg, use_sci)
%CELL_STR One eigenvalue cell, bold for the reference row.
%
%   Digits come from the problem's FMT_MODE and FMT_N: a fixed number of
%   decimals, or a fixed number of significant digits, the decimals then varying
%   with the magnitude of the value. Plain fixed notation either way -- a cell
%   whose magnitude leaves no room for a decimal simply shows all its integer
%   digits, which only happens where a discretisation has broken down anyway.
%
%   A run holding fewer than i eigenvalues -- one whose N is below the index
%   asked for -- gets "---", the marker the DOF cell of the reference uses.
    if i > numel(row.eigs)
        s = '---';
        return;
    end
    % A scientific row carries two significant digits, not FMT_N. Six of them
    % plus a \times 10^{-11} typesets about as wide as the seventeen-decimal
    % fixed form it replaces, and one such cell sets the width of its column for
    % the whole table; two digits bring the cell back to the width of an
    % ordinary one. The row is the near-zero ground state, where the digits past
    % the second carry nothing anyway.
    SCI_DIGITS = 2;

    v = row.eigs(i);
    if use_sci
        if v == 0
            s = '0';
        else
            e = floor(log10(abs(v)));
            s = sprintf('%.*f \\times 10^{%d}', SCI_DIGITS - 1, v / 10^e, e);
        end
        if row.bold
            s = sprintf('\\mathbf{%s}', s);
        end
        s = sprintf('$%s$', s);
        return;
    end
    switch cfg.fmt_mode
        case 'decimals'
            ndec = cfg.fmt_n;
        case 'significant'
            if v == 0
                ndec = cfg.fmt_n - 1;
            else
                ndec = max(cfg.fmt_n - 1 - floor(log10(abs(v))), 0);
            end
        otherwise
            error('head_paper:fmt', 'unknown fmt_mode %s', cfg.fmt_mode);
    end
    s = sprintf('%.*f', ndec, v);
    if row.bold
        s = sprintf('\\textbf{%s}', s);
    end
end


function s = index_list(idx)
%INDEX_LIST Eigenvalue indices as "$\lambda_{50}$, $\lambda_{100}$ and ...".
    parts = arrayfun(@(k) sprintf('$\\lambda_{%d}$', k), idx, 'UniformOutput', false);
    if numel(parts) == 1
        s = parts{1};
    else
        s = sprintf('%s and %s', strjoin(parts(1:end-1), ', '), parts{end});
    end
end


function s = latex_name(name)
%LATEX_NAME Escape underscores for \texttt{} in text mode.
    s = strrep(name, '_', '\_');
end
