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
%     - two LaTeX snippets into results_paper/eigenvalues_head/, holding the
%       same booktabs table in the two layouts of WRITE_LATEX_TRANSPOSED:
%       <problem>_eigenvalues_head_transposed.tex splits the runs over two
%       subfloats and needs subfig, <problem>_..._transposed_one_table.tex puts
%       them all in one tabular and does not. Each is in the
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
    SPLIT = 3;                       % methods in the first of the two subtables

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
        write_latex_transposed(tex, cfg, rows, NEIG, EXTRA, SPLIT);
        fprintf('Wrote %s\n', tex);

        % The same data in one tabular, for a document that does not load subfig.
        tex = fullfile(out_dir, ...
              sprintf('%s_eigenvalues_head_transposed_one_table.tex', cfg.name));
        write_latex_transposed(tex, cfg, rows, NEIG, EXTRA, []);
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
end


function write_latex_transposed(tex, cfg, rows, NEIG, extras, split_after)
%WRITE_LATEX_TRANSPOSED Transposed booktabs table: one row per eigenvalue.
%
%   One column per run (the reference plus every method and size), one row per
%   eigenvalue index, so that a single eigenvalue can be read across all
%   discretisations. The method name spans its runs as a \multicolumn group head,
%   the DOF and timing become the two leading body rows, and the reference column
%   stays bold.
%
%   SPLIT_AFTER deals the runs over two \subfloat blocks of one float, the first
%   carrying that many methods and the second the rest, both repeating the
%   reference block so that either can be read on its own. Pass [] to put every
%   run in one tabular instead; that layout needs no subfig, at the cost of a
%   table wide enough that \resizebox has to scale it well down.
%
%   EXTRAS are eigenvalue indices past the leading NEIG, each written on a row of
%   its own after an empty row, so that the jump in the index is visible. A
%   non-empty EXTRAS also wraps each tabular in \resizebox, the deeper
%   eigenvalues being the widest cells and so setting every column width. Pass []
%   for the leading indices alone.
    fid = fopen(tex, 'w');
    if fid == -1
        error('head_paper:tex', 'Could not open %s for writing.', tex);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

    split = ~isempty(split_after);
    t = caption_bits(cfg);
    if split
        [parts, part_labels] = split_rows(rows, split_after);
    else
        parts = {rows};
    end

    fprintf(fid, ['%% Eigenvalue comparison for the %s problem \\texttt{%s}: ', ...
        '%s,\n'], cfg.pretty, latex_name(cfg.name), t.method_list_plain);
    if split
        fprintf(fid, ['%% %s, first %d eigenvalues, self-computed with timing.\n', ...
            '%% Transposed layout: one column per run, one row per eigenvalue,\n', ...
            '%% split over two subfloats that both repeat the reference block.\n'], ...
            t.lev_note, NEIG);
    else
        fprintf(fid, ['%% %s, first %d eigenvalues, self-computed with timing.\n', ...
            '%% Transposed layout: one column per run, one row per eigenvalue,\n', ...
            '%% every run in one tabular -- the same data as the subfloat version,\n', ...
            '%% for a document that does not load subfig.\n'], ...
            t.lev_note, NEIG);
    end
    if ~isempty(extras)
        fprintf(fid, ['%% Extended with the eigenvalues %s;\n', ...
            '%% a tabular is scaled down if it exceeds the text width.\n'], ...
            index_list(extras));
    end
    fprintf(fid, '%%\n%% Requires in the preamble:\n');
    fprintf(fid, ['%%   \\usepackage{booktabs}\n', ...
                  '%%   \\usepackage{amsmath}\n%%   \\usepackage{listings}\n']);
    if split
        fprintf(fid, '%%   \\usepackage{subfig}           %% \\subfloat\n');
    end
    if ~isempty(extras)
        fprintf(fid, '%%   \\usepackage{graphicx}          %% \\resizebox\n');
    end
    fprintf(fid, '\n');

    fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
    if split
        for k = 1:numel(parts)
            if k > 1
                % Stacks the second subfloat under the first rather than beside it.
                fprintf(fid, '  \\\\\n');
            end
            fprintf(fid, ['  \\subfloat[%s', ...
                          '\\label{tab:eigenvalues_head_%s_transposed_%c}]{%%\n'], ...
                    part_labels{k}, cfg.name, 'a' + k - 1);
            write_tabular(fid, cfg, parts{k}, NEIG, extras);
            fprintf(fid, '  }\n');
        end
        label_suffix = '';
    else
        write_tabular(fid, cfg, parts{1}, NEIG, extras);
        % Distinct label: the two layouts hold the same numbers and may well be
        % \input into the same document.
        label_suffix = '_one_table';
    end

    % The operator, the interval and the potential are left to the text that
    % introduces the problem; the caption opens with the problem name alone.
    fprintf(fid, ['  \\caption{%s problem. ', ...
        'Numerical eigenvalues computed by each discretisation %s with varying ', ...
        'degrees of freedom ($\\mathtt{DOF}$); compared with the reference ', ...
        'eigenvalues from \\texttt{MATSLISE}.}\n'], ...
        cfg.pretty, t.method_list);
    fprintf(fid, '  \\label{tab:eigenvalues_head_%s_transposed%s}\n', ...
            cfg.name, label_suffix);
    fprintf(fid, '\\end{table}\n');
end


function [parts, labels] = split_rows(rows, split_after)
%SPLIT_ROWS Rows for each subfloat, and the method list each one carries.
%
%   The reference block goes to both, so that a subfloat can be read without the
%   other; the method groups are dealt out in order, the first SPLIT_AFTER of
%   them to the first subfloat.
    is_truth = strcmp({rows.group}, 'truth');
    truth = rows(is_truth);
    rest  = rows(~is_truth);

    groups = unique({rest.group}, 'stable');
    sel = {groups(1:min(split_after, numel(groups))), ...
           groups(min(split_after, numel(groups)) + 1:end)};

    parts = cell(1, 2);
    labels = cell(1, 2);
    for k = 1:2
        parts{k} = [truth, rest(ismember({rest.group}, sel{k}))];
        tt = cellfun(@(m) sprintf('\\texttt{%s}', m), sel{k}, ...
                     'UniformOutput', false);
        if numel(tt) == 1
            labels{k} = sprintf('%s.', tt{1});
        else
            labels{k} = sprintf('%s and %s.', strjoin(tt(1:end-1), ', '), tt{end});
        end
    end
end


function write_tabular(fid, cfg, rows, NEIG, extras)
%WRITE_TABULAR One subfloat's tabular: the group head, DOF, timing, eigenvalues.
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

    % The size and \tabcolsep changes scope to the tabular (grouped) so that the
    % captions keep the normal body size.
    fprintf(fid, '    {\\scriptsize\n    \\setlength{\\tabcolsep}{3pt}\n');
    if ~isempty(extras)
        % Scaled DOWN only: \width is the natural width of the tabular, so a
        % table that already fits is passed through at its own size rather than
        % blown up to fill the text width. Trailing %% so that the line break
        % adds no space before the tabular.
        fprintf(fid, ['    \\resizebox{\\ifdim\\width>\\textwidth\\textwidth', ...
                      '\\else\\width\\fi}{!}{%%\n']);
    end
    fprintf(fid, '    \\begin{tabular}{%s}\n      \\toprule\n', colspec);

    % Group head: method name spanning its runs, underlined by \cmidrule.
    fprintf(fid, '     ');
    for g = 1:numel(starts)
        lab = rows(starts(g)).method_label;
        if ~strcmp(rows(starts(g)).group, 'truth')
            lab = sprintf('\\texttt{%s}', lab);
        end
        fprintf(fid, ' & \\multicolumn{%d}{c}{%s}', stops(g) - starts(g) + 1, lab);
    end
    fprintf(fid, ' \\\\\n     ');
    for g = 1:numel(starts)
        fprintf(fid, ' \\cmidrule(lr){%d-%d}', starts(g) + 1, stops(g) + 1);
    end
    fprintf(fid, '\n');

    % Run metadata, then the eigenvalues.
    fprintf(fid, '      $\\mathtt{DOF}$');
    fprintf(fid, ' & %s', rows.dof);
    fprintf(fid, ' \\\\\n      Time (s)');
    fprintf(fid, ' & %s', rows.time);
    fprintf(fid, ' \\\\\n      \\midrule\n');

    for i = 1:NEIG
        write_eig_row(fid, rows, i, cfg);
    end
    % The deeper eigenvalues, each after an empty row: the index jumps from one
    % to the next, and a rule would read as a new block of the same sequence.
    for i = extras(:)'
        fprintf(fid, '     %s \\\\\n', repmat(' &', 1, ncol));
        write_eig_row(fid, rows, i, cfg);
    end

    if isempty(extras)
        fprintf(fid, '      \\bottomrule\n    \\end{tabular}}\n');
    else
        fprintf(fid, '      \\bottomrule\n    \\end{tabular}}}\n');   % tabular, \resizebox, size group
    end
end


function write_eig_row(fid, rows, i, cfg)
%WRITE_EIG_ROW One eigenvalue row of the transposed table: index, then each run.
    fprintf(fid, '      $\\lambda_{%d}$', i);
    for j = 1:numel(rows)
        fprintf(fid, ' & %s', cell_str(rows(j), i, cfg));
    end
    fprintf(fid, ' \\\\\n');
end


function s = cell_str(row, i, cfg)
%CELL_STR One eigenvalue cell, bold for the reference row.
%
%   In 'decimals' mode a cell carries FMT_N digits after the point.
%
%   In 'significant' mode it carries FMT_N significant digits, the decimals
%   varying with the magnitude, while fixed notation can show exactly that many
%   -- which it can while the value has at most FMT_N digits before the point
%   and no long run of leading zeros after it. Outside that range fixed notation
%   turns ugly in one of two ways, and the cell is written in scientific
%   notation to SCI_DIGITS significant digits instead:
%
%     too large  a value of 1.2e9 has no room for a decimal and prints every one
%               of its ten integer digits, four more than were asked for
%     too small  a value of 1e-11 needs sixteen decimals before its sixth
%               significant digit appears
%
%   Either way the cell would be far wider than an ordinary one, and a column is
%   as wide as its widest cell. The large case is where a discretisation has
%   broken down at the top of its own spectrum, the small case the near-zero
%   Coffey--Evans ground state; in both the digits past the second carry
%   nothing, hence SCI_DIGITS = 2.
%
%   A run holding fewer than i eigenvalues -- one whose N is below the index
%   asked for -- gets "---", the marker the DOF cell of the reference uses.
    MAX_DECIMALS = 8;
    SCI_DIGITS = 2;

    if i > numel(row.eigs)
        s = '---';
        return;
    end
    v = row.eigs(i);

    if strcmp(cfg.fmt_mode, 'decimals')
        s = bold_fixed(sprintf('%.*f', cfg.fmt_n, v), row.bold);
        return;
    end
    if ~strcmp(cfg.fmt_mode, 'significant')
        error('head_paper:fmt', 'unknown fmt_mode %s', cfg.fmt_mode);
    end

    if v == 0
        s = bold_fixed(sprintf('%.*f', cfg.fmt_n - 1, v), row.bold);
        return;
    end

    e = floor(log10(abs(v)));
    ndec = cfg.fmt_n - 1 - e;
    if ndec >= 0 && ndec <= MAX_DECIMALS
        s = bold_fixed(sprintf('%.*f', ndec, v), row.bold);
        return;
    end

    s = sprintf('%.*f \\times 10^{%d}', SCI_DIGITS - 1, v / 10^e, e);
    if row.bold
        s = sprintf('\\mathbf{%s}', s);
    end
    s = sprintf('$%s$', s);
end


function s = bold_fixed(s, bold)
%BOLD_FIXED A fixed-notation cell, bold for the reference row.
    if bold
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
