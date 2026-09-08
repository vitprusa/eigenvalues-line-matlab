function plot_dof_sweep_paper(name)
%PLOT_DOF_SWEEP_PAPER Paper DOF-sweep figure for a Sturm-Liouville test problem.
%
%   Line (one-dimensional) counterpart of the plane project's
%   experiments_paper/eigenvalues_dof_sweep_paper. It GENERATES ITS OWN DATA:
%   every method in src/ is run at DOF = 100, 300 and 500, and the resulting
%   spectra are plotted against the eigenvalue index, over the MATSLISE
%   reference of data/.
%
%   The methods are split over two panels, so that no panel carries more curves
%   than can be told apart by colour, and both repeat the reference:
%
%     (a)  DST, FD, CDM
%     (b)  MBCDM, Numerov, Numerov+AC, LGCC
%
%   Drawn with no title, colour coded: the method is encoded by colour and the
%   DOF level by line style; the MATSLISE reference is a thick solid black line.
%   Both axes are linear. A run with DOF below the plotted index range simply
%   stops where its spectrum ends.
%
%   Writes into results_paper/eigenvalues_dof_sweep/
%     - <problem>_dof_sweep_colour_a.eps and _b.eps, the two panels, and
%     - <problem>_dof_sweep_colour.tex, a \subfloat float holding both, and
%       (the _colour marks these as the colour figures, as in the plane project,
%       leaving the plain names free for a black-and-white variant)
%     - cache/<problem>_<method>_N<N>-eigenvalues.csv, the spectrum of every run
%       with its DOF count and measured time. A run whose CSV is present is read
%       back rather than recomputed; delete one to compute it again.
%
%   Called with no argument (or an empty one) it does both problems; pass a
%   problem name (paine, coffey_evans) to do just that one.
%
%   Requires dst/idst and Chebfun.
arguments
    % problem name, empty for every configured problem
    name (1, :) char = ''
end

    NEIG = 500;                     % eigenvalue indices drawn, 1 .. NEIG

    here         = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(here));
    addpath(fullfile(project_root, 'src'));

    out_dir   = fullfile(project_root, 'results_paper', 'eigenvalues_dof_sweep');
    cache_dir = fullfile(out_dir, 'cache');
    if ~exist(cache_dir, 'dir')
        mkdir(cache_dir);
    end

    if isempty(name)
        names = problem_names();
    else
        names = {name};
    end

    for ni = 1:numel(names)
        cfg = problem_config(names{ni}, project_root);
        fprintf('=== %s ===\n', cfg.name);

        reference = read_matslise_csv(cfg.reference_csv);
        reference = reference(1:min(NEIG, numel(reference)));

        results = compute_sweep(cfg, cache_dir);

        for g = 1:numel(cfg.groups)
            eps_file = fullfile(out_dir, sprintf('%s_dof_sweep_colour_%c.eps', ...
                                cfg.name, 'a' + g - 1));
            draw_panel(eps_file, cfg, results, cfg.groups{g}, reference, NEIG);
            fprintf('Wrote plot %s\n', eps_file);
        end

        tex = fullfile(out_dir, sprintf('%s_dof_sweep_colour.tex', cfg.name));
        write_latex(tex, cfg, NEIG);
        fprintf('Wrote %s\n', tex);
    end
end


function names = problem_names()
%PROBLEM_NAMES The configured problems, in the order they are done.
    names = {'paine', 'coffey_evans'};
end


function cfg = problem_config(name, project_root)
%PROBLEM_CONFIG Interval, potential, DOF schedule, panels and reference.
%
%   DOF is the size of the matrix solved, which is N for every method here, so
%   the sweep is over N, and every method walks the same three levels.
    DOFS = [100 300 500];

    methods = struct( ...
        'label', {'DST', 'FD', 'CDM', 'MBCDM', 'Numerov', 'Numerov+AC', 'LGCC'}, ...
        'key',   {'dst', 'fd', 'chebyshev', 'chebyshev_mapped', ...
                  'numerov', 'numerov_corrected', 'legendre_galerkin'}, ...
        'fun',   {@eig_dst, @eig_fd, @eig_chebyshev, @eig_chebyshev_mapped, ...
                  @eig_numerov, @eig_numerov_corrected, @eig_legendre_galerkin}, ...
        'N',     {DOFS, DOFS, DOFS, DOFS, DOFS, DOFS, DOFS});

    switch name
        case 'paine'
            cfg = struct( ...
                'name', 'paine', 'pretty', 'Paine', ...
                'q_text', '$q(x) = \mathrm{e}^{x}$', ...
                'a', 0, 'b', pi, 'q', @q_paine, 'methods', methods, ...
                'reference_csv', fullfile(project_root, 'data', ...
                                          'paine_matslise.csv'));
        case 'coffey_evans'
            cfg = struct( ...
                'name', 'coffey_evans', 'pretty', 'Coffey--Evans', ...
                'q_text', ['$q(x) = -2\beta\cos 2x + \beta^{2}\sin^{2} 2x$, ', ...
                           '$\beta = 30$'], ...
                'a', -pi/2, 'b', pi/2, 'q', @q_coffey_evans, 'methods', methods, ...
                'reference_csv', fullfile(project_root, 'data', ...
                                          'coffey_evans_matslise.csv'));
        otherwise
            error('dof_sweep:unknownProblem', ...
                  'Unknown problem "%s" (configured: %s).', name, ...
                  strjoin(problem_names(), ', '));
    end

    % Panels, by method key. Kept to at most four curve families each, that
    % being how many colours stay distinguishable in one set of axes.
    cfg.groups = {{'dst', 'fd', 'chebyshev'}, ...
                  {'chebyshev_mapped', 'numerov', 'numerov_corrected', ...
                   'legendre_galerkin'}};
end


function results = compute_sweep(cfg, cache_dir)
%COMPUTE_SWEEP Every method at every DOF level; cache the spectra; collect them.
    results = struct();
    for mi = 1:numel(cfg.methods)
        m = cfg.methods(mi);
        runs = {};
        for k = 1:numel(m.N)
            N   = m.N(k);
            csv = fullfile(cache_dir, sprintf('%s_%s_N%d-eigenvalues.csv', ...
                           cfg.name, m.key, N));
            if exist(csv, 'file')
                [evals, dofs, tsec] = read_head_csv(csv);
                fprintf('  %-11s N = %-5d  (cached)\n', m.label, N);
            else
                ensure_warm(m, cfg);
                t = tic;
                evals = m.fun(N, cfg.a, cfg.b, cfg.q);
                tsec = toc(t);
                evals = sort(real(evals(:)), 'ascend');
                dofs = numel(evals);
                write_head_csv(csv, cfg, m.label, N, dofs, tsec, evals);
                fprintf('  %-11s N = %-5d  dofs = %-6d  time = %8.2f s\n', ...
                        m.label, N, dofs, tsec);
            end
            runs{end+1} = struct('dofs', dofs, 'evals', evals, 'time', tsec); %#ok<AGROW>
        end
        results.(m.key) = runs;
    end
end


function ensure_warm(m, cfg)
%ENSURE_WARM Discarded runs of method M, before its first timing on a problem.
%
%   Covers MATLAB's JIT compilation and, for the Chebyshev and Legendre-Galerkin
%   methods, Chebfun's own load. Done at a small fixed size: the first scheduled
%   level of the Legendre-Galerkin method is not cheap, and a warm-up is by
%   definition thrown away.
    persistent warmed
    if isempty(warmed)
        warmed = struct();
    end
    key = sprintf('%s_%s', m.key, cfg.name);
    if isfield(warmed, key)
        return;
    end
    warmed.(key) = true;   % set first, so that a failed warm-up is not retried
    try
        for i = 1:2
            m.fun(40, cfg.a, cfg.b, cfg.q);
        end
    catch
        % Ignored: the real run reports its own failure.
    end
end


function draw_panel(eps_file, cfg, results, keys, reference, NEIG)
%DRAW_PANEL One panel: the reference and the sweep of the methods in KEYS.
%
%   Colour coded, as the plane project's colour variant: the method is encoded
%   by colour and the DOF level by line style, with a width bump on each cycle
%   of the styles so that a fifth level would stay distinct from the first. The
%   MATSLISE reference is a thick solid black line, drawn first so that it lies
%   under the method curves: several of them track it to ten significant digits
%   and a reference drawn on top would hide them completely.
%
%   Both axes are linear, as in the plane project. The ordinate is clipped just
%   above the reference: a run that breaks down at the top of its own spectrum
%   leaves the axes rather than setting their scale.
    colours = method_colours();
    styles  = {':', '-.', '--', '-'};

    % Render all text with the LaTeX interpreter, i.e. in the standard LaTeX
    % Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 1000 620], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');

    main = axes(fig);
    hs = draw_all(main, cfg, results, keys, reference, colours, styles, NEIG, true);
    xlabel(main, 'eigenvalue index $k$');
    ylabel(main, '$\lambda_k$');
    % No title: the problem and the panel are carried by the file name and the
    % subfloat caption.
    % The handles are passed explicitly: they are in the order the legend wants,
    % one column per method and the reference last, which is not the order the
    % curves are drawn in.
    legend(main, hs, 'Location', 'northwest', 'NumColumns', numel(keys) + 1, ...
           'FontSize', 7, 'Box', 'off');
    grid(main, 'on'); box(main, 'on');
    xlim(main, [0 NEIG]);
    % Clipped to the physical band: the runs that break down at the top of their
    % own spectrum run off the top of the axes.
    ylim(main, [0, 1.5 * reference(end)]);

    try
        exportgraphics(fig, eps_file, 'ContentType', 'vector');
    catch
        print(fig, eps_file, '-depsc2', '-painters');
    end
    close(fig);
end


function c = method_colours()
%METHOD_COLOURS Method key -> colour.
%
%   DST and FD keep the blue and orange the plane project gives them, and the
%   Chebyshev differentiation matrix the purple of the plane's Chebyshev
%   collocation; the rest follow MATLAB's default order. The panels carry three
%   and four methods, so no panel has to separate more than four colours.
    c = struct( ...
        'dst',               [0.0000 0.4470 0.7410], ...
        'fd',                [0.8500 0.3250 0.0980], ...
        'chebyshev',         [0.4940 0.1840 0.5560], ...
        'chebyshev_mapped',  [0.4660 0.6740 0.1880], ...
        'numerov',           [0.6350 0.0780 0.1840], ...
        'numerov_corrected', [0.3010 0.7450 0.9330], ...
        'legendre_galerkin', [0.9290 0.6940 0.1250]);
end


function hs = draw_all(ax, cfg, results, keys, reference, colours, styles, nmax, pad_legend)
%DRAW_ALL The curves of one panel onto one axes; returns the legend handles.
%
%   The reference goes down first, so that the method curves lie over it rather
%   than under it. The handles come back in legend order instead -- the methods,
%   each padded to a full column, and the reference last.
    labels = containers.Map({cfg.methods.key}, {cfg.methods.label});
    maxcount = max(cellfun(@(k) numel(results.(k)), keys));

    hold(ax, 'on');

    idx = 1:min(nmax, numel(reference));
    h_ref = plot(ax, idx, reference(idx), 'k-', 'LineWidth', 2.6, ...
                 'DisplayName', '$\mathtt{MATSLISE}$ (reference)');

    hs = gobjects(0);
    for mi = 1:numel(keys)
        key = keys{mi};
        runs = results.(key);
        for k = 1:numel(runs)
            r = runs{k};
            idx = 1:min(nmax, numel(r.evals));
            si = mod(k - 1, numel(styles)) + 1;               % cycle the 4 line styles
            lw = 2.0 + 1.0 * floor((k - 1) / numel(styles));  % thicker on each extra cycle
            hs(end+1) = plot(ax, idx, r.evals(idx), styles{si}, ...
                'Color', colours.(key), 'LineWidth', lw, ...
                'DisplayName', sprintf('$\\mathtt{%s}$ ($\\mathtt{DOF}$ = %d)', tt_label(labels(key)), r.dofs)); %#ok<AGROW>
        end
        % Pad this method's legend column to maxcount with invisible blank rows,
        % so the column-major legend keeps one column per method.
        if pad_legend
            for p = 1:(maxcount - numel(runs))
                hs(end+1) = plot(ax, NaN, NaN, 'LineStyle', 'none', ...
                                 'Marker', 'none', 'DisplayName', ' '); %#ok<AGROW>
            end
        end
    end
    hs(end+1) = h_ref;

    hold(ax, 'off');
end


function write_latex(tex, cfg, NEIG)
%WRITE_LATEX The \subfloat float holding both panels.
    fid = fopen(tex, 'w');
    if fid == -1
        error('dof_sweep:tex', 'Could not open %s for writing.', tex);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

    labels = {cfg.methods.label};
    keys   = {cfg.methods.key};
    group_labels = cell(1, numel(cfg.groups));
    for g = 1:numel(cfg.groups)
        tt = cellfun(@(k) sprintf('\\texttt{%s}', labels{strcmp(keys, k)}), ...
                     cfg.groups{g}, 'UniformOutput', false);
        group_labels{g} = sprintf('%s and %s', strjoin(tt(1:end-1), ', '), tt{end});
    end

    % Every method at once, for the caption of the whole float.
    tt_all = cellfun(@(m) sprintf('\\texttt{%s}', m), labels, 'UniformOutput', false);
    method_list = sprintf('%s and %s', strjoin(tt_all(1:end-1), ', '), tt_all{end});

    fprintf(fid, ['%% Eigenvalue DOF sweep for the %s problem \\texttt{%s}: the ', ...
        'eigenvalue index $k$\n%% against $\\lambda_k$, at several DOF levels ', ...
        'per method, over the\n%% \\texttt{MATSLISE} reference. Indices 1 to ', ...
        '%d.\n'], cfg.pretty, latex_name(cfg.name), NEIG);
    fprintf(fid, ['%% The methods are split over two panels, both repeating ', ...
        'the reference.\n%%\n%% Requires in the preamble:\n']);
    fprintf(fid, ['%%   \\usepackage{graphicx}\n', ...
                  '%%   \\usepackage{subfig}\n', ...
                  '%%   \\usepackage{epstopdf}   %% the figures are EPS; needed under pdflatex,\n', ...
                  '%%                           %% which then wants -shell-escape. Not needed for\n', ...
                  '%%                           %% latex/dvips, which reads EPS directly.\n%%\n']);
    fprintf(fid, ['%% The figures live next to this file, in\n', ...
        '%% results_paper/eigenvalues_dof_sweep/. If you \\input this snippet from\n', ...
        '%% elsewhere, point graphicx at that folder, e.g.\n', ...
        '%%   \\graphicspath{{results_paper/eigenvalues_dof_sweep/}}\n', ...
        '%% and drop the leading path from the \\includegraphics arguments below.\n\n']);

    fprintf(fid, '\\begin{figure}[htbp]\n  \\centering\n');
    for g = 1:numel(cfg.groups)
        if g > 1
            fprintf(fid, '  \\\\\n');
        end
        fprintf(fid, ['  \\subfloat[%s.', ...
                      '\\label{fig:dof_sweep_%s_colour_%c}]{%%\n', ...
                      '    \\includegraphics[width=0.8\\textwidth]{%s_dof_sweep_colour_%c}}\n'], ...
                group_labels{g}, cfg.name, 'a' + g - 1, cfg.name, 'a' + g - 1);
    end
    % The operator, the interval and the potential are left to the text that
    % introduces the problem, and the method names to the table of methods; the
    % caption opens with the problem name alone.
    fprintf(fid, ['  \\caption{%s problem. ', ...
        'Numerical eigenvalues computed by each discretisation %s with varying ', ...
        'degrees of freedom ($\\mathtt{DOF}$); compared with the reference ', ...
        'eigenvalues from \\texttt{MATSLISE}.}\n'], ...
        cfg.pretty, method_list);
    fprintf(fid, '  \\label{fig:dof_sweep_%s_colour}\n', cfg.name);
    fprintf(fid, '\\end{figure}\n');
end


function evals = read_matslise_csv(path)
%READ_MATSLISE_CSV Reference eigenvalues from a data/ CSV.
    fid = fopen(path, 'r');
    if fid == -1
        error('dof_sweep:reference', ...
              ['Could not open reference CSV %s. Generate it with ', ...
               'data/generate_matslise_reference.m.'], path);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    evals = [];
    col = NaN;
    while true
        line = fgetl(fid);
        if ~ischar(line); break; end
        s = strtrim(line);
        if isempty(s) || s(1) == '#'; continue; end
        parts = strtrim(strsplit(s, ','));
        if isnan(col)
            col = find(strcmp(parts, 'eigenvalue'), 1);
            if isempty(col)
                error('dof_sweep:reference', 'No "eigenvalue" column in %s.', path);
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
%READ_HEAD_CSV Read a cached spectrum back: eigenvalues, DOF count, time.
    fid = fopen(path, 'r');
    if fid == -1
        error('dof_sweep:csvread', 'Could not open %s.', path);
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
        error('dof_sweep:csv', 'Could not open %s for writing.', csv);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, ['# Problem: %s (%s, paper DOF-sweep run, ', ...
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


function s = latex_name(name)
%LATEX_NAME Escape underscores for \texttt{} in text mode.
    s = strrep(name, '_', '\_');
end


function s = tt_label(name)
%TT_LABEL Method name for \mathtt{} in the math mode of a figure legend.
%   The + of Numerov+AC is braced, so that it stays an ordinary symbol and
%   keeps the tight spacing of the name instead of being set as a binary
%   operator.
    s = strrep(name, '+', '{+}');
end
