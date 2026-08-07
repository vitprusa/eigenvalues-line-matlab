# Discrete sine transform (DST) based computation of Sturm–Liouville eigenvalues

Discrete sine transform (DST) based method for the eigenvalues of the
one-dimensional Sturm–Liouville operator on an interval, and a comparison
against the standard alternatives.

## Overview

The core idea is to discretise the second-derivative operator with homogeneous
Dirichlet boundary conditions using the **discrete sine transform (DST)**, which
diagonalises that operator on an interval. See

> Fusi, Lorenzo, Oliver Křenek, Vít Průša, Casey Rodriguez, Rebecca Tozzi, and Martin Vejvoda. "Discrete versus continuous—Linear lattice models and their exact continuous counterparts." International Journal of Engineering Science 224 (2026): 104530, [10.1016/j.ijengsci.2026.104530](https://doi.org/10.1016/j.ijengsci.2026.104530)

for reference and thorough discussion.

The problem solved here is the regular Sturm–Liouville problem in Liouville
normal form,

```
-y'' + q(x) y = lambda y,    y(a) = y(b) = 0,
```

that is `p = w = 1`, with Dirichlet conditions at both endpoints. Differentiation
is diagonal in DST space, so on the uniform interior grid of `N` points the
operator is `idst(kvec .* dst(y))` with `kvec = (k*pi/(b-a))^2`, and the
multiplication by `q` is diagonal on the same grid. The matrix is formed by
applying the operator to the identity and handed to `eig`.

The repository also contains finite-difference, Numerov, Chebyshev and
Legendre–Galerkin implementations under `src/`, used to cross-check the DST
results, and reference eigenvalues computed with MATSLISE under `data/`.

This code is a refactoring of the MATLAB codes in the GitHub repository
<https://github.com/krenekoliver/Masters-thesis>, and builds on Oliver Křenek's
master thesis "Rate-type models for metamaterials" (Charles University, Faculty
of Mathematics and Physics, 2026),
[hdl.handle.net/20.500.11956/209791](http://hdl.handle.net/20.500.11956/209791).

## Requirements

- MATLAB.
  - `dst`/`idst` — used by `eig_dst`.
  - **Chebfun** — only for `eig_chebyshev`, `eig_chebyshev_mapped` and
    `eig_legendre_galerkin`, which use `chebpts`, `diffmat`, `chebfun`,
    `legpoly` and `legcoeffs`; expected on the MATLAB path.
- **MATSLISE** — only to regenerate the reference data in `data/`, not to run
  the solvers. Expected as a directory `MATSLISE2` sibling to the repository
  root; see `data/README.md`.

## Getting started

Put the source folder on the path and call any solver:

```matlab
addpath('src')

lambdas = eig_dst(1000, 0, pi, @q_paine);
lambdas = eig_fd(10000, -pi/2, pi/2, @q_coffey_evans);
```

Every solver has the same signature,

```matlab
lambdas = eig_<method>(N, a, b, q)
```

taking the matrix size `N`, the interval endpoints `a` and `b`, and the
potential `q` as a function handle, and returning `N` eigenvalues sorted in
ascending order. Eigenvalues are indexed from one: `lambdas(1)` is the ground
state.

The methods carry no timing and write no files; they return the eigenvalues and
nothing else.

## Methods (`src/`)

| Function | Method | Notes |
|----------|--------|-------|
| `eig_dst` | discrete sine transform | the method of interest; needs `dst`/`idst` |
| `eig_fd` | second-order finite differences | three-point stencil |
| `eig_chebyshev` | Chebyshev differentiation matrix (CDM) | `N + 2` Chebyshev points, boundary rows removed; needs Chebfun |
| `eig_chebyshev_mapped` | mapped barycentric Chebyshev (MBCDM) | Kosloff–Tal-Ezer map, optional fifth argument `alpha` in `(0, 1]`, default `0.999`; needs Chebfun |
| `eig_numerov` | Numerov | generalised eigenvalue problem `A y = lambda B y` |
| `eig_numerov_corrected` | Numerov with the asymptotic correction | Paine–de Hoog–Anderssen correction on top of `eig_numerov` |
| `eig_legendre_galerkin` | Legendre–Galerkin (LGCC) | basis `phi_k = P_k - P_{k+2}`; needs Chebfun |

`eig_chebyshev_mapped` is the only one taking an extra argument:

```matlab
lambdas = eig_chebyshev_mapped(1000, 0, pi, @q_paine, 0.999);
```

The same seven methods are described for the article in
`results_paper/latex_helpers/methods_description_table.tex`, see
[LaTeX helpers](#latex-helpers-results_paperlatex_helpers) below.

## Test problems (`src/q_*.m`)

The potentials are function files, passed to the solvers as handles:

| Function | `q(x)` | Interval |
|----------|--------|----------|
| `q_paine` | `exp(x)` | `[0, pi]` |
| `q_coffey_evans` | `-2*beta*cos(2x) + beta^2*sin(2x)^2`, `beta = 30` by default | `[-pi/2, pi/2]` |

`q_coffey_evans` takes `beta` as an optional second argument, so a different
value is used by wrapping it:

```matlab
lambdas = eig_dst(1001, -pi/2, pi/2, @(x) q_coffey_evans(x, 50));
```

## Reference data (`data/`)

`data/` holds 500 reference eigenvalues per problem, computed with MATSLISE from
the command line, as `paine_matslise.csv` and `coffey_evans_matslise.csv`. Each
CSV opens with a `#`-commented header recording the problem, the method, the
citations, the exact calls, the wall clock time under a `WALL_CLOCK` key, and
the column meanings; there is no timestamp, so re-running leaves the files
unchanged unless the numbers do. Read them with:

```matlab
T = readtable('data/paine_matslise.csv', 'CommentStyle', '#');
```

Both index conventions are carried in the data, since MATSLISE numbers the
ground state 0 while this project numbers it 1: the `k` column counts from one
and `matslise_index` counts from zero. Regenerate with

```matlab
generate_matslise_reference(500, 1e-12)
```

run from `data/`. See `data/README.md` for the method, the citations, and where
MATSLISE has to be located.

## Paper experiments (`experiments_paper/`)

Two experiments, each with a driver that generates its own data and a shell
runner beside it, writing into a directory of the same name under
`results_paper/`.

### Eigenvalue comparison (`eigenvalues_head`)

`experiments_paper/eigenvalues_head/` produces the article's eigenvalue
comparison. It generates its own data: for each test problem it runs every
method in `src/` at `N = 500` and `N = 1000`, finite differences additionally at
`N = 10000`, and times each run from a warm start, a throwaway call absorbing
the JIT compilation and Chebfun's load. From a shell:

```bash
experiments_paper/eigenvalues_head/run_eigenvalues_head_paper.sh
experiments_paper/eigenvalues_head/run_eigenvalues_head_paper.sh paine
```

The optional argument restricts the run to one problem (`paine`,
`coffey_evans`).

It writes one transposed table per problem into
`results_paper/eigenvalues_head/`: rows are the eigenvalue indices and columns
the runs — the MATSLISE reference plus every method and size — so that a single
eigenvalue reads across all discretisations, with DOF and timing as the two
leading body rows. The table carries `lambda_1` to `lambda_8` and then
`lambda_50`, `lambda_100`, `lambda_400` and `lambda_500`.

The spectra are cached as CSV under `results_paper/eigenvalues_head/cache/`,
each carrying its DOF count and measured time in the header, and are read back
rather than recomputed: the first run takes about two and a half minutes, a
redraw of the table four seconds. Delete a cached CSV to compute that run again.

The number of digits is set per problem in `problem_configs`, as `fmt_mode` and
`fmt_n`: `'decimals'` prints `fmt_n` digits after the point, `'significant'`
prints `fmt_n` significant digits by varying the decimals with the magnitude.
Both problems are set to six significant digits.

Fixed notation is used while it can show exactly that many digits, which it can
while the value has at most six digits before the point and no long run of
leading zeros after it. A cell outside that range is written in scientific
notation to two significant digits instead. There are two ways out of the range:

- **too large** — `1.2e9` has no room for a decimal, so fixed notation prints
  all ten of its integer digits, four more than were asked for;
- **too small** — `1e-11` needs sixteen decimals before its sixth significant
  digit appears.

Either way the cell would be far wider than an ordinary one, and a column is as
wide as its widest cell. The test is applied per cell, not per row, so a single
run falling out of range does not change the notation of the rest of its row.

### DOF sweep (`eigenvalues_dof_sweep`)

`experiments_paper/eigenvalues_dof_sweep/` produces the article's DOF-sweep
figures. It also generates its own data: for each test problem it runs every
method in `src/` at `DOF = 100, 300` and `500` — DOF being the size of the
matrix solved, which is `N` for every method here — and draws the resulting
spectra against the eigenvalue index, indices 1 to 500, over the MATSLISE
reference of `data/`. From a shell:

```bash
experiments_paper/eigenvalues_dof_sweep/run_dof_sweep_paper.sh
experiments_paper/eigenvalues_dof_sweep/run_dof_sweep_paper.sh paine
```

With no argument it does both problems; the optional argument restricts the run
to one (`paine`, `coffey_evans`).

For each problem it writes two figures and a snippet holding both into
`results_paper/eigenvalues_dof_sweep/`, the methods split so that no figure has
to separate more than four colours:

| figure | curves |
| --- | --- |
| `<problem>_dof_sweep_colour_a.eps` | reference, `DST`, `FD`, `CDM` |
| `<problem>_dof_sweep_colour_b.eps` | reference, `MBCDM`, `Numerov`, `Numerov+AC`, `LGCC` |

Both repeat the reference, so either can be read on its own. The figures carry
no title, the problem and the panel being named by the file name and the
caption. The method is encoded by colour and the DOF level by line style; the
reference is a thick solid black line, drawn before the method curves so that it
lies under them rather than over them, several of the methods tracking it too
closely to remain visible otherwise. Both axes are linear, with the ordinate
clipped just above the reference, so a run that breaks down at the top of its
own spectrum leaves the axes rather than setting their scale. The `_colour`
suffix marks these as the colour figures, leaving the plain names free for a
black-and-white variant.

The spectra are cached as CSV under `results_paper/eigenvalues_dof_sweep/cache/`
in the same format as the comparison experiment, and are likewise read back
rather than recomputed: a problem takes about half a minute from cold and a
redraw of the figures costs nothing. Adding a problem is a case of
`problem_config` and an entry of `problem_names`.

## LaTeX helpers (`results_paper/latex_helpers`)

Snippets written by hand rather than by an experiment, kept next to the
generated ones so the article can `\input` them the same way:

| file | content |
| --- | --- |
| `methods_description_table.tex` | a `booktabs` table describing all seven discretisations of `src/`, one row per method |

The table carries the labels the figures and the comparison tables use (`DST`,
`FD`, `CDM`, `MBCDM`, `Numerov`, `Numerov+AC`, `LGCC`), and cites Chebfun, the
Paine–de Hoog–Anderssen correction and the compact Legendre basis. Its label is
`tab:methods_description_line`, distinct from the `tab:methods_description` of
the companion table of the plane project, since both snippets end up in the same
document. Nothing regenerates it: keep it in step with `src/eig_*.m` by hand.

## Validation

No automated test harness is committed. The paper experiments above are the
standing cross-method check: they recompute every method against the MATSLISE
reference at each resolution, and their cached CSVs hold the full spectrum, the
DOF count and the measured time of every run. The numbers are left to be read
off the tables, the figures and the cache.

Three further checks were run against the solvers as they stand:

- All seven methods agree with the MATSLISE reference to the digits printed on
  both test problems: Paine `4.896669  10.045190  16.019267  23.266271
  32.263707`, Coffey–Evans `0  117.946  231.665  231.665  231.665`.
- The Coffey–Evans spectrum was confirmed independently by a tridiagonal
  finite-difference solve at `N = 80000` outside MATLAB, which reproduces the
  level at `117.9463` between the ground state and the well-known triple
  cluster at `231.6649`.
- `eig_numerov_corrected` reproduces the `q = 0` spectrum to a relative
  `3e-12` on `[0, pi]`, `[-pi/2, pi/2]` and `[0, 2]`, against `3.9e-1` for
  the uncorrected `eig_numerov` — the correction is exact for `q = 0` by
  construction, so this also pins down the index alignment.

The generated tables were checked with `pdflatex` against an `amsart` preamble
matching the sibling project's (`a4paper`, `geometry scale=0.9`): both compile
with no overfull or underfull boxes.

## Repository layout

```
src/                    one function per method, plus the two potentials
src/eig_dst.m           discrete sine transform
src/eig_fd.m            second-order finite differences
src/eig_chebyshev.m     Chebyshev differentiation matrix
src/eig_chebyshev_mapped.m   mapped barycentric Chebyshev
src/eig_numerov.m       Numerov
src/eig_numerov_corrected.m  Numerov with the asymptotic correction
src/eig_legendre_galerkin.m  Legendre-Galerkin
src/q_paine.m           potential exp(x)
src/q_coffey_evans.m    Coffey-Evans potential, beta = 30 by default
data/                   MATSLISE reference eigenvalues (500 per problem) as
                        <problem>_matslise.csv, the generator
                        generate_matslise_reference.m that produced them, and a
                        README covering the method, the citations and where
                        MATSLISE must be located
experiments_paper/eigenvalues_head/ the article's eigenvalue comparison:
                        make_eigenvalues_head_paper_tables.m generates the data
                        and writes the tables, run_eigenvalues_head_paper.sh
                        runs it headless
experiments_paper/eigenvalues_dof_sweep/ the article's DOF sweep:
                        plot_dof_sweep_paper.m generates the data and draws the
                        figures, run_dof_sweep_paper.sh runs it headless
results_paper/eigenvalues_head/ one transposed table per problem,
                        <problem>_eigenvalues_head_transposed.tex, plus cache/
                        holding the spectrum of every run as
                        <problem>_<method>_N<N>-eigenvalues.csv with its DOF
                        count and measured time; the CSVs are the cache the
                        table is drawn from, delete one to recompute that run
results_paper/eigenvalues_dof_sweep/ two figures per problem,
                        <problem>_dof_sweep_colour_a.eps and _b.eps, and
                        <problem>_dof_sweep_colour.tex holding both, plus
                        cache/ in the same format as above
results_paper/latex_helpers/ hand-written LaTeX snippets, currently
                        methods_description_table.tex describing the seven
                        methods of src/; not regenerated by any experiment
src_old/                the original scripts this code was refactored from, kept
                        unchanged for reference: one script per method and
                        problem pair (SLP_* for q = exp(x), CE_* for
                        Coffey-Evans) with the interval hardcoded and the CPU
                        timing and text-file writers inlined
```

The repository tracks the generated artifacts, so the results are available
without rerunning anything: the MATSLISE reference under `data/`, and the
article's tables and the spectra behind them under `results_paper/`. They can be
regenerated with `data/generate_matslise_reference.m`, which needs MATSLISE, and
`experiments_paper/eigenvalues_head/run_eigenvalues_head_paper.sh`, which needs
only MATLAB and Chebfun.

## Authors

The original scripts, kept in `src_old/`, were written by Oliver Křenek as part
of his master thesis, supervised by Vít Průša (<vit.prusa@matfyz.cuni.cz>). Vít
Průša is responsible for the conceptualisation of the work.

The refactored solvers in `src/`, the MATSLISE reference data and its generator
in `data/`, the paper experiment in `experiments_paper/`, and the documentation
were written by Claude Code (Claude Opus 5).

## License

The whole software is distributed under the BSD 3-Clause License.
