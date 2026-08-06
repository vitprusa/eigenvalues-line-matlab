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
| `eig_legendre_galerkin` | Legendre–Galerkin (LGCC) | basis `phi_k = P_k - P_{k+2}`; needs Chebfun, by far the slowest |

`eig_chebyshev_mapped` is the only one taking an extra argument:

```matlab
lambdas = eig_chebyshev_mapped(1000, 0, pi, @q_paine, 0.999);
```

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

## Validation

No automated test harness is committed. The following checks were run against
the solvers as they stand:

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
src_old/                the original scripts this code was refactored from, kept
                        unchanged for reference: one script per method and
                        problem pair (SLP_* for q = exp(x), CE_* for
                        Coffey-Evans) with the interval hardcoded and the CPU
                        timing and text-file writers inlined
```

The repository tracks the generated reference data, so `data/` is available
without rerunning anything; it can be regenerated with
`generate_matslise_reference.m`, which needs MATSLISE.

## Authors

The original scripts, kept in `src_old/`, were written by Oliver Křenek as part
of his master thesis, supervised by Vít Průša (<vit.prusa@matfyz.cuni.cz>). Vít
Průša is responsible for the conceptualisation of the work.

The refactored solvers in `src/`, the MATSLISE reference data and its generator
in `data/`, and the documentation were written by Claude Code (Claude Opus 5).

## License

The whole software is distributed under the BSD 3-Clause License.
