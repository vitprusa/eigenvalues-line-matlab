# Reference eigenvalues

High-accuracy eigenvalues for the two test problems used in this project,
computed with MATSLISE and used as the reference against which the methods in
`../src` are compared.

## What is computed

Both problems are regular Sturm–Liouville problems in Liouville normal form,

```
-y'' + q(x) y = lambda y,    y(a) = y(b) = 0,
```

that is `p = w = 1`, with Dirichlet conditions at both endpoints. For each
problem the first 500 eigenvalues are computed.

| file | problem | q(x) | interval |
| --- | --- | --- | --- |
| `paine_matslise.csv` | Paine problem 1 | `exp(x)` | `[0, pi]` |
| `coffey_evans_matslise.csv` | Coffey–Evans, `beta = 30` | `-2*beta*cos(2x) + beta^2*sin(2x)^2` | `[-pi/2, pi/2]` |

The potentials are not re-typed here: `generate_matslise_reference.m` passes the
function handles `@q_paine` and `@q_coffey_evans` from `../src` straight to
MATSLISE, so the reference and the methods under test always use the same `q`.
`q_coffey_evans` is called without a second argument, so `beta` takes its
default value of 30.

## Library: MATSLISE

The reference values come from **MATSLISE 2** (Ledoux, Van Daele, Vanden Berghe,
Department of Applied Mathematics, Computer Science and Statistics, Ghent
University), a MATLAB package for Sturm–Liouville problems. It uses CP methods
— a piecewise constant reference potential plus perturbation corrections — on an
automatically constructed mesh, and returns an error estimate and a status flag
for every eigenvalue.

### Reference

Cite the toolbox as (Chicago style):

> Ledoux, Veerle, and Marnix Van Daele. "Matslise 2.0: A Matlab Toolbox for
> Sturm–Liouville Computations." *ACM Transactions on Mathematical Software* 42,
> no. 4 (2016): 29:1–29:18. https://doi.org/10.1145/2839299.

The original package, which version 2 revises and which the bundled help files
cite, is:

> Ledoux, Veerle, Marnix Van Daele, and Guido Vanden Berghe. "MATSLISE: A MATLAB
> Package for the Numerical Solution of Sturm–Liouville and Schrödinger
> Equations." *ACM Transactions on Mathematical Software* 31, no. 4 (2005):
> 532–554. https://doi.org/10.1145/1114268.1114273.

The software itself is distributed on SourceForge, as `MATSLISE2.zip`:

> Ledoux, Veerle, and Marnix Van Daele. *MATSLISE*. Version 2. MATLAB package.
> SourceForge, 2015. https://sourceforge.net/projects/matslise/.

The project home page is <https://matslise.sourceforge.io/>; the download page is
<https://sourceforge.net/projects/matslise/files/>. MATSLISE is free for
non-commercial use. All three references are repeated in the header of every
generated CSV.

### Command-line use

MATSLISE ships with a GUI, but it is used here purely from the command line:

```matlab
one = @(x) ones(size(x));
o   = slp(one, @q_coffey_evans, one, -pi/2, pi/2, 1, 0, 1, 0);
E   = computeEigenvalues(o, 0, 499, 1e-12, true);
```

The requested input tolerance is `1e-12`. The exact calls are also recorded in
the header of each CSV.

### Where MATSLISE must be located

`generate_matslise_reference.m` resolves the library **relative to its own
location**, as a directory named `MATSLISE2` that is a sibling of
`eigenvalues-line-matlab`:

```
<parent>/
├── MATSLISE2/                     <- unpack MATSLISE here
│   └── source/                       (added to the path with genpath)
└── eigenvalues-line-matlab/
    ├── src/                       <- q_paine.m, q_coffey_evans.m, solvers
    └── data/                      <- this directory
        ├── README.md
        ├── generate_matslise_reference.m
        ├── paine_matslise.csv
        └── coffey_evans_matslise.csv
```

So from this directory MATSLISE is at `../../MATSLISE2`, and only its `source`
subdirectory is put on the MATLAB path. The script adds it itself; there is no
need to run MATSLISE's own `startup.m` first, and the `examples` directory is
deliberately left off the path. If MATSLISE lives anywhere else, adjust the
`matslise_dir` line near the top of `generate_matslise_reference.m`.

## Regenerating

From this directory:

```
matlab -batch "generate_matslise_reference(500, 1e-12)"
```

Both arguments are optional and default to `500` and `1e-12`. The script
overwrites both CSV files. A throwaway MATSLISE call is made before the timed
ones so that the recorded times exclude JIT compilation and class loading.

## Output format

Each CSV opens with a `#`-commented header giving the problem, the method, the
references, the exact commands used, the measured time, and the environment,
followed by a normal header row and 500 data rows. The header carries no
timestamp, so re-running the script leaves the files byte-identical unless the
numbers or the timings actually change.

The timing is written under greppable keys, so it can be recovered later without
reading the whole header:

```
$ grep WALL_CLOCK *.csv
coffey_evans_matslise.csv:#   WALL_CLOCK = 0.227644 s
coffey_evans_matslise.csv:#   WALL_CLOCK_PER_EIGENVALUE = 0.000455288 s
paine_matslise.csv:#   WALL_CLOCK = 0.154234 s
paine_matslise.csv:#   WALL_CLOCK_PER_EIGENVALUE = 0.000308468 s
```

`WALL_CLOCK` covers the single `computeEigenvalues` call that returns all 500
eigenvalues, mesh construction included, measured from a warm start.

Read the data with:

```matlab
T = readtable('paine_matslise.csv', 'CommentStyle', '#');
```

Columns:

| column | meaning |
| --- | --- |
| `k` | eigenvalue index counting **from one**, the convention of the article |
| `matslise_index` | index counting **from zero**, as returned by MATSLISE |
| `eigenvalue` | `lambda_k` |
| `estimated_error` | MATSLISE error estimate for that eigenvalue |
| `status` | MATSLISE status flag; `0` means correct calculation |

Both index columns are present because the two conventions differ by one and
the mismatch is easy to introduce. MATSLISE numbers the ground state 0; the
article numbers it 1. Requesting indices `0 … 499` therefore yields
`k = 1 … 500`.

## Results as generated

MATLAB R2026a, GLNXA64. All 500 status flags were `0` in both runs.

| problem | `WALL_CLOCK` | s / eigenvalue | max abs. estimated error |
| --- | --- | --- | --- |
| Paine 1 | 0.154 s | 3.1e-4 | 9.8e-9 |
| Coffey–Evans | 0.228 s | 4.6e-4 | 4.6e-8 |

The eigenvalues are reproducible to the last digit printed, but the timings are
not: the Coffey–Evans run has been observed between 0.23 s and 0.47 s on the same
machine. Treat `WALL_CLOCK` as an order-of-magnitude figure, and re-time on a
quiet machine if the number is going into the article.

Spot values: Paine `k=1` → 4.89666937996706, `k=500` → 250007.047629702;
Coffey–Evans `k=1` → -4.24e-11, `k=2` → 117.946307662133, `k=3,4,5` →
231.66492923… (the near-triple cluster), `k=500` → 250450.103049251.
