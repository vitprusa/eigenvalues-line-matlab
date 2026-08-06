#!/usr/bin/env bash
# Recompute the paper eigenvalue-comparison spectra for the one-dimensional
# Sturm-Liouville test problems: every method in src/ at N = 500 and N = 1000,
# finite differences additionally at N = 10000, timing every run, and write one
# transposed LaTeX table per problem into results_paper/eigenvalues_head/.
# The spectra are cached as CSV in results_paper/eigenvalues_head/cache/ and are
# read back rather than recomputed, so regenerating the table alone costs no
# computation; delete a cached CSV to compute that run again.
# Generates its own data; the MATSLISE reference in data/ is read, never
# recomputed. Requires dst/idst and Chebfun.
#
# Pass a problem name (paine, coffey_evans) to regenerate just that one.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
name="${1:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); make_eigenvalues_head_paper_tables('${name}')"
