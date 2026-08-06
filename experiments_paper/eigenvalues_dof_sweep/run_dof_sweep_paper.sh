#!/usr/bin/env bash
# Generate the paper DOF-sweep figures (no titles, colour coded) for the
# one-dimensional Sturm-Liouville test problems: every method in src/ at
# DOF = 100, 300 and 500, plotted against the eigenvalue index over the
# MATSLISE reference of data/.
#
# Writes the two panels as EPS and one \subfloat snippet holding both into
# results_paper/eigenvalues_dof_sweep/. The spectra are cached as CSV in
# results_paper/eigenvalues_dof_sweep/cache/ and are read back rather than
# recomputed; delete a cached CSV to compute that run again.
#
# Does both problems; pass a problem name (paine, coffey_evans) to restrict the
# run to one.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
name="${1:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); plot_dof_sweep_paper('${name}')"
