# FORE-reset-benchmark

Reproducible Python implementation of the FORE zero-crossing reset benchmark, with weighted cocycle and EDMD estimates for spectral stability analysis.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20682315.svg)](https://doi.org/10.5281/zenodo.20682315)

## Overview

This repository contains the Python code used to reproduce the numerical benchmarks associated with the manuscript:

**A Weighted Spectral Criterion for Uniform Exponential Stability of Homogeneous Zero-Crossing Reset Maps**

The repository includes two reproducible benchmarks:

1. **FORE zero-crossing reset benchmark**  
   This script computes the full state-dependent post-reset map of the FORE benchmark, evaluates finite-horizon weighted cocycle estimates, computes projected-data EDMD spectral diagnostics, and generates the CSV table and figure used in the manuscript.

2. **Synthetic \(S^2\) weighted-reset benchmark**  
   This script defines a controlled three-dimensional homogeneous induced map on \(\mathbb R^3\), with angular dynamics on \(S^2\), to illustrate a parameter-dependent spectral stability transition. The example shows how finite-horizon cocycle and weighted EDMD diagnostics detect a crossing of the stability threshold \(1\), while the one-step gain \(\max 
ho_\mu\) remains above \(1\).

Both benchmarks are designed to illustrate the weighted spectral-radius mechanism for homogeneous reset-induced maps. The numerical estimates are diagnostic approximations and are not claimed to be rigorous interval-certified upper bounds.

## Repository contents

The repository should contain at least the following files:

- `run_fore_benchmark.py`  
  Reproduces the FORE zero-crossing reset benchmark.

- `run_synthetic_s2_benchmark.py`  
  Reproduces the synthetic \(S^2\) weighted-reset benchmark.

- `README.md`  
  This file.

- `LICENSE`  
  MIT License file.

Recommended optional files before public deposit:

- `CITATION.cff` or `codemeta.json`  
  Machine-readable citation metadata for the software repository.

## Generated output files

### FORE benchmark

Running

```bash
python run_fore_benchmark.py
```

generates:

- `fore_benchmark_results.csv`  
  CSV table containing, for each dwell-time value \(	au_m\), the terminal finite-horizon cocycle estimate, the weighted EDMD estimate, the maximum one-step radial gain, the frozen-branch reference value, return-time ranges, and radial-gain ranges.

- `nahs_fore_benchmark.pdf`  
  Vector figure used in the manuscript.

- `nahs_fore_benchmark.png`  
  Raster version of the same figure.

### Synthetic \(S^2\) benchmark

Running

```bash
python run_synthetic_s2_benchmark.py
```

generates:

- `synthetic_s2_benchmark_results.csv`  
  CSV table containing, for each scaling value \(\mu\), the finite-horizon cocycle estimate, the weighted EDMD estimate, the maximum one-step radial gain, and the stability interpretation.

- `synthetic_s2_weighted_benchmark.pdf`  
  Vector figure used in the manuscript.

## Software requirements

The scripts were written for Python 3 and require the following Python packages:

- `numpy`
- `scipy`
- `matplotlib`

The scripts only use standard-library modules in addition to the packages listed above, such as `csv`, `dataclasses`, `pathlib`, and `typing`.

A typical installation command is:

```bash
python -m pip install numpy scipy matplotlib
```

## How to run

From the directory containing the scripts, run:

```bash
python run_fore_benchmark.py
python run_synthetic_s2_benchmark.py
```

### Portability note

The distributed scripts may have their default output directory set for the original sandbox environment, for example:

```python
out_dir: Path = Path("/mnt/data")
```

or

```python
OUTDIR = Path("/mnt/data")
```

For normal local, GitHub, or Zenodo use, change these lines to:

```python
out_dir: Path = Path(".")
```

or

```python
OUTDIR = Path(".")
```

or to any desired local output directory before running the scripts.

## Reproducibility information

### FORE benchmark

The FORE script uses deterministic uniform grids and does not use random sampling. Therefore, no random seed is required.

Default numerical parameters are defined in the `Config` dataclass:

- `taus = (0.10, 0.20, 0.40, 0.60, 0.80)`
- `n_tab = 200`
- `n_init = 300`
- `horizon = 100`
- `edmd_degree = 12`
- `edmd_inits = 80`
- `edmd_steps = 20`
- `time_step = 0.01`
- `tmax = 12.0`

The benchmark matrices are defined inside the script as:

```python
A = [[0, 0, 1],
     [1, -0.2, 1],
     [0, -1, -1]]

J = [[1, 0],
     [0, 1],
     [0, 0]]

C = [0, -1, 0]
```

The state-dependent post-reset map computed by the script is

```text
g_tau(z) = - J^T exp(A I_tau(z)) J z,
```

where `I_tau(z)` is the first admissible reset time `t >= tau` satisfying

```text
C exp(A t) J z >= 0.
```

The finite-horizon cocycle estimate reported by the script is the terminal horizon quantity `||W^K||^(1/K)` with `K = 100`. It is intentionally not the maximum over `k <= K`, because the latter measures short transient amplification and may exceed one even when the asymptotic weighted spectral radius is below one.

### Synthetic \(S^2\) benchmark

The synthetic script uses deterministic pseudo-random sampling with fixed seeds.

The induced homogeneous map is

```text
g_mu(z) = mu * rho0(s) * Pi0(s) * ||z||,  s = z / ||z||.
```

The radial gain satisfies

```text
rho_mu(s) = mu * rho0(s),
```

while the angular map `Pi0` is independent of `mu`. Consequently, the weighted operator satisfies

```text
W_{rho_mu} = mu W_{rho0}.
```

For this reason, the script computes the base diagnostics once for `mu = 1` and rescales them across the parameter sweep. This is a mathematical property of the constructed example, not a numerical shortcut.

Default numerical parameters include:

- `mus = (0.90, 1.00, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40)`
- `n_init = 5000` sampled directions for the finite-horizon cocycle diagnostic
- `horizon = 100`
- polynomial EDMD dictionary with all monomials in `(s1, s2, s3)` of total degree at most `3`
- `20` scalar EDMD observables
- `700` initial directions for EDMD
- `70` steps per EDMD trajectory
- pseudoinverse singular-value cutoff `1e-10`

## Interpretation of outputs

### FORE CSV columns

The columns in `fore_benchmark_results.csv` are:

- `tau_m`: minimum dwell time.
- `r_cocycle_K100`: terminal finite-horizon weighted cocycle estimate.
- `transient_max_rho`: maximum one-step radial gain on the angular grid.
- `r_edmd_deg12`: dominant-modulus eigenvalue estimate obtained by weighted EDMD with Fourier degree 12.
- `frozen_exp_minus_2tau`: frozen regular-branch reference value `exp(-2 tau_m)`.
- `min_return_time`: smallest computed return time on the angular grid.
- `max_return_time`: largest computed return time on the angular grid.
- `min_rho`: smallest one-step radial gain on the angular grid.
- `max_rho`: largest one-step radial gain on the angular grid.

### Synthetic \(S^2\) CSV columns

The columns in `synthetic_s2_benchmark_results.csv` are:

- `mu`: scalar radial-gain scaling parameter.
- `rhat_100`: finite-horizon weighted cocycle estimate with horizon `100`.
- `rhat_EDMD`: dominant-modulus weighted EDMD estimate.
- `max_rho`: maximum sampled one-step radial gain.
- `interpretation`: qualitative stability interpretation based on the diagnostic estimates.

## Expected output summary

With the default parameters, the FORE benchmark should show `r_cocycle_K100` and `r_edmd_deg12` below one for the tested `tau_m` values. The one-step gain `max_rho` may be above one; this indicates transient single-step radial amplification and is not the same as the asymptotic weighted spectral-radius estimate.

With the default parameters, the synthetic \(S^2\) benchmark should show a transition of the diagnostic estimates across the stability threshold `1` as `mu` increases. The one-step gain `max_rho` remains above one throughout the tested range, illustrating that one-step amplification and long-run weighted spectral stability are different notions.

## Numerical status

The quantities reported by both scripts are reproducible numerical diagnostics. They are not advertised as rigorous interval-certified upper bounds on the exact spectral radius. The manuscript discusses how certified finite-horizon upper bounds could be obtained in the one-dimensional projective case by validated interval enclosures.

## Archived version

The archived version of this repository is available on Zenodo:

```text
DOI: 10.5281/zenodo.20682315
```

https://doi.org/10.5281/zenodo.20682315

If you add `run_synthetic_s2_benchmark.py` after an earlier Zenodo release, create a new GitHub release and archive that release in Zenodo. If the manuscript cites a version-specific DOI, update the manuscript to cite the DOI corresponding to the release that contains both scripts. If the manuscript cites the concept DOI, make sure the latest Zenodo version contains both scripts.

## Recommended citation

Suggested software citation format:

```text
Mulero-Martínez, J. I. FORE-reset-benchmark: Reproducible Python implementation of the FORE zero-crossing reset benchmark, with weighted cocycle and EDMD estimates for spectral stability analysis. Zenodo, 2026. DOI: 10.5281/zenodo.20682315.
```

## Open-code deposit notes

For the UPCT open-code workflow, the repository should be public before linking or archiving it through Zenodo. Before final deposit, check that the repository contains a license, this README file, and the exact script versions used to generate the numerical results in the manuscript.

## License

This repository is released under the MIT License. See the `LICENSE` file for details.

## Contact

Author/contact: Juan I. Mulero-Martínez, Department of Automation, Electrical Engineering, and Electronic Technology, Technical University of Cartagena (UPCT), Spain. juan.mulero@upct.es

## Version history

- `v1.0`: Initial public deposit version corresponding to `run_fore_benchmark.py`.
- `v1.1`: Added `run_synthetic_s2_benchmark.py` for the synthetic \(S^2\) weighted-reset benchmark.
