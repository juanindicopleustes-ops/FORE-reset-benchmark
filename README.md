# FORE-reset-benchmark

Reproducible Python implementation of the FORE zero-crossing reset benchmark, with weighted cocycle and EDMD estimates for spectral stability analysis.

## Overview

This repository contains the Python code used to reproduce the numerical benchmark associated with the manuscript:

**Homogeneous zero-crossing reset systems and weighted spectral stability criterion**

The code computes the full state-dependent post-reset map of the FORE benchmark, evaluates finite-horizon weighted cocycle estimates, computes projected-data EDMD spectral diagnostics, and generates the CSV table and figure used in the manuscript.

The benchmark is designed to illustrate the weighted spectral-radius mechanism for homogeneous zero-crossing reset maps. The numerical estimates are diagnostic approximations and are not claimed to be rigorous certified upper bounds.

## Repository contents

The repository should contain at least the following files:

- `run_fore_benchmark.py`  
  Main Python script. It computes the numerical estimates and generates the output CSV file and benchmark figure.

- `README.md`  
  This file.

- `LICENSE`  
  MIT License file.

Recommended optional files before public deposit:

- `CITATION.cff` or `codemeta.json`  
  Machine-readable citation metadata for the code repository.

## Generated output files

Running the script generates:

- `fore_benchmark_results.csv`  
  CSV table containing, for each dwell-time value \(\tau_m\), the terminal finite-horizon cocycle estimate, the weighted EDMD estimate, the maximum one-step radial gain, the frozen-branch reference value, return-time ranges, and radial-gain ranges.

- `nahs_fore_benchmark.pdf`  
  Vector figure used in the manuscript.

- `nahs_fore_benchmark.png`  
  Raster version of the same figure.

## Software requirements

The script was written for Python 3 and requires the following Python packages:

- `numpy`
- `scipy`
- `matplotlib`

The script only uses standard-library modules in addition to the packages listed above: `csv`, `dataclasses`, `pathlib`, and `__future__` annotations.

A typical installation command is:

```bash
python -m pip install numpy scipy matplotlib
```

## How to run

From the directory containing `run_fore_benchmark.py`, run:

```bash
python run_fore_benchmark.py
```

### Portability note

The distributed script was generated in a sandbox environment and its default output directory may be set in the `Config` dataclass as:

```python
out_dir: Path = Path("/mnt/data")
```

For normal local, GitHub, or Zenodo use, change this line to:

```python
out_dir: Path = Path(".")
```

or to any desired local output directory before running the script.

## Reproducibility information

The script uses deterministic uniform grids and does not use random sampling. Therefore, no random seed is required.

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

## Interpretation of outputs

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

The quantities `r_cocycle_K100` and `r_edmd_deg12` are diagnostic numerical approximations to the weighted spectral-radius mechanism discussed in the manuscript. They are not advertised as rigorous certified upper bounds.

## Expected output summary

With the default parameters, the generated table should show `r_cocycle_K100` and `r_edmd_deg12` below one for the tested `tau_m` values. The one-step gain `max_rho` may be above one; this indicates transient single-step radial amplification and is not the same as the asymptotic weighted spectral-radius estimate.

## Recommended citation

Until the repository receives a DOI, cite the accompanying manuscript and the public GitHub repository URL. After Zenodo publication, cite the software using the Zenodo DOI assigned to the archived release.

Suggested software citation format:

```text
Author(s). FORE-reset-benchmark: Reproducible Python implementation of the FORE zero-crossing reset benchmark. Version v1.0. Zenodo. DOI: [to be completed after deposit].
```

## Open-code deposit notes

For the UPCT open-code workflow, the repository should be public before linking or archiving it through Zenodo. Before final deposit, check that the repository contains a license, this README file, and the exact script version used to generate the numerical results in the manuscript.

## License

This repository is released under the MIT License.

MIT License

Copyright (c) 2026 [Author(s)]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contact

Author/contact: [add name, affiliation, and email before deposit]

## Version history

- `v1.0`: Initial public deposit version corresponding to `run_fore_benchmark.py`.
