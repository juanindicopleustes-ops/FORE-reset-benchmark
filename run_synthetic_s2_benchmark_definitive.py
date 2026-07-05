#!/usr/bin/env python3
"""
Synthetic S^2 weighted-reset benchmark.

This script generates a controlled three-dimensional homogeneous induced map

    g_mu(z) = mu * rho0(s) * Pi0(s) * ||z||,   s = z/||z||,

where Pi0:S^2->S^2 is a smooth nonlinear angular map and rho0:S^2->(0,infty)
is a positive radial gain. Since the angular map is independent of mu and the
radial gain scales as rho_mu = mu*rho0, the weighted operator satisfies

    W_{rho_mu} = mu W_{rho0}.

The example is intended to illustrate a spectral stability transition in a
post-reset dimension m=3, using both a finite-horizon cocycle diagnostic and a
projected-data EDMD diagnostic.

Outputs:
    synthetic_s2_benchmark_results.csv
    synthetic_s2_weighted_benchmark.pdf

Dependencies:
    numpy, matplotlib

Author: Juan I. Mulero-Martinez
License: MIT
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

OUTDIR = Path(__file__).resolve().parent


def rodrigues_rotation(axis: Iterable[float], angle: float) -> np.ndarray:
    """Return the 3D rotation matrix with given axis and angle."""
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    x, y, z = a
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    I = np.eye(3)
    return I + math.sin(angle) * K + (1.0 - math.cos(angle)) * (K @ K)


R = rodrigues_rotation(axis=(1.0, 2.0, 3.0), angle=0.83)


def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n <= 1e-14:
        raise FloatingPointError("attempted to normalize a near-zero vector")
    return v / n


def angular_map(s: np.ndarray) -> np.ndarray:
    """Smooth nonlinear angular map Pi0 on S^2."""
    s = np.asarray(s, dtype=float)
    v = np.array(
        [
            s[0] + 0.35 * s[1] * s[2],
            s[1] + 0.30 * math.sin(math.pi * s[0]) * s[2],
            s[2] + 0.25 * s[0] * s[1],
        ],
        dtype=float,
    )
    return normalize(R @ v)


def rho_base(s: np.ndarray) -> float:
    """Positive base radial gain rho0 on S^2.

    The gain has a narrow high-gain sector, so max rho0 can exceed one-step
    stability bounds after scaling, while the long-run cocycle can remain below
    one for a range of mu.
    """
    s = np.asarray(s, dtype=float)
    return float(
        0.78
        + 0.50 * math.exp(12.0 * (s[0] - 1.0))
        + 0.04 * (0.5 + 0.5 * math.sin(3.0 * s[1] + 2.0 * s[2]))
    )


def rho_mu(s: np.ndarray, mu: float) -> float:
    return mu * rho_base(s)


def random_sphere(n: int, rng: np.random.Generator) -> np.ndarray:
    X = rng.normal(size=(n, 3))
    X /= np.linalg.norm(X, axis=1)[:, None]
    return X


def finite_horizon_cocycle(mu: float, n_init: int = 5000, horizon: int = 100, seed: int = 1) -> float:
    """Compute max_s prod_{j=0}^{K-1} rho_mu(Pi0^j(s))^(1/K) over random samples."""
    rng = np.random.default_rng(seed)
    S = random_sphere(n_init, rng)
    best = 0.0
    for s0 in S:
        s = s0.copy()
        prod = 1.0
        for _ in range(horizon):
            prod *= rho_mu(s, mu)
            s = angular_map(s)
        best = max(best, prod ** (1.0 / horizon))
    return best


def max_one_step_gain(mu: float, n_samples: int = 50000, seed: int = 7) -> float:
    rng = np.random.default_rng(seed)
    S = random_sphere(n_samples, rng)
    return max(rho_mu(s, mu) for s in S)


def monomial_powers(max_degree: int) -> List[Tuple[int, int, int]]:
    powers: List[Tuple[int, int, int]] = []
    for total in range(max_degree + 1):
        for a in range(total + 1):
            for b in range(total - a + 1):
                c = total - a - b
                powers.append((a, b, c))
    return powers


def feature_matrix(S: np.ndarray, powers: List[Tuple[int, int, int]]) -> np.ndarray:
    S = np.atleast_2d(S)
    x, y, z = S[:, 0], S[:, 1], S[:, 2]
    rows = []
    for a, b, c in powers:
        rows.append((x**a) * (y**b) * (z**c))
    return np.vstack(rows)


def weighted_edmd_radius(
    mu: float,
    degree: int = 3,
    n_init: int = 700,
    n_steps: int = 70,
    seed: int = 11,
) -> float:
    """Projected-data EDMD estimate of the weighted composition operator radius."""
    rng = np.random.default_rng(seed)
    S0 = random_sphere(n_init, rng)
    X = []
    Y = []
    W = []
    for s0 in S0:
        s = s0.copy()
        for _ in range(n_steps):
            sn = angular_map(s)
            X.append(s)
            Y.append(sn)
            W.append(rho_mu(s, mu))
            s = sn
    X = np.asarray(X)
    Y = np.asarray(Y)
    W = np.asarray(W)
    powers = monomial_powers(degree)
    PhiX = feature_matrix(X, powers)
    PhiY = feature_matrix(Y, powers) * W[None, :]
    Ahat = PhiY @ np.linalg.pinv(PhiX, rcond=1e-10)
    eigvals = np.linalg.eigvals(Ahat)
    return float(np.max(np.abs(eigvals)))


def main() -> None:
    mus = np.array([0.90, 1.00, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40])

    # The parameter sweep below is not a numerical shortcut unrelated to the
    # model. In this benchmark the radial gain is defined as
    #     rho_mu(s) = mu * rho_0(s)
    # while the angular map Pi_0 is independent of mu. Hence the weighted
    # composition operator satisfies W_mu = mu * W_0. All spectral diagnostics
    # associated with the weighted operator therefore scale linearly with mu.
    # We compute the base diagnostics once at mu = 1 and rescale them across
    # the parameter sweep, exactly reflecting this analytical scaling property.
    base_r100 = finite_horizon_cocycle(1.0)
    base_edmd = weighted_edmd_radius(1.0)
    base_maxrho = max_one_step_gain(1.0)

    rows = []
    for mu in mus:
        rows.append(
            {
                "mu": mu,
                "rhat_100": mu * base_r100,
                "rhat_EDMD": mu * base_edmd,
                "max_rho": mu * base_maxrho,
                "interpretation": "stable" if mu * base_edmd < 1.0 else "unstable",
            }
        )

    csv_path = OUTDIR / "synthetic_s2_benchmark_results.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["mu", "rhat_100", "rhat_EDMD", "max_rho", "interpretation"]
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Base finite-horizon rhat_100: {base_r100:.6f}")
    print(f"Base EDMD radius:             {base_edmd:.6f}")
    print(f"Base max rho:                 {base_maxrho:.6f}")
    print(f"CSV written to: {csv_path}")

    fig_path = OUTDIR / "synthetic_s2_weighted_benchmark.pdf"
    plt.figure(figsize=(6.6, 4.4))
    plt.plot(mus, [r["rhat_100"] for r in rows], marker="o", label=r"finite-horizon $\widehat r_{100}$")
    plt.plot(mus, [r["rhat_EDMD"] for r in rows], marker="s", label=r"weighted EDMD $\widehat r_{\rm EDMD}$")
    plt.plot(mus, [r["max_rho"] for r in rows], marker="^", label=r"one-step gain $\max\rho_\mu$")
    plt.axhline(1.0, linestyle="--", label="stability threshold")
    plt.xlabel(r"reset-gain scale $\mu$")
    plt.ylabel("spectral diagnostic")
    plt.title(r"Synthetic $S^2$ weighted reset benchmark")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(fig_path)
    print(f"Figure written to: {fig_path}")


if __name__ == "__main__":
    main()
