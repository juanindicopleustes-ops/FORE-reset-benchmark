"""Reproducible numerical benchmark for the FORE zero-crossing reset map.

Outputs:
  * fore_benchmark_results_v3.csv
  * nahs_fore_benchmark_v3.pdf
  * nahs_fore_benchmark_v3.png

The code computes the state-dependent post-reset map
    g_tau(z) = -J^T exp(A I_tau(z)) J z,
where I_tau(z) is the first admissible reset time t >= tau satisfying
    C exp(A t) J z >= 0.
The reported spectral-radius estimates are diagnostic, not certified bounds.
"""
from __future__ import annotations
import csv
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from numpy.linalg import norm
from scipy.linalg import eig
import matplotlib.pyplot as plt

A = np.array([[0.0, 0.0, 1.0], [1.0, -0.2, 1.0], [0.0, -1.0, -1.0]])
J = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
C = np.array([0.0, -1.0, 0.0])

@dataclass(frozen=True)
class Config:
    taus: tuple[float, ...] = (0.10, 0.20, 0.40, 0.60, 0.80)
    n_tab: int = 200
    n_init: int = 300
    horizon: int = 100
    edmd_degree: int = 12
    edmd_inits: int = 80
    edmd_steps: int = 20
    time_step: float = 0.01
    tmax: float = 12.0
    out_dir: Path = Path("/mnt/data")

class ReturnMapEvaluator:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.t_grid = np.arange(0.0, cfg.tmax + cfg.time_step, cfg.time_step)
        lam, V = eig(A)
        Vinv = np.linalg.inv(V)
        E = np.exp(np.outer(self.t_grid, lam))  # Nt x 3
        # M(t)=J^T V diag(exp(lam t)) Vinv J = sum_i exp(lam_i t) a_i b_i^T
        self.r_grid = np.empty((len(self.t_grid), 2))
        self.M_grid = np.empty((len(self.t_grid), 2, 2))
        CV = C @ V
        VinvJ = Vinv @ J
        for k in range(len(self.t_grid)):
            D = E[k]
            self.r_grid[k] = np.real((CV * D) @ VinvJ)
            self.M_grid[k] = np.real(J.T @ (V * D) @ VinvJ)

    def interp_row(self, t: float) -> np.ndarray:
        return np.array([np.interp(t, self.t_grid, self.r_grid[:, j]) for j in range(2)])

    def interp_M(self, t: float) -> np.ndarray:
        M = np.empty((2, 2))
        for i in range(2):
            for j in range(2):
                M[i, j] = np.interp(t, self.t_grid, self.M_grid[:, i, j])
        return M

    def first_return_time(self, z: np.ndarray, tau: float) -> float:
        i0 = int(np.ceil(tau / self.cfg.time_step))
        f_tau = float(self.interp_row(tau) @ z)
        if f_tau >= 0.0:
            return float(tau)
        vals = self.r_grid[i0:, :] @ z
        hit = np.flatnonzero(vals >= 0.0)
        if len(hit) == 0:
            raise RuntimeError(f"No crossing before {self.cfg.tmax} for tau={tau}, z={z}")
        ir = i0 + int(hit[0])
        t_right = float(self.t_grid[ir]); f_right = float(vals[int(hit[0])])
        t_left = max(float(tau), float(self.t_grid[ir - 1])); f_left = float(self.interp_row(t_left) @ z)
        if abs(f_right - f_left) < 1e-15:
            return t_right
        alpha = -f_left / (f_right - f_left)
        return float(t_left + alpha * (t_right - t_left))

    def g(self, z: np.ndarray, tau: float) -> tuple[np.ndarray, float]:
        T = self.first_return_time(z, tau)
        return -self.interp_M(T) @ z, T

def angle_of(v: np.ndarray) -> float:
    return float(np.mod(np.arctan2(v[1], v[0]), 2.0 * np.pi))

def build_lookup(tau: float, cfg: Config, ev: ReturnMapEvaluator):
    theta = np.linspace(0.0, 2.0 * np.pi, cfg.n_tab, endpoint=False)
    rho = np.empty_like(theta); theta_next = np.empty_like(theta); returns = np.empty_like(theta)
    for i, th in enumerate(theta):
        z = np.array([np.cos(th), np.sin(th)])
        gz, T = ev.g(z, tau)
        rho[i] = norm(gz); theta_next[i] = angle_of(gz); returns[i] = T
    return theta, rho, theta_next, returns

def periodic_interp(x, grid, values):
    # Fast periodic linear interpolation on a uniform grid over [0, 2*pi).
    n = len(values)
    u = (np.mod(x, 2.0 * np.pi) / (2.0 * np.pi)) * n
    i = int(np.floor(u)) % n
    a = u - np.floor(u)
    return (1.0 - a) * values[i] + a * values[(i + 1) % n]

def unwrap_map(values):
    return np.unwrap(values)

def cocycle_estimate(grid, rho, theta_next, cfg):
    """Finite-horizon Gelfand estimate at the terminal horizon K.

    We intentionally report ||W^K||^(1/K), not max_{k<=K} ||W^k||^(1/k),
    because the latter measures short transient amplification and may exceed one
    even when the asymptotic spectral radius is below one.
    """
    theta_next_u = unwrap_map(theta_next)
    best_log = -np.inf
    for th0 in np.linspace(0.0, 2.0 * np.pi, cfg.n_init, endpoint=False):
        th = th0; logprod = 0.0
        for _ in range(cfg.horizon):
            r = float(periodic_interp(th, grid, rho))
            logprod += np.log(max(r, 1e-300))
            th = float(periodic_interp(th, grid, theta_next_u))
        best_log = max(best_log, logprod / cfg.horizon)
    return float(np.exp(best_log))

def fourier_features(theta, degree):
    rows = [np.ones_like(theta)]
    for k in range(1, degree + 1):
        rows.extend([np.cos(k * theta), np.sin(k * theta)])
    return np.vstack(rows)

def edmd_estimate(grid, rho, theta_next, cfg):
    theta_next_u = unwrap_map(theta_next)
    xs, ys, ws = [], [], []
    for th0 in np.linspace(0.0, 2.0 * np.pi, cfg.edmd_inits, endpoint=False):
        th = th0
        for _ in range(cfg.edmd_steps):
            r = float(periodic_interp(th, grid, rho)); th_y = float(periodic_interp(th, grid, theta_next_u))
            xs.append(np.mod(th, 2.0 * np.pi)); ys.append(np.mod(th_y, 2.0 * np.pi)); ws.append(r)
            th = th_y
    X = fourier_features(np.asarray(xs), cfg.edmd_degree)
    Y = fourier_features(np.asarray(ys), cfg.edmd_degree) * np.asarray(ws)[None, :]
    W = Y @ np.linalg.pinv(X, rcond=1e-10)
    return float(np.max(np.abs(np.linalg.eigvals(W))))

def main():
    cfg = Config(); cfg.out_dir.mkdir(parents=True, exist_ok=True); ev = ReturnMapEvaluator(cfg)
    rows = []
    for tau in cfg.taus:
        print(f"tau={tau:.2f}")
        grid, rho, theta_next, returns = build_lookup(tau, cfg, ev)
        row = {"tau_m": tau,
               "r_cocycle_K100": cocycle_estimate(grid, rho, theta_next, cfg),
               "transient_max_rho": float(np.max(rho)),
               "r_edmd_deg12": edmd_estimate(grid, rho, theta_next, cfg),
               "frozen_exp_minus_2tau": float(np.exp(-2.0 * tau)),
               "min_return_time": float(np.min(returns)),
               "max_return_time": float(np.max(returns)),
               "min_rho": float(np.min(rho)),
               "max_rho": float(np.max(rho))}
        rows.append(row); print(row)
    csv_path = cfg.out_dir / "fore_benchmark_results_v3.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys())); writer.writeheader(); writer.writerows(rows)
    taus = np.array([r["tau_m"] for r in rows])
    plt.figure(figsize=(6.6, 4.2))
    plt.plot(taus, [r["r_cocycle_K100"] for r in rows], marker="o", label=r"terminal cocycle $\widehat r_{100}$")
    plt.plot(taus, [r["r_edmd_deg12"] for r in rows], marker="s", label=r"weighted EDMD $\widehat r_{\rm EDMD}$")
    plt.plot(taus, [r["max_rho"] for r in rows], marker="x", linestyle=":", label=r"one-step gain $\max \rho$")
    plt.plot(taus, [r["frozen_exp_minus_2tau"] for r in rows], marker="^", label=r"frozen branch $e^{-2\tau_m}$")
    plt.axhline(1.0, linestyle="--", linewidth=1.0, label="stability threshold")
    plt.axvline(0.6145, linestyle=":", linewidth=1.0, label="reported LMI/SOS threshold")
    plt.xlabel(r"minimum dwell time $\tau_m$"); plt.ylabel("spectral-radius estimate")
    plt.ylim(0.0, 1.25); plt.grid(True, alpha=0.3); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(cfg.out_dir / "nahs_fore_benchmark_v3.pdf"); plt.savefig(cfg.out_dir / "nahs_fore_benchmark_v3.png", dpi=200)
    print(f"Wrote {csv_path}")
if __name__ == "__main__": main()
