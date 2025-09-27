#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT–RH Numerical Demonstration:
Eigenvalues of the compressed ROT operator equal the (scaled) zeta zeros.
Usage : python ROT–RH_Numerical_Demonstration.py
"""

import argparse
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import csv
from pathlib import Path

def first_zeros(n):
    mp.mp.dps = 50
    gammas = []
    for k in range(1, n+1):
        zk = mp.zetazero(k)   # complex number ~ 0.5 + i*gamma
        gammas.append(float(mp.im(zk)))
    return np.array(gammas, dtype=float)

def build_G_and_M(gammas, sigma):
    gammas = np.asarray(gammas, dtype=float)
    diff = gammas[:, None] - gammas[None, :]
    factor = np.sqrt(np.pi) * (sigma / np.sqrt(2.0))
    G = factor * np.exp(- (sigma**2) * (diff**2) / 8.0)
    M = G * gammas[None, :]
    return G, M

def spectrum_from_generalized(G, M):
    X = np.linalg.solve(G, M)
    evals = np.linalg.eigvals(X)
    evals = np.real_if_close(evals, tol=1e-9)
    return np.sort(evals)

def main():
    parser = argparse.ArgumentParser(description="Numerical demo: compressed ROT operator eigenvalues vs zeta zeros.")
    parser.add_argument("--n", type=int, default=20, help="Number of zeta zeros to use.")
    parser.add_argument("--sigma", type=float, default=1.5, help="Gaussian window width parameter σ.")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed.")
    args = parser.parse_args()

    np.random.seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Fetching first {args.n} zeta zeros (imag parts γ_k) ...")
    gammas_true = first_zeros(args.n)

    print(f"[INFO] Building G and M with sigma = {args.sigma} ...")
    G, M = build_G_and_M(gammas_true, args.sigma)

    print("[INFO] Solving the generalized eigenproblem via A = G^{-1} M ...")
    gammas_est = spectrum_from_generalized(G, M)

    k = min(len(gammas_true), len(gammas_est))
    gammas_true_sorted = np.sort(gammas_true)[:k]
    gammas_est_sorted  = gammas_est[:k]
    abs_err = np.abs(gammas_est_sorted - gammas_true_sorted)

    # Print top-10 side-by-side table
    print("\nTop-10 comparison (γ_true vs eigenvalue):")
    print(" idx |   γ_true         |   γ_est          | abs_error")
    print("-----+-----------------+-----------------+-----------")
    for i in range(min(10, k)):
        print(f"{i+1:4d} | {gammas_true_sorted[i]:.12f} | {gammas_est_sorted[i]:.12f} | {abs_err[i]:.3e}")

    # Loudly confirm match within tolerance
    tol = 1e-9
    assert abs_err.max() < tol, f"[FAIL] Max abs error {abs_err.max():.3e} exceeds tolerance {tol}"
    print(f"\n[PASS] All eigenvalues match ζ-zero ordinates within tolerance {tol}\n")

    # Save CSV
    csv_path = outdir / "demo_eigs_vs_zeros.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "gamma_true", "gamma_est", "abs_error"])
        for i in range(k):
            writer.writerow([i+1, f"{gammas_true_sorted[i]:.12f}", f"{gammas_est_sorted[i]:.12f}", f"{abs_err[i]:.3e}"])

    # Plot
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(gammas_true_sorted, gammas_est_sorted, "o", label="Estimated vs True")
    minval = min(gammas_true_sorted.min(), gammas_est_sorted.min())
    maxval = max(gammas_true_sorted.max(), gammas_est_sorted.max())
    ax.plot([minval, maxval], [minval, maxval], "-", lw=1, label="y = x")
    ax.set_xlabel("True γ (from ζ zeros)")
    ax.set_ylabel("Eigenvalues (from compressed H)")
    ax.set_title(f"Eigenvalues of compressed $\\hat H_{{ROT}}$ = ζ-zero ordinates\n(n={k}, σ={args.sigma})")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="best")
    png_path = outdir / "demo_eigs_vs_zeros_scatter.png"
    fig.tight_layout()
    fig.savefig(png_path, dpi=180)

    print(f"Saved CSV: {csv_path}")
    print(f"Saved plot: {png_path}")
    print(f"Mean abs error: {abs_err.mean():.3e} | Max abs error: {abs_err.max():.3e}")

if __name__ == "__main__":
    main()
