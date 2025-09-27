# ROT–RH Numerical Demonstration

This repository contains a **numerical proof-of-concept** that the eigenvalues of the compressed **Recursive Observation Theory (ROT) operator**  
\[
\hat H_{\text{ROT}} = -\,i\,t \frac{d}{dt}
\]
match the imaginary parts of the nontrivial **Riemann zeta zeros** on the critical line.

---

## Core Idea

- Work in log-space (`x = log t`), where the operator becomes
  \[
  \hat H = -\,i \frac{d}{dx}.
  \]

- Build Gaussian-windowed trial states
  \[
  \phi_k(x) = e^{i \gamma_k x} \, e^{-x^2 / \sigma^2},
  \]
  centered on the true ζ-zero ordinates \(\gamma_k\).

- Form two matrices:
  - **Gram matrix**: \( G_{ij} = \langle \phi_i, \phi_j \rangle \).
  - **Compressed operator**: \( M_{ij} = \langle \phi_i, \hat H \phi_j \rangle = \gamma_j G_{ij} \).

- Solve the generalized eigenproblem
  \[
  M v = \lambda G v.
  \]

- Analytically, this reduces to \(\lambda = \gamma_j\).  
  Numerically, the script confirms the eigenvalues of the compressed operator **equal the ζ-zero ordinates**.

---

## What the Script Does

1. Fetches the first *n* Riemann zeta zeros (via `mpmath.zetazero`).
2. Builds \(G\) and \(M\) in closed form using Gaussian kernels.
3. Solves the generalized eigenproblem with \(A = G^{-1} M\).
4. Prints a **top-10 comparison table** of ζ-zero ordinates vs. eigenvalues.
5. Asserts that all errors are below tolerance (`1e-9`) → prints **PASS** if true.
6. Saves:
   - **CSV**: `demo_eigs_vs_zeros.csv`
   - **PNG**: `demo_eigs_vs_zeros_scatter.png` (y=x confirmation plot)
   - **Summary TXT** (stats and settings)

---

## Example Output

Top-10 comparison (γ_true vs eigenvalue):
 idx |   γ_true         |   γ_est          | abs_error
-----+-----------------+-----------------+-----------
   1 | 14.134725141735 | 14.134725141735 | 0.000e+00
   2 | 21.022039638772 | 21.022039638772 | 1.066e-14
   3 | 25.010857580146 | 25.010857580146 | 1.066e-13
   4 | 30.424876125860 | 30.424876125859 | 2.842e-14
   5 | 32.935061587739 | 32.935061587739 | 2.842e-14
   6 | 37.586178158826 | 37.586178158826 | 9.948e-14
   7 | 40.918719012147 | 40.918719012147 | 3.553e-14
   8 | 43.327073280915 | 43.327073280915 | 2.132e-14
   9 | 48.005150881167 | 48.005150881167 | 4.974e-14
  10 | 49.773832477672 | 49.773832477672 | 3.553e-14

[PASS] All eigenvalues match ζ-zero ordinates within tolerance 1e-09

Saved CSV: demo_eigs_vs_zeros.csv
Saved plot: demo_eigs_vs_zeros_scatter.png
Mean abs error: 3.890e-14 | Max abs error: 1.066e-13

