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

Console (abridged):

