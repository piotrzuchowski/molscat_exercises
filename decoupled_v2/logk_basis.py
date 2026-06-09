#!/usr/bin/env python3
"""Conditioned basis for log k(T) expansions in the variable x = T^(-1/6).

You expand   log k(T) = sum_{n=0}^{deg} a_n * x^n,   x = T^(-1/6).
Over a narrow window (e.g. 5-50 K -> x in [0.52, 0.77]) the monomials x^n are
nearly collinear: cond ~ 1e7 and the a_n are perfectly correlated, so they are
useless as ML targets even though they reconstruct the curve fine.

This module fits the SAME polynomial in an orthogonal Chebyshev basis on the
window (cond ~ 4, decorrelated coefficients c_m) -- use the c_m as regression
targets/features -- and converts losslessly back to the monomial a_n for
reporting.  The Chebyshev <-> monomial map is a fixed linear change of basis.

IMPORTANT: build ONE LogKBasis from your global temperature grid and reuse it
for every curve.  Because the T-grid is the same for all systems, the domain is
identical, so the c_m are directly comparable across the dataset (which is what
makes them valid regression targets).
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial import Chebyshev, Polynomial


def x_of_T(T):
    """Expansion variable x = T^(-1/6) (physical: -C6/R^6 capture scaling)."""
    return np.asarray(T, float) ** (-1.0 / 6.0)


class LogKBasis:
    """Fit / evaluate / convert log k(T) expansions in x = T^(-1/6)."""

    def __init__(self, T_grid, deg=5):
        self.T_grid = np.asarray(T_grid, float)
        self.x_grid = x_of_T(self.T_grid)
        self.deg = int(deg)
        # Fixed window -> Chebyshev maps [xmin, xmax] to [-1, 1] internally.
        self.domain = [float(self.x_grid.min()), float(self.x_grid.max())]

    # ----- fitting: curve -> decorrelated Chebyshev coefficients ----------- #
    def fit(self, logk, weights=None):
        """log k values on T_grid  ->  Chebyshev coeffs c (length deg+1).

        `weights` are per-point (e.g. 1/sigma); pass None for an unweighted fit.
        These c_m are the regression targets.
        """
        logk = np.asarray(logk, float)
        return Chebyshev.fit(self.x_grid, logk, self.deg,
                             domain=self.domain, w=weights).coef

    def fit_many(self, logk_matrix, weights=None):
        """(n_curves, n_T) matrix of log k  ->  (n_curves, deg+1) of c_m."""
        logk_matrix = np.atleast_2d(np.asarray(logk_matrix, float))
        return np.vstack([self.fit(row, weights) for row in logk_matrix])

    # ----- evaluation: coefficients -> curve ------------------------------- #
    def eval(self, c, T=None):
        """Evaluate log k from Chebyshev coeffs c at temperatures T (default grid)."""
        x = self.x_grid if T is None else x_of_T(T)
        return Chebyshev(np.asarray(c, float), domain=self.domain)(x)

    def k_of_T(self, c, T=None):
        """Rate k(T) = exp(log k) from Chebyshev coeffs c."""
        return np.exp(self.eval(c, T))

    # ----- lossless basis change (for reporting) --------------------------- #
    def to_monomial(self, c):
        """Chebyshev coeffs c  ->  monomial a_n with  log k = sum_n a_n x^n."""
        cheb = Chebyshev(np.asarray(c, float), domain=self.domain)
        poly = cheb.convert(kind=Polynomial, domain=self.domain, window=self.domain)
        return poly.coef

    def from_monomial(self, a):
        """Monomial a_n (powers of x)  ->  Chebyshev coeffs c (legacy fits)."""
        poly = Polynomial(np.asarray(a, float), domain=self.domain, window=self.domain)
        cheb = poly.convert(kind=Chebyshev, domain=self.domain)
        return cheb.coef

    # ----- diagnostics ----------------------------------------------------- #
    def condition_numbers(self):
        """(cond_monomial, cond_chebyshev) of the design matrices on this grid."""
        xt = (2 * self.x_grid - sum(self.domain)) / (self.domain[1] - self.domain[0])
        M = np.vander(self.x_grid, self.deg + 1, increasing=True)
        C = np.polynomial.chebyshev.chebvander(xt, self.deg)
        return np.linalg.cond(M), np.linalg.cond(C)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    # Self-test: round-trips and conditioning on a 5-50 K grid.
    T = np.linspace(5.0, 50.0, 15)
    basis = LogKBasis(T, deg=5)

    # synthetic generalized-Arrhenius curve: log k = b0 + b1 x + b2 x^2
    x = x_of_T(T)
    logk_true = -23.0 + 3.1 * x - 1.7 * x ** 2 + 0.4 * x ** 3
    c = basis.fit(logk_true)

    logk_fit = basis.eval(c)
    a = basis.to_monomial(c)              # report-friendly a_n in x
    c_back = basis.from_monomial(a)       # inverse

    cm, cc = basis.condition_numbers()
    print(f"cond(monomial)  = {cm:.2e}")
    print(f"cond(chebyshev) = {cc:.2e}")
    print(f"fit reproduces curve : {np.allclose(logk_fit, logk_true, atol=1e-9)}")
    print(f"cheb -> a_n -> cheb  : {np.allclose(c_back, c, atol=1e-9)}")
    print(f"a_n -> curve matches : "
          f"{np.allclose(sum(a[n]*x**n for n in range(len(a))), logk_true, atol=1e-9)}")
    print(f"Chebyshev coeffs c_m = {np.array2string(c, precision=4)}")
    print(f"monomial  coeffs a_n = {np.array2string(a, precision=4)}")
