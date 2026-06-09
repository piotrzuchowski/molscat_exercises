"""Orthonormal (Chebyshev) basis for log k(T) curves on a fixed global T grid.

log k is expanded in x = T^(-1/6). The raw monomial basis {x^n} is severely
ill-conditioned over a narrow x-range, so its coefficients are correlated and
unstable to regress. This builds an orthonormal basis over the supplied T grid
(Chebyshev design matrix in the mapped variable, then QR), giving coefficients
c_m that are DECORRELATED (Q^T Q = I) and ordered by significance -> safe to
regress with independent per-mode models. Mapping back to the monomial a_n
convention (log k = sum_n a_n T^(-n/6)) is lossless.

Conventions: works in natural log (np.log), matching the usual a_n reporting.
"""

import numpy as np
from numpy.polynomial import chebyshev as Ch, polynomial as Poly


class LogKBasis:
    def __init__(self, T, deg=5):
        self.T = np.asarray(T, float)
        self.deg = int(deg)
        self.x = self.T ** (-1.0 / 6.0)            # expansion variable
        self.domain = (float(self.x.min()), float(self.x.max()))
        a, b = self.domain
        u = (2.0 * self.x - (a + b)) / (b - a)     # map x-> [-1,1]
        V = Ch.chebvander(u, self.deg)             # (nT, deg+1): T_m(u)
        Q, R = np.linalg.qr(V)
        s = np.sign(np.diag(R)); s[s == 0] = 1.0   # fix sign so c_m comparable
        self.Q = Q * s
        self.R = R * s[:, None]

    # --- curves -> decorrelated coefficients ------------------------------ #
    def fit_many(self, logk):
        """logk: (n_curves, nT) natural-log values -> C: (n_curves, deg+1)."""
        return np.asarray(logk, float) @ self.Q

    def fit(self, logk_row):
        return np.asarray(logk_row, float) @ self.Q

    # --- coefficients -> curves ------------------------------------------- #
    def eval(self, c):
        """c: (deg+1,) or (n, deg+1) -> log k at the grid T."""
        return np.asarray(c, float) @ self.Q.T

    def _cheb_coeffs(self, c):
        return np.linalg.solve(self.R, np.asarray(c, float))

    def k_of_T(self, c, T_fine):
        """Evaluate k (not log) at arbitrary temperatures from coefficients."""
        b = self._cheb_coeffs(c)
        a, bb = self.domain
        xf = np.asarray(T_fine, float) ** (-1.0 / 6.0)
        uf = (2.0 * xf - (a + bb)) / (bb - a)
        return np.exp(Ch.chebval(uf, b))

    # --- lossless conversion to the monomial a_n convention --------------- #
    def to_monomial(self, c):
        """a_n such that  log k = sum_n a_n * T^(-n/6)  (natural log)."""
        b = self._cheb_coeffs(c)
        a, bb = self.domain
        series = Ch.Chebyshev(b, domain=[a, bb], window=[-1.0, 1.0])
        return series.convert(kind=Poly.Polynomial).coef
