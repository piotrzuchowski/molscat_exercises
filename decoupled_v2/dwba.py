#!/usr/bin/env python3
r"""Semiclassical (WKB) distorted-wave Born approximation for rotational
(de)excitation of a linear rigid rotor by an atom -- the cheap analytic cousin
of the full close-coupling MOLSCAT run, intended as a fast estimator and as a
physics-based mean function for a k(T) regressor.

Model potential (Born-Mayer + dispersion, unnormalized Legendre):
    V(R,theta) = V0(R) + V1(R) P1 + V2(R) P2
    V0 = A e^{-aR} - C6 R^{-6}            (isotropic -> distorting potential)
    V1 = eta*A e^{-aR} - a1*C6 R^{-6}     (lambda=1 -> drives |dj|=1)
    V2 = a2*(A e^{-aR} - C6 R^{-6})       (lambda=2 -> drives |dj|=0,2)

Method (first order in the anisotropy V_lambda, exact angular coupling):
  * distorted waves phi_{k l}(R) are the regular JWKB solutions in V0 with the
    Langer substitution l(l+1) -> (l+1/2)^2, normalized to sin(kR - l pi/2 + d).
  * DWBA K-matrix:  K^J_{jl,j'l'} = -(2mu/hbar^2) f_lambda(jl,j'l';J)
                       * (k_j k_j')^{-1/2} \int phi_{k_j l} V_lambda phi_{k_j' l'} dR
  * S = 1 + 2iK to first order -> |S_{fi}|^2 = 4|K|^2 (f != i)
  * sigma(j->j') = (4 pi / (k_j^2 (2j+1))) sum_J (2J+1) sum_{l l'} |K^J|^2
  * k(T) by Maxwell-Boltzmann average over the initial collision energy.

f_lambda is the Percival-Seaton coefficient
    f = (-1)^{j+j'-J} [(2j+1)(2j'+1)(2l+1)(2l'+1)]^{1/2}
        ( j j' L; 0 0 0)( l l' L; 0 0 0){ j j' L; l' l J }.

Units: R [Ang], energies [cm^-1], mu [amu], B [cm^-1]; sigma [Ang^2], k [cm^3/s].
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from sympy.physics.wigner import wigner_3j, wigner_6j

HARTREE2CM = 219474.6313702
BOHR2ANG = 0.529177210903
HBAR2_2MU = 16.85762920          # hbar^2/2mu = HBAR2_2MU/mu  [cm^-1 Ang^2]

# k(T) unit constants (SI -> cm^3/s)
KB_CM, CM2J, KB_J, AMU, ANG2 = 0.6950348, 1.98644586e-23, 1.380649e-23, 1.66053907e-27, 1e-20


# --------------------------------------------------------------------------- #
@lru_cache(maxsize=None)
def f_lambda(j, l, jp, lp, J, lam):
    """Percival-Seaton angular coefficient (cached). Returns a float (may be 0)."""
    t1 = wigner_3j(j, jp, lam, 0, 0, 0)
    if t1 == 0:
        return 0.0
    t2 = wigner_3j(l, lp, lam, 0, 0, 0)
    if t2 == 0:
        return 0.0
    s6 = wigner_6j(j, jp, lam, lp, l, J)
    if s6 == 0:
        return 0.0
    pref = (-1) ** (j + jp - J)
    val = pref * np.sqrt((2 * j + 1) * (2 * jp + 1) * (2 * l + 1) * (2 * lp + 1))
    return float(val * t1 * t2 * s6)


def _k_of_T(E_cm, sigma_ang2, T, mu_amu):
    """Boltzmann-average sigma(E) [Ang^2] -> k(T) [cm^3/s]; T may be an array."""
    E, s = E_cm * CM2J, sigma_ang2 * ANG2
    T = np.atleast_1d(np.asarray(T, float))
    out = np.empty_like(T)
    for i, Ti in enumerate(T):
        kT = KB_J * Ti
        integ = trapezoid(s * E * np.exp(-E / kT), E)
        out[i] = np.sqrt(8.0 / (np.pi * mu_amu * AMU)) * kT ** -1.5 * integ * 1e6
    return out


# --------------------------------------------------------------------------- #
class DWBA:
    """WKB distorted-wave Born estimator for one potential + reduced mass."""

    def __init__(self, A_au, a_au, C6_au, a1, a2, eta, B_cm1, mu_amu,
                 R=(1.4, 30.0, 0.01)):
        self.A = A_au * HARTREE2CM             # cm^-1
        self.a = a_au / BOHR2ANG               # Ang^-1
        self.C6 = C6_au * HARTREE2CM * BOHR2ANG ** 6   # cm^-1 Ang^6
        self.a1, self.a2, self.eta = a1, a2, eta
        self.B, self.mu = B_cm1, mu_amu

        self.Rg = np.arange(*R)
        exp = self.A * np.exp(-self.a * self.Rg)
        disp = self.C6 * self.Rg ** -6
        self.V0 = exp - disp
        self.V1 = eta * exp - a1 * disp
        self.V2 = a2 * (exp - disp)
        self.U0 = self.V0 * mu_amu / HBAR2_2MU          # 2mu/hbar^2 * V0 [Ang^-2]
        self.coef = mu_amu / HBAR2_2MU                  # 2mu/hbar^2  [1/(cm^-1 Ang^2)]

    # ----- WKB distorted wave (regular, amplitude-1 asymptotically) -------- #
    def _wave(self, k, l, cache):
        key = (round(float(k), 10), l)
        w = cache.get(key)
        if w is not None:
            return w
        R = self.Rg
        kappa2 = k * k - (l + 0.5) ** 2 / R ** 2 - self.U0
        neg = np.where(kappa2 <= 0.0)[0]
        start = (neg[-1] + 1) if neg.size else 0     # outermost turning point
        phi = np.zeros_like(R)
        if start < R.size - 2:
            kappa = np.sqrt(np.clip(kappa2[start:], (0.05 * k) ** 2, None))
            phase = cumulative_trapezoid(kappa, R[start:], initial=0.0)
            phi[start:] = np.sqrt(k / kappa) * np.sin(phase + 0.25 * np.pi)
        cache[key] = phi
        return phi

    def _kchan(self, Ekin):
        return np.sqrt(Ekin * self.mu / HBAR2_2MU)

    # ----- state-to-state cross section ------------------------------------ #
    def sigma(self, j, jp, Ecoll, Jmax=80, tol=1e-4):
        """sigma(j->j') [Ang^2] at initial collision energy Ecoll [cm^-1]."""
        lam = abs(j - jp)
        if lam not in (1, 2):
            return 0.0                                   # only V1, V2 drive these
        Vlam = self.V1 if lam == 1 else self.V2
        Etot = Ecoll + self.B * j * (j + 1)
        Ekin_f = Etot - self.B * jp * (jp + 1)
        if Ekin_f <= 0.0:                                # final channel closed
            return 0.0
        ki, kf = self._kchan(Ecoll), self._kchan(Ekin_f)

        cache = {}
        sig = 0.0
        run = 0.0
        quiet = 0
        for J in range(0, Jmax + 1):
            sJ = 0.0
            for l in range(abs(J - j), J + j + 1):
                phi_i = self._wave(ki, l, cache)
                for lp in (l - lam, l, l + lam) if lam == 2 else (l - 1, l + 1):
                    if lp < abs(J - jp) or lp > J + jp:
                        continue
                    f = f_lambda(j, l, jp, lp, J, lam)
                    if f == 0.0:
                        continue
                    phi_f = self._wave(kf, lp, cache)
                    I = trapezoid(phi_i * Vlam * phi_f, self.Rg)
                    K = -self.coef * f * I / np.sqrt(ki * kf)
                    sJ += K * K
            sig += (2 * J + 1) * sJ
            run = max(run, (2 * J + 1) * sJ)
            quiet = quiet + 1 if (2 * J + 1) * sJ < tol * run else 0
            if quiet >= 4 and J > j + jp:
                break
        return 4.0 * np.pi / (ki * ki * (2 * j + 1)) * sig

    # ----- thermal rate ---------------------------------------------------- #
    def rate(self, j, jp, T, E_grid=None):
        """k(j->j'; T) [cm^3/s]. T scalar or array."""
        if E_grid is None:
            E_grid = self._default_Egrid()
        sig = np.array([self.sigma(j, jp, E) for E in E_grid])
        return _k_of_T(E_grid, sig, T, self.mu)

    def _default_Egrid(self, T_hi=50.0):
        E_max = 6.0 * KB_CM * T_hi
        return np.unique(np.concatenate([np.geomspace(0.05, 50.0, 35),
                                         np.linspace(50.0, E_max, 12)]))


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    import time
    import pandas as pd

    df = pd.read_csv("potentials_library.csv")
    row = df.iloc[0]
    d = DWBA(row["A"], row["a"], row["C6"], row["a1"], row["a2"], row["eta"],
             row["B_cm1"], row["mu_H2"])
    print(f"id={int(row['id'])}  B={row['B_cm1']:.4g} cm-1  "
          f"eps_90={row['eps_90']:.1f} cm-1  mu={row['mu_H2']:.3f}")

    t0 = time.time()
    for (j, jp) in [(1, 0), (2, 0), (2, 1)]:
        s = d.sigma(j, jp, 50.0)
        print(f"  sigma({j}->{jp}; 50 cm-1) = {s:.4e} Ang^2  (lambda={abs(j-jp)})")
    T = np.array([5.0, 10.0, 20.0, 50.0])
    k10 = d.rate(1, 0, T)
    print(f"  k(1->0) at T={T} K :")
    print("   ", np.array2string(k10, precision=3))
    # detailed-balance check on the reverse (0->1 excitation) via direct calc
    k01 = d.rate(0, 1, T)
    db = k10 * (3.0 / 1.0) * np.exp(-2.0 * d.B / (KB_CM * T))   # g_j=3, g_jp=1
    print(f"  detailed-balance k(0->1) pred vs direct (ratio):")
    print("   ", np.array2string(k01 / db, precision=3))
    print(f"  [{time.time() - t0:.1f}s]")
