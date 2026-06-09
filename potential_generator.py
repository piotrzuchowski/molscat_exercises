#!/usr/bin/env python3
"""Fictitious potential generator for MOLSCAT ITYPE=1 (atom + linear rigid rotor).

Sobol-samples a family of Born-Mayer + dispersion isotropic potentials with two
leading anisotropies, screens each one for physical sanity, characterizes the
survivors, and writes a CSV library of MOLSCAT-ready coefficients plus
descriptors. See potential_generator_plan.md for the full specification.

Potential model (unnormalized Legendre, MOLSCAT convention), P1 = cos t,
P2 = 3 cos^2 t - 1:
    V0(R) = A*exp(-a*R) - C6*R**-6        (lambda=0, isotropic)
    V1(R) = eta*A*exp(-a*R) - a1*C6*R**-6 (lambda=1, odd anisotropy)
    V2(R) = a2*(A*exp(-a*R) - C6*R**-6)   (lambda=2, quadrupolar)
    V(R, t) = V0 + V1*P1 + V2*P2
            = A*exp(-a*R) * (1 + eta*P1 + a2*P2)
              - C6*R**-6  * (1 + a1*P1  + a2*P2)

The repulsive wall and the dispersion tail carry *independent* odd
anisotropies (eta and a1), so unlike the old single-g(theta) scaling the
well position, wall radius, and depth are all angle-dependent and must be
characterized per angle on the full V(R, theta).

Internal working units are ATOMIC UNITS: R in Bohr, energy in Hartree, so the
dispersion coefficient C6 is in au (Hartree*Bohr^6), where it is easy to anchor
to known systems. Descriptors are reported in spectroscopy units (eps in cm^-1,
R in Angstrom). B in cm^-1, M in amu.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm, qmc

# Unit conversions (CODATA-ish).
HARTREE2CM = 219474.6313702   # Hartree -> cm^-1
BOHR2ANG = 0.529177210903     # Bohr -> Angstrom

# --------------------------------------------------------------------------- #
# CONFIG: ranges, sample count, bounds, units, flags (retune freely)          #
# --------------------------------------------------------------------------- #
# C6 anchored to atomic dispersion coefficients via Slater-Kirkwood combination:
#   C6(Ne-He) ~ 3.0 au  ...  C6(Ne-Xe) ~ 39 au   (the "most likely" vdW window).
CONFIG = {
    # Sampling
    "sobol_m": 11,           # request 2**m candidates (2048 -> ~1200 accepted)
    "seed": 12345,           # Sobol scramble / reproducibility

    # Raw parameter ranges  (lo, hi), all in ATOMIC UNITS.
    # A/a tuned so the derived well-depth bulk lands in the vdW window 15-80 cm^-1
    # (shallow He-Ne-like end ~10, deep Ne-Xe/polar tail toward a few hundred).
    "A":  (150.0, 5000.0),   # log-uniform, repulsive prefactor   [Hartree]
    "a":  (2.2,   2.9),      # log-uniform, repulsive steepness   [Bohr^-1]
    "C6": (3.0,   40.0),     # log-uniform, dispersion coeff      [au, Hartree*Bohr^6]
    "a1": (-0.4,  0.4),      # linear, odd anisotropy of the -C6/R^6 tail
    "a2": (-0.4,  0.4),      # linear, quadrupolar anisotropy (both terms)
    "eta_sigma": 0.2,        # Gaussian (mean 0) odd anisotropy of A*exp(-aR) wall
    "B":  (1.0e-3, 10.0),    # log-uniform, rotational constant   [cm^-1]
    "M":  (8.0,   100.0),    # log-uniform, molecular mass [amu]; >= LiH (~8)

    # Radial grid / root finding  (Bohr)
    "R_min": 4.0,            # inner edge MOLSCAT would integrate from [Bohr]
    "R_max": 30.0,           # outer edge [Bohr]
    "n_grid": 4000,          # radial scan resolution

    # Screening thresholds
    "g_floor": 0.05,         # angular factor must exceed this for all theta
    "eps_min": 5.0,          # min acceptable well depth [cm^-1]
    "eps_max": 500.0,        # max acceptable well depth [cm^-1]

    # Reduced-mass partners
    "m_He": 4.002602,        # He mass [amu]
    "m_H2": 2.01588,         # H2 mass [amu]

    # Optional rough B-M anticorrelation rejection (off by default)
    "couple_BM": False,

    # Output
    "library_csv": "potentials_library.csv",
    "log_csv": "potentials_log.csv",
}

# Special angles (degrees) and their cos(theta).
ANGLES = (0.0, 90.0, 180.0)
COSTHETA = {0.0: 1.0, 90.0: 0.0, 180.0: -1.0}


# --------------------------------------------------------------------------- #
# Core math                                                                   #
# --------------------------------------------------------------------------- #
def v0(R, A, a, C6):
    """Isotropic Born-Mayer + dispersion potential (no damping)."""
    return A * np.exp(-a * R) - C6 * R ** -6


def dv0_dR(R, A, a, C6):
    """Analytic radial derivative of v0."""
    return -a * A * np.exp(-a * R) + 6.0 * C6 * R ** -7


def angular_factor(costheta, a1, a2):
    """Unnormalized-Legendre angular factor 1 + a1*P1 + a2*P2 (dispersion tail)."""
    p1 = costheta
    p2 = 3.0 * costheta ** 2 - 1.0
    return 1.0 + a1 * p1 + a2 * p2


def angular_rep(costheta, eta, a2):
    """Angular factor multiplying the repulsive A*exp(-aR) wall: 1 + eta*P1 + a2*P2."""
    p1 = costheta
    p2 = 3.0 * costheta ** 2 - 1.0
    return 1.0 + eta * p1 + a2 * p2


def v_theta(R, A, a, C6, a1, a2, eta, costheta):
    """Full potential at a fixed angle (au): wall and tail carry separate odd anisotropy."""
    return (A * np.exp(-a * R) * angular_rep(costheta, eta, a2)
            - C6 * R ** -6 * angular_factor(costheta, a1, a2))


def dv_theta_dR(R, A, a, C6, a1, a2, eta, costheta):
    """Analytic radial derivative of v_theta."""
    return (-a * A * np.exp(-a * R) * angular_rep(costheta, eta, a2)
            + 6.0 * C6 * R ** -7 * angular_factor(costheta, a1, a2))


# --------------------------------------------------------------------------- #
# Root finding / characterization of the isotropic V0                         #
# --------------------------------------------------------------------------- #
def _radial_grid(cfg):
    return np.linspace(cfg["R_min"], cfg["R_max"], cfg["n_grid"])


def find_minima(A, a, C6, a1, a2, eta, costheta, cfg, Rgrid=None):
    """Return list of outer minima (R_m) of V(R,theta) on the grid, refined by Brent."""
    if Rgrid is None:
        Rgrid = _radial_grid(cfg)
    args = (A, a, C6, a1, a2, eta, costheta)
    d = dv_theta_dR(Rgrid, *args)
    # A minimum is a sign change of dV/dR from negative to positive.
    sign_change = np.where((d[:-1] < 0.0) & (d[1:] >= 0.0))[0]
    minima = []
    for i in sign_change:
        try:
            Rm = brentq(dv_theta_dR, Rgrid[i], Rgrid[i + 1], args=args)
        except ValueError:
            continue
        minima.append(Rm)
    return minima


def find_wall_radius(A, a, C6, a1, a2, eta, costheta, R_inner, R_outer):
    """Inner V=0 crossing (sigma-like radius) bracketed by [R_inner, R_outer]."""
    args = (A, a, C6, a1, a2, eta, costheta)
    f_in, f_out = v_theta(R_inner, *args), v_theta(R_outer, *args)
    if f_in * f_out > 0.0:
        return np.nan
    return brentq(v_theta, R_inner, R_outer, args=args)


# --------------------------------------------------------------------------- #
# Screening rules (return (ok: bool, reason: str|None))                       #
# --------------------------------------------------------------------------- #
def check_single_well(minima):
    """Rule 1: exactly one outer minimum on the grid."""
    if len(minima) == 0:
        return False, "no_outer_minimum"
    if len(minima) > 1:
        return False, "multiple_minima"
    return True, None


def check_no_collapse(A, a, C6, a1, a2, eta, costheta, Rm, cfg, Rgrid=None):
    """Rule 2: V(R,theta) monotonically decreasing on [R_min, Rm] with a repulsive wall.

    No barrier maximum inside the integration window (the -C6/R^6 term does not
    overtake the exponential before R_min), and V(R_min) > 0 so a real
    repulsive wall (V=0 crossing) exists -- checked at this angle.
    """
    args = (A, a, C6, a1, a2, eta, costheta)
    if v_theta(cfg["R_min"], *args) <= 0.0:
        return False, "no_repulsive_wall"
    if Rgrid is None:
        Rgrid = _radial_grid(cfg)
    inner = Rgrid[Rgrid <= Rm]
    if inner.size and np.any(dv_theta_dR(inner, *args) > 0.0):
        return False, "inner_collapse"
    return True, None


def check_angular_physicality(a1, a2, eta, cfg):
    """Rule 3: both radial-term angular factors stay above a positive floor for all theta.

    The dispersion factor (1+a1*P1+a2*P2) must stay positive so the tail remains
    attractive; the wall factor (1+eta*P1+a2*P2) must stay positive so the
    A*exp(-aR) term remains repulsive at every angle.
    """
    c = np.linspace(-1.0, 1.0, 201)
    if np.min(angular_factor(c, a1, a2)) <= cfg["g_floor"]:
        return False, "angular_inversion_tail"
    if np.min(angular_rep(c, eta, a2)) <= cfg["g_floor"]:
        return False, "angular_inversion_wall"
    return True, None


def check_well_depth(eps_by_angle, cfg):
    """Rule 4: well depth at every special angle within sane bounds."""
    for eps in eps_by_angle.values():
        if not (cfg["eps_min"] <= eps <= cfg["eps_max"]):
            return False, "well_depth_out_of_range"
    return True, None


# --------------------------------------------------------------------------- #
# Full characterization of a surviving candidate                              #
# --------------------------------------------------------------------------- #
def characterize(params, cfg):
    """Compute all descriptors at 0/90/180 plus anisotropy contrasts.

    Returns (descriptors: dict, reject_reason: str|None). With independent wall
    and tail anisotropies the well shape varies with angle, so the minimum
    position, wall radius, and depth are found separately on each angular cut.
    """
    A, a, C6, a1, a2, eta = (params[k]
                             for k in ("A", "a", "C6", "a1", "a2", "eta"))
    Rgrid = _radial_grid(cfg)

    ok, reason = check_angular_physicality(a1, a2, eta, cfg)
    if not ok:
        return None, reason

    desc, eps = {}, {}
    for ang in ANGLES:
        ct = COSTHETA[ang]
        minima = find_minima(A, a, C6, a1, a2, eta, ct, cfg, Rgrid)
        ok, reason = check_single_well(minima)
        if not ok:
            return None, f"{reason}@{int(ang)}"
        Rm = minima[0]

        ok, reason = check_no_collapse(A, a, C6, a1, a2, eta, ct, Rm, cfg, Rgrid)
        if not ok:
            return None, f"{reason}@{int(ang)}"

        Rwall = find_wall_radius(A, a, C6, a1, a2, eta, ct, cfg["R_min"], Rm)
        if not np.isfinite(Rwall):
            return None, f"no_repulsive_wall@{int(ang)}"

        eps[ang] = -v_theta(Rm, A, a, C6, a1, a2, eta, ct) * HARTREE2CM
        desc[f"Rm_{int(ang)}"] = Rm * BOHR2ANG
        desc[f"Rwall_{int(ang)}"] = Rwall * BOHR2ANG

    ok, reason = check_well_depth(eps, cfg)
    if not ok:
        return None, reason

    desc["eps_0"], desc["eps_90"], desc["eps_180"] = eps[0.0], eps[90.0], eps[180.0]
    desc["d_eps_odd"] = eps[0.0] - eps[180.0]
    desc["d_eps_even"] = 0.5 * (eps[0.0] + eps[180.0]) - eps[90.0]
    return desc, None


# --------------------------------------------------------------------------- #
# Sampling                                                                     #
# --------------------------------------------------------------------------- #
def _scale_log(u, lo, hi):
    return 10.0 ** (np.log10(lo) + u * (np.log10(hi) - np.log10(lo)))


def _scale_lin(u, lo, hi):
    return lo + u * (hi - lo)


def _scale_norm(u, sigma, mu=0.0):
    """Map a uniform Sobol coordinate to N(mu, sigma) via the inverse normal CDF."""
    u = np.clip(u, 1e-12, 1.0 - 1e-12)       # guard against ppf(0/1) = +-inf
    return mu + sigma * norm.ppf(u)


def sample_sobol(cfg):
    """Scaled Sobol design over the 8-D parameter cube -> DataFrame of params."""
    sampler = qmc.Sobol(d=8, scramble=True, seed=cfg["seed"])
    u = sampler.random_base2(cfg["sobol_m"])  # shape (2**m, 8)
    df = pd.DataFrame({
        "A":   _scale_log(u[:, 0], *cfg["A"]),
        "a":   _scale_log(u[:, 1], *cfg["a"]),
        "C6":  _scale_log(u[:, 2], *cfg["C6"]),
        "a1":  _scale_lin(u[:, 3], *cfg["a1"]),
        "a2":  _scale_lin(u[:, 4], *cfg["a2"]),
        "eta": _scale_norm(u[:, 5], cfg["eta_sigma"]),
        "B":   _scale_log(u[:, 6], *cfg["B"]),
        "M":   _scale_log(u[:, 7], *cfg["M"]),
    })
    return df


def reduced_masses(M, cfg):
    """Reduced masses with He and H2 partners (never sample mu directly)."""
    mu_He = cfg["m_He"] * M / (cfg["m_He"] + M)
    mu_H2 = cfg["m_H2"] * M / (cfg["m_H2"] + M)
    return mu_He, mu_H2


# --------------------------------------------------------------------------- #
# Generation loop                                                              #
# --------------------------------------------------------------------------- #
def generate(cfg):
    """Sample -> screen -> characterize. Returns (library_df, full_log_df)."""
    samples = sample_sobol(cfg)
    lib_rows, log_rows = [], []

    for idx, row in samples.iterrows():
        params = row.to_dict()
        mu_He, mu_H2 = reduced_masses(params["M"], cfg)

        base = {
            "id": idx,
            "A": params["A"], "a": params["a"], "C6": params["C6"],
            "a1": params["a1"], "a2": params["a2"], "eta": params["eta"],
            "B_cm1": params["B"], "M_amu": params["M"],
            "mu_He": mu_He, "mu_H2": mu_H2,
        }

        if cfg["couple_BM"]:
            # Reject implausibly heavy-but-stiff rotors (rough B ~ M^-1 trend).
            if params["B"] * params["M"] > 50.0:
                log_rows.append({**base, "accepted": False,
                                 "reject_reason": "implausible_BM"})
                continue

        desc, reason = characterize(params, cfg)
        if desc is None:
            log_rows.append({**base, "accepted": False, "reject_reason": reason})
            continue

        full = {**base, **desc}
        lib_rows.append(full)
        log_rows.append({**full, "accepted": True, "reject_reason": ""})

    lib_df = pd.DataFrame(lib_rows)
    log_df = pd.DataFrame(log_rows)
    return lib_df, log_df


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def main(cfg=CONFIG):
    lib_df, log_df = generate(cfg)
    n_total = len(log_df)
    n_acc = len(lib_df)

    lib_df.to_csv(cfg["library_csv"], index=False)
    log_df.to_csv(cfg["log_csv"], index=False)

    print(f"Sampled {n_total} candidates; accepted {n_acc} "
          f"({100.0 * n_acc / max(n_total, 1):.1f}%).")
    if n_total > n_acc:
        reasons = log_df.loc[~log_df["accepted"], "reject_reason"].value_counts()
        print("Rejection reasons:")
        for reason, count in reasons.items():
            print(f"  {reason:24s} {count}")
    print(f"Library -> {cfg['library_csv']}")
    print(f"Full log -> {cfg['log_csv']}")
    return lib_df, log_df


if __name__ == "__main__":
    main()
