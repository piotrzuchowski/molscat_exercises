#!/usr/bin/env python3
"""decoupled_v3 potential generator — V1 sampled by MAGNITUDE, not signed ratio.

Physics: under molecular-axis reversal theta -> pi-theta, P_lambda -> (-1)^lambda
P_lambda, so V1 -> -V1 while V0, V2 are unchanged. The SIGN of V1 is therefore a
gauge of the (arbitrary) theta=0 end, not a physical descriptor. We sample the
strength of the odd anisotropy, rho1 = |V1|/|V0|, and give V1 its own independent
radial form (exp wall + R^-7 dipole-dispersion tail):

    V(R,theta) = V0(R) + s1 * |V1|(R) * P1(cos t) + V2(R) * P2(cos t)

    V0(R) =      A e^{-aR}        -      C6 R^{-6}           (isotropic)
    V2(R) = zeta A e^{-aR}        - a2   C6 R^{-6}           (even, decoupled wall/tail)
    |V1|(R) = A1 e^{-b1 R} - C7 R^{-7}                        (odd MAGNITUDE)

|V1| is built from physical parameters so it has a well of depth eps1 = rho1*eps0
at R_{m,1}, with reduced steepness x1 = b1 R_{m,1} > 7:
    b1 = x1 / R_{m,1}
    A1 = 7 eps1 e^{x1} / (x1 - 7)
    C7 = eps1 x1 R_{m,1}^7 / (x1 - 7)
s1 = +-1 is an INTERNAL sign (geometry only); it is NOT an ML feature. P1=cos t,
P2 = 3cos^2 t - 1 (unnormalised, matching v2). Units: A,a,C6 atomic units.
ML descriptors are absolute: |d_eps_odd|, rel_odd = |d_eps_odd|/(eps_90+floor).
"""
import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import qmc, norm

HARTREE2CM = 219474.6313702
BOHR2ANG = 0.529177210903

CONFIG = {
    "sobol_m": 15,            # 2**15 = 32768 candidates (-> aim ~10k accepted)
    "seed": 2026,
    # base (even) parameters, atomic units
    "A":  (150.0, 5000.0),   # log-uniform repulsive prefactor [Hartree]
    "a":  (2.2,   2.9),      # log-uniform steepness [1/Bohr]
    "C6": (3.0,   40.0),     # log-uniform dispersion [au]
    "a2": (-0.4,  0.4),      # even tail anisotropy (linear)
    "zeta_sigma": 0.2,       # even wall anisotropy (Gaussian, mean 0)
    "B":  (1.0e-3, 10.0),    # log-uniform rotational constant [cm^-1]
    "M":  (8.0,   100.0),    # log-uniform mass [amu]
    # --- odd anisotropy (V1) magnitude sampling ---
    "rho1_edges": (1.0e-3, 0.05, 0.5, 1.0),   # mixture breakpoints
    "rho1_frac":  (0.40, 0.40, 0.20),          # weak / mid / strong fractions
    "m1": (0.90, 1.15),      # R_{m,1} = m1 * R_m(V0)  (V1 well near the V0 well)
    "x1": (9.0, 15.0),       # reduced steepness b1*Rm1 (linear; default range)
    # radial grid / screening (Bohr, cm^-1)
    "R_min": 4.0, "R_max": 30.0, "n_grid": 4000,
    "g_floor": 0.05,         # even angular factors must exceed this
    "eps_min": 5.0, "eps_max": 500.0,
    "eps_floor": 1.0e-6,     # cm^-1, for rel_odd denominator
    "m_He": 4.002602, "m_H2": 2.01588,
    "library_csv": "potentials_library.csv",
    "log_csv": "potentials_log.csv",
}

ANGLES = (0.0, 90.0, 180.0)
COSTHETA = {0.0: 1.0, 90.0: 0.0, 180.0: -1.0}


# --------------------------------------------------------------------------- #
# Core radial functions (R in Bohr, energy in Hartree)                        #
# --------------------------------------------------------------------------- #
def v0(R, A, a, C6):
    return A * np.exp(-a * R) - C6 * R ** -6


def dv0_dR(R, A, a, C6):
    return -a * A * np.exp(-a * R) + 6.0 * C6 * R ** -7


def v1_mag(R, A1, b1, C7):
    """Positive-magnitude odd component A1 e^{-b1 R} - C7 R^-7."""
    return A1 * np.exp(-b1 * R) - C7 * R ** -7


def v_theta(R, A, a, C6, a2, zeta, A1, b1, C7, s1, ct):
    """Full V(R,theta); ct = cos(theta)."""
    p1 = ct
    p2 = 3.0 * ct ** 2 - 1.0
    even = A * np.exp(-a * R) * (1.0 + zeta * p2) - C6 * R ** -6 * (1.0 + a2 * p2)
    odd = s1 * v1_mag(R, A1, b1, C7) * p1
    return even + odd


def dv_theta_dR(R, A, a, C6, a2, zeta, A1, b1, C7, s1, ct):
    p1 = ct
    p2 = 3.0 * ct ** 2 - 1.0
    deven = (-a * A * np.exp(-a * R) * (1.0 + zeta * p2)
             + 6.0 * C6 * R ** -7 * (1.0 + a2 * p2))
    dodd = s1 * (-b1 * A1 * np.exp(-b1 * R) + 7.0 * C7 * R ** -8) * p1
    return deven + dodd


# --------------------------------------------------------------------------- #
# Sampling helpers                                                            #
# --------------------------------------------------------------------------- #
def _log(u, lo, hi):
    return 10.0 ** (np.log10(lo) + u * (np.log10(hi) - np.log10(lo)))


def _lin(u, lo, hi):
    return lo + u * (hi - lo)


def _rho1(u, cfg):
    """Mixture: 40% log[1e-3,0.05], 40% log[0.05,0.5], 20% uniform[0.5,1.0]."""
    e0, e1, e2, e3 = cfg["rho1_edges"]
    f0, f1, _ = cfg["rho1_frac"]
    out = np.empty_like(u)
    m0 = u < f0
    m1 = (u >= f0) & (u < f0 + f1)
    m2 = u >= f0 + f1
    out[m0] = _log(u[m0] / f0, e0, e1)
    out[m1] = _log((u[m1] - f0) / f1, e1, e2)
    out[m2] = _lin((u[m2] - (f0 + f1)) / (1.0 - f0 - f1), e2, e3)
    return out


def sample_sobol(cfg):
    s = qmc.Sobol(d=10, scramble=True, seed=cfg["seed"])
    u = s.random_base2(cfg["sobol_m"])
    rng = np.random.default_rng(cfg["seed"])
    return pd.DataFrame({
        "A":    _log(u[:, 0], *cfg["A"]),
        "a":    _log(u[:, 1], *cfg["a"]),
        "C6":   _log(u[:, 2], *cfg["C6"]),
        "a2":   _lin(u[:, 3], *cfg["a2"]),
        "zeta": norm.ppf(np.clip(u[:, 4], 1e-9, 1 - 1e-9)) * cfg["zeta_sigma"],
        "B":    _log(u[:, 5], *cfg["B"]),
        "M":    _log(u[:, 6], *cfg["M"]),
        "rho1": _rho1(u[:, 7], cfg),
        "m1":   _lin(u[:, 8], *cfg["m1"]),
        "x1":   _lin(u[:, 9], *cfg["x1"]),
        "s1":   rng.choice([-1.0, 1.0], size=len(u)),
    })


# --------------------------------------------------------------------------- #
# V1 construction + characterization                                          #
# --------------------------------------------------------------------------- #
def v0_minimum(A, a, C6, cfg):
    """Return (Rm0 [Bohr], eps0 [Hartree, positive]) or (nan, nan)."""
    res = minimize_scalar(lambda R: v0(R, A, a, C6),
                          bounds=(2.5, cfg["R_max"]), method="bounded")
    Rm0 = float(res.x)
    V = v0(Rm0, A, a, C6)
    if V >= 0.0:
        return np.nan, np.nan
    return Rm0, -V


def build_v1(p, Rm0, eps0):
    """(A1, b1, C7) [au] for |V1| with depth eps1=rho1*eps0 at Rm1=m1*Rm0."""
    Rm1 = p["m1"] * Rm0
    eps1 = p["rho1"] * eps0
    x1 = p["x1"]
    b1 = x1 / Rm1
    A1 = 7.0 * eps1 * np.exp(x1) / (x1 - 7.0)
    C7 = eps1 * x1 * Rm1 ** 7 / (x1 - 7.0)
    return A1, b1, C7, Rm1, eps1


def _grid(cfg):
    return np.linspace(cfg["R_min"], cfg["R_max"], cfg["n_grid"])


def characterize(p, cfg):
    """Return (descriptor dict, reason|None)."""
    A, a, C6, a2, zeta = p["A"], p["a"], p["C6"], p["a2"], p["zeta"]
    # even angular factors must stay positive
    c = np.linspace(-1.0, 1.0, 201)
    p2 = 3.0 * c ** 2 - 1.0
    if np.min(1.0 + a2 * p2) <= cfg["g_floor"]:
        return None, "angular_inversion_tail"
    if np.min(1.0 + zeta * p2) <= cfg["g_floor"]:
        return None, "angular_inversion_wall"

    Rm0, eps0 = v0_minimum(A, a, C6, cfg)
    if not np.isfinite(Rm0):
        return None, "no_isotropic_well"
    A1, b1, C7, Rm1, eps1 = build_v1(p, Rm0, eps0)
    s1 = p["s1"]
    args = (A, a, C6, a2, zeta, A1, b1, C7, s1)

    Rgrid = _grid(cfg)
    desc, eps = {}, {}
    for ang in ANGLES:
        ct = COSTHETA[ang]
        aa = args + (ct,)
        # single outer minimum
        d = dv_theta_dR(Rgrid, *aa)
        sc = np.where((d[:-1] < 0.0) & (d[1:] >= 0.0))[0]
        mins = []
        for i in sc:
            try:
                mins.append(brentq(dv_theta_dR, Rgrid[i], Rgrid[i + 1], args=aa))
            except ValueError:
                pass
        if len(mins) == 0:
            return None, f"no_outer_minimum@{int(ang)}"
        if len(mins) > 1:
            return None, f"multiple_minima@{int(ang)}"
        Rm = mins[0]
        # repulsive wall + no inner collapse
        if v_theta(cfg["R_min"], *aa) <= 0.0:
            return None, f"no_repulsive_wall@{int(ang)}"
        inner = Rgrid[Rgrid <= Rm]
        if inner.size and np.any(dv_theta_dR(inner, *aa) > 0.0):
            return None, f"inner_collapse@{int(ang)}"
        # wall radius (V=0 crossing)
        f_in, f_out = v_theta(cfg["R_min"], *aa), v_theta(Rm, *aa)
        if f_in * f_out > 0.0:
            return None, f"no_repulsive_wall@{int(ang)}"
        Rwall = brentq(v_theta, cfg["R_min"], Rm, args=aa)
        eps[ang] = -v_theta(Rm, *aa) * HARTREE2CM
        desc[f"Rm_{int(ang)}"] = Rm * BOHR2ANG
        desc[f"Rwall_{int(ang)}"] = Rwall * BOHR2ANG

    for e in eps.values():
        if not (cfg["eps_min"] <= e <= cfg["eps_max"]):
            return None, "well_depth_out_of_range"

    for ang in ANGLES:
        desc[f"eps_{int(ang)}"] = eps[ang]
    d_odd = eps[0.0] - eps[180.0]
    desc["d_eps_odd_signed"] = d_odd                       # diagnostic only
    desc["abs_d_eps_odd"] = abs(d_odd)                     # ML descriptor
    desc["rel_odd"] = abs(d_odd) / (eps[90.0] + cfg["eps_floor"])
    desc["d_eps_even"] = 0.5 * (eps[0.0] + eps[180.0]) - eps[90.0]
    # V1 build parameters (needed by the deck builder)
    desc.update({"eps0_iso_cm1": eps0 * HARTREE2CM, "Rm0_iso_ang": Rm0 * BOHR2ANG,
                 "A1_au": A1, "b1_au": b1, "C7_au": C7,
                 "Rm1_ang": Rm1 * BOHR2ANG, "eps1_cm1": eps1 * HARTREE2CM})
    return desc, None


# --------------------------------------------------------------------------- #
def generate(cfg):
    samples = sample_sobol(cfg)
    lib, log = [], []
    for idx, row in samples.iterrows():
        p = row.to_dict()
        mu_He = cfg["m_He"] * p["M"] / (cfg["m_He"] + p["M"])
        mu_H2 = cfg["m_H2"] * p["M"] / (cfg["m_H2"] + p["M"])
        base = {"id": idx, "A": p["A"], "a": p["a"], "C6": p["C6"],
                "a2": p["a2"], "zeta": p["zeta"], "rho1": p["rho1"],
                "m1": p["m1"], "x1": p["x1"], "s1": p["s1"],
                "B_cm1": p["B"], "M_amu": p["M"],
                "mu_He": mu_He, "mu_H2": mu_H2}
        desc, reason = characterize(p, cfg)
        if desc is None:
            log.append({**base, "accepted": False, "reject_reason": reason})
            continue
        full = {**base, **desc}
        lib.append(full)
        log.append({**full, "accepted": True, "reject_reason": ""})
    return pd.DataFrame(lib), pd.DataFrame(log)


def main(cfg=CONFIG):
    lib, log = generate(cfg)
    lib.to_csv(cfg["library_csv"], index=False)
    log.to_csv(cfg["log_csv"], index=False)
    n_tot, n_acc = len(log), len(lib)
    print(f"Sampled {n_tot} candidates; accepted {n_acc} "
          f"({100.0 * n_acc / max(n_tot, 1):.1f}%).")
    if n_tot > n_acc:
        print("Rejections:")
        for r, c in log.loc[~log["accepted"], "reject_reason"].value_counts().items():
            print(f"  {r:28s} {c}")
    print(f"Library -> {cfg['library_csv']}  ({n_acc} rows)")
    return lib, log


if __name__ == "__main__":
    main()
