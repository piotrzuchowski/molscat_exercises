import numpy as np
"""
helper functions for energy grid and k(T) integration 

"""
KB_CM = 0.6950348        # cm^-1 / K
CM2J  = 1.98644586e-23   # 1 cm^-1 -> J (hc)
KB_J  = 1.380649e-23     # J / K
AMU   = 1.66053907e-27   # kg
ANG2  = 1.0e-20          # Ang^2 -> m^2

def energy_grid(B, jmax, T_hi=50.0, E_min=0.05, E_split=50.0, n_log=40, n_lin=12):
    """Absolute collision-energy grid [cm^-1] for k(T) up to T_hi (E_max = 6 kT_hi)."""
    E_max = 6.0 * KB_CM * T_hi                       # ~208.5 cm^-1 at 50 K
    E = np.concatenate([np.geomspace(E_min, E_split, n_log),
                        np.linspace(E_split, E_max, n_lin + 1)[1:]])
    thr = np.array([B * j * (j + 1) for j in range(1, jmax + 1)])
    thr = thr[(thr > E_min) & (thr < E_max)]          # rotor openings in range
    E = np.concatenate([E, 0.98 * thr, thr, 1.02 * thr])
    return np.unique(E)

def k_of_T(E_cm, sigma_ang2, T, mu_amu):
    """Boltzmann-averaged rate k(T) [cm^3/s].  E in cm^-1, sigma in Ang^2, T in K."""
    E, s = E_cm * CM2J, sigma_ang2 * ANG2
    kT, mu = KB_J * T, mu_amu * AMU
    integ = np.trapz(s * E * np.exp(-E / kT), E)
    return np.sqrt(8.0 / (np.pi * mu)) * kT**-1.5 * integ * 1.0e6   # m^3/s -> cm^3/s