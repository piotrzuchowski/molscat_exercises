import numpy as np

def set_jmax(B, E, cap=30):
    """
    Set the rotational basis Jmax from the rotational constant.

    base = sqrt(E / B)        # rough J reachable at energy E
    Jmax = base + 5           # margin above the populated range
    capped at `cap` (default 30).

    Parameters
    ----------
    B : float   rotational constant (same energy units as E, e.g. cm^-1)
    E : float   reference energy (e.g. kT, or a cutoff energy)
    cap : int   hard upper limit on Jmax

    Returns
    -------
    int : chosen Jmax
    """
    base = np.sqrt(E / B)
    jmax = int(round(base)) + 5
    return min(jmax, cap)