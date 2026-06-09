#!/usr/bin/env python3
"""Parse MOLSCAT state-to-state integral cross sections into inelastic tables.

Reads one or more molscat-basic .out files (ISIGPR=2 output) and, for each,
extracts SIG(F,I) in Angstrom^2 at every energy. Builds, per initial level,
the inelastic cross section sigma_inel(I) = sum_{F!=I} SIG(F,I), and writes:
  * <stem>_xs.csv         : long table E, F, I, j_F, j_I, SIG, inelastic flag
  * <stem>_inel.csv       : wide table E vs sigma_inel(from each open initial j)
and prints the inelastic-from-ground table to stdout.

Usage: parse_xs.py file1.out [file2.out ...]
       parse_xs.py            # all *.out in cwd
"""

import glob
import re
import sys
import numpy as np
import pandas as pd

# data line: ENERGY JTOTL JSTEP JTOTU F I SIG   (E and SIG in E-notation, 5 ints)
NUM = r"[-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?"
LINE = re.compile(
    rf"^\s*({NUM})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+({NUM})\s*$")


def parse_file(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            m = LINE.match(ln)
            if not m:
                continue
            E, jtl, jst, jtu, F, I, sig = m.groups()
            rows.append((float(E.replace("D", "E")), int(F), int(I),
                         float(sig.replace("D", "E"))))
    df = pd.DataFrame(rows, columns=["E", "F", "I", "SIG"])
    if df.empty:
        return df
    # level index (1-based) -> j  (JLEVEL = 0,1,2,... in energy order)
    df["j_F"] = df["F"] - 1
    df["j_I"] = df["I"] - 1
    df["inelastic"] = df["F"] != df["I"]
    return df


def summarize(df):
    """Wide table: rows=energy, cols=sigma_inel from each initial level."""
    inel = df[df["inelastic"]].groupby(["E", "I"])["SIG"].sum().reset_index()
    wide = inel.pivot(index="E", columns="I", values="SIG").fillna(0.0)
    wide.columns = [f"inel_from_j{c-1}" for c in wide.columns]
    return wide.sort_index()


def main(argv):
    files = argv[1:] or sorted(glob.glob("*.out"))
    if not files:
        print("no .out files given/found")
        return
    for path in files:
        stem = path.rsplit(".", 1)[0]
        df = parse_file(path)
        if df.empty:
            print(f"\n=== {path}: no cross-section data found ===")
            continue
        df.to_csv(f"{stem}_xs.csv", index=False)
        wide = summarize(df)
        wide.to_csv(f"{stem}_inel.csv")

        print(f"\n=== {path} ===")
        print(f"energies: {df['E'].nunique()}   "
              f"levels: {df['I'].nunique()}   "
              f"-> {stem}_xs.csv, {stem}_inel.csv")
        # show inelastic-from-ground (j=0) where nonzero
        g = df[(df["I"] == 1) & df["inelastic"]]
        tot = g.groupby("E")["SIG"].sum()
        tot = tot[tot > 0]
        if tot.empty:
            print("  (no open inelastic channels from j=0 in this grid)")
            continue
        print("  Inelastic cross sections from ground state j=0 (Angstrom^2):")
        # break down by final j
        piv = (g.pivot_table(index="E", columns="j_F", values="SIG",
                             aggfunc="sum").fillna(0.0))
        piv.columns = [f"0->{c}" for c in piv.columns]
        piv["total_inel"] = piv.sum(axis=1)
        piv = piv[piv["total_inel"] > 0]
        with pd.option_context("display.float_format", lambda v: f"{v:.4g}",
                               "display.width", 200,
                               "display.max_columns", 30):
            print(piv.to_string())


if __name__ == "__main__":
    main(sys.argv)
