#!/usr/bin/env python3
"""Driver: produce a library of exactly 1000 accepted potentials via Sobol sampling.

Uses the model, ranges, and screening rules defined in potential_generator.py
(no assumptions changed). Generates the full 2**sobol_m Sobol design, screens it,
and keeps the first 1000 accepted candidates in Sobol order. The original Sobol
index is preserved as `id` for provenance and clean active-learning extension.
"""

import potential_generator as pg

N_TARGET = 1000


def main():
    cfg = dict(pg.CONFIG)  # default assumptions, sobol_m=11 -> 2048 candidates
    lib_df, log_df = pg.generate(cfg)

    n_acc = len(lib_df)
    if n_acc < N_TARGET:
        raise SystemExit(
            f"Only {n_acc} accepted from {len(log_df)} candidates; "
            f"increase CONFIG['sobol_m'] to reach {N_TARGET}.")

    lib_df = lib_df.iloc[:N_TARGET].reset_index(drop=True)

    lib_df.to_csv(cfg["library_csv"], index=False)
    log_df.to_csv(cfg["log_csv"], index=False)

    print(f"Sobol candidates: {len(log_df)}; total accepted: {n_acc}; "
          f"written to library: {len(lib_df)}.")
    print(f"Library -> {cfg['library_csv']}")
    print(f"Full log -> {cfg['log_csv']}")
    return lib_df, log_df


if __name__ == "__main__":
    main()
