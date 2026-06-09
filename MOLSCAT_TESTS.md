# MOLSCAT scattering tests — run description

Atom + linear-rigid-rotor (MOLSCAT 2020, `ITYPE=1`) close-coupling cross
sections over a library of fictitious but physically screened potentials.
This document is self-contained so the batch can be moved to the **bem**
cluster (WCSS, PBS Pro) and run there.

---

## 1. What the calculation is

For each potential in `potentials_library.csv` we build one MOLSCAT input deck
and compute state-to-state integral cross sections (rotor J = 0,1,2,… ) at three
collision energies tied to the rotor constant B.

### Potential model

Unnormalized Legendre expansion (MOLSCAT convention), `P1 = cosθ`,
`P2 = 3cos²θ − 1`:

```
V(R,θ) = V0(R) + V1(R)·P1 + V2(R)·P2
V0(R)  = A·exp(−aR) − C6·R⁻⁶                (λ=0, isotropic)
V1(R)  = η·A·exp(−aR) − a1·C6·R⁻⁶          (λ=1, odd: wall and tail differ)
V2(R)  = a2·[A·exp(−aR) − C6·R⁻⁶]          (λ=2, quadrupolar)
```

Equivalently the wall and dispersion tail carry **independent** odd anisotropy:

```
V(R,θ) = A·exp(−aR)·(1 + η·P1 + a2·P2) − C6·R⁻⁶·(1 + a1·P1 + a2·P2)
```

Parameters `A, a, C6` are in atomic units; `a1, a2` linear in [−0.4, 0.4];
`η` Gaussian (mean 0, σ = 0.2); `B` in cm⁻¹; molecular mass `M` in 8–100 amu.
The library was Sobol-sampled and screened (single well, no inner collapse,
angular positivity, well depth 5–500 cm⁻¹) by `potential_generator.py`.

### Per-case scattering settings (`molscat_run.py` → `CONFIG`)

| Setting | Value | Meaning |
|---|---|---|
| `mu_col` | `mu_H2` | reduced mass uses the H₂ partner |
| `energy_factors` | (4, 16, 64) | collision energies = factor × B (cm⁻¹) |
| `RMAX_bohr` | 30.0 | outer propagation edge |
| `RMIN_bohr` | 4.0 | nominal inner edge (adaptively deepened, see below) |
| `RMIN_floor_bohr` | 2.5 | hard inner limit, never start inside this |
| `wall_safety` | 30 | require V0(RMIN) ≥ 30 × E_max |
| `DR` | 0.01 Å | log-derivative propagator step |
| `IPROPS` | 6 | diabatic modified log-derivative (Manolopoulos) |
| `JMAX` | 10 | rotor basis size |
| `ltop_buffer` | 12 | extra partial waves above semiclassical l_max |

### MOLSCAT deck specifics (why these flags)

- `IRMSET = 0` — disables MOLSCAT's inner turning-point search (FINDRM); the
  undamped −C6/R⁶ term collapses at small R, so we set RMIN ourselves.
- `BCYOMN = 1.0E6` — hard-wall log-derivative boundary. The λ=0 coupling carries
  a (−1)^J sign, so odd-J channels look "locally open" at RMIN and would trip the
  WKB initialiser; at RMIN the wavefunction is ≈ 0 (V0 ≫ E) so the hard wall is
  exact (verified: σ(0→0) matches a single-channel run).
- **Adaptive RMIN** (`choose_rmin`): RMIN is pushed inward until
  `V0(RMIN) ≥ wall_safety × E_max` (worst channel J=0, L=0), then the inward
  extension is re-screened for monotonic repulsion (no hole). A potential whose
  wall is too short fails *loudly* (`wall_too_short` / `inner_collapse`) rather
  than silently masking an open channel.

---

## 2. Files

| File | Role | Needed on cluster |
|---|---|---|
| `molscat_run.py` | driver: builds decks, runs MOLSCAT, parses σ | **yes** |
| `potentials_library.csv` | 1089 potentials + descriptors | **yes** |
| `potential_generator.py` | regenerates the library (not needed to run) | optional |
| `molscat_xsec.csv` | output: parsed cross sections | produced |
| `molscat_runs/` | output: per-case `.inp` decks and raw `.out` | produced |
| `plot_*.py` | local visualization only | no |

Only `molscat_run.py` + `potentials_library.csv` + the MOLSCAT binary are
required to run the batch.

---

## 3. Running it

### Driver invocation

```bash
python molscat_run.py            # validation batch: first 5 cases (CONFIG["n_cases"])
python molscat_run.py 1089       # all potentials
python molscat_run.py 200        # first 200
```

`n_cases` is the only positional argument. Outputs:
`molscat_xsec.csv` (one row per state-to-state cross section: id, B_cm1, mu,
energy, I, F, J_I, J_F, kind, sigma[Å²]) and decks/raw output under
`molscat_runs/`.

### One thing to change for the cluster

`CONFIG["binary"]` in `molscat_run.py` is currently the **local** path:

```python
"binary": "/Users/piotrzuchowski/develop/molscat/source_code/molscat-basic",
```

Point it at the MOLSCAT executable on bem (build it there, or copy it), e.g.

```python
"binary": "/home/pzuch007/molscat/source_code/molscat-basic",
```

The driver feeds each deck to the binary on stdin and reads stdout, so no other
MOLSCAT environment is required.

---

## 4. Moving to bem (WCSS, SLURM)

You already have an SSH alias (`Host bem`, user `pzuch007`, ProxyJump
`ui.wcss.pl`) and your own MOLSCAT job scripts there. From this directory:

```bash
# copy the inputs (driver + library + this doc)
rsync -avz molscat_run.py potentials_library.csv MOLSCAT_TESTS.md \
      bem:~/molscat_tests/

# log in
ssh bem
cd ~/molscat_tests
```

Two ways to run, depending on how your existing scripts are wired:

**(a) Let the Python driver call the MOLSCAT binary** (simplest). Point
`CONFIG["binary"]` at the bem executable (the same one your scripts launch) and
submit a single job:

```bash
#!/bin/bash
#SBATCH --job-name=molscat_lib
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=molscat_%j.out

module load python            # or: source ~/molscat_env/bin/activate
python molscat_run.py 1089
```

`sbatch run.slurm`; monitor with `squeue -u pzuch007`.

**(b) Reuse your existing per-deck MOLSCAT scripts.** The driver can just emit
the decks without running them — generate all `molscat_runs/pot_*.inp`, then feed
them to your own SLURM submitter. (Right now `molscat_run.py` runs each deck
inline; if you prefer "decks only", say so and I'll add a `--decks-only` flag.)

> MOLSCAT is single-threaded and processes one deck at a time, so a 1-core job
> per case is the natural unit. To spread the 1089 cases across the queue, use a
> SLURM **job array** (`#SBATCH --array=0-1088`) with one id per task. That
> needs the driver to accept an explicit id (or id range) instead of only a
> leading count — ask and I'll add `--ids` / `--start/--stop`.

### Python deps on bem

`molscat_run.py` needs `numpy`, `pandas`, `scipy`. Either load a cluster Python
module that provides them, or create a local environment once:

```bash
python -m venv ~/molscat_env
source ~/molscat_env/bin/activate
pip install numpy pandas scipy
```

### Bringing results back

```bash
rsync -avz bem:~/molscat_tests/molscat_xsec.csv ./
rsync -avz bem:~/molscat_tests/molscat_runs/   ./molscat_runs_cluster/
```

---

## 5. Reading the results

`molscat_xsec.csv` columns:

| column | meaning |
|---|---|
| `id` | potential id (matches `potentials_library.csv`) |
| `B_cm1` | rotor constant (cm⁻¹) |
| `mu` | reduced mass used (amu, H₂ partner) |
| `energy` | collision energy (cm⁻¹) = factor × B |
| `I`, `F` | initial / final pair-state index (1-based) |
| `J_I`, `J_F` | initial / final rotor J (= index − 1) |
| `kind` | `elastic` (I=F) or `inelastic` |
| `sigma` | integral cross section (Å²) |

A case that fails screening or MOLSCAT is reported on stdout with a reason
(`wall_too_short`, `inner_collapse`, `returncode=…`, `no_cross_sections_parsed`)
and simply omitted from the CSV.
