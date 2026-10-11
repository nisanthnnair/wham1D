# 1D WHAM

Weighted Histogram Analysis Method for stitching umbrella-sampling histograms
along a single collective variable (CV). The umbrella bias is

```text
(1/2) k (cv - cv0)^2
```

## Folder structure

```text
wham1d/
├── wham_1D.f90              # 1D histograms -> 1D free energy
├── input.example            # sample input (non-periodic CV, iwrap=0)
├── input.example.periodic   # sample input (2pi-periodic CV, iwrap=1)
├── whaminput.example        # sample umbrella list
├── README.md
```

## Build

Requires `mpif90` (Open MPI / MPICH). On macOS with Homebrew gfortran, the
Makefile sets `SDKROOT` automatically.

```bash
mpif90 wham_1D.f90  -o wham1d.x
```

## Run directory

From a working directory that contains the input files and histogram paths:

```bash
mpirun -np 1 /path/to/wham1d.x
# or
mpirun -np N /path/to/wham1d.x
```

Output file: `free_energy` with columns `cv  F(cv)  P(cv)` (kcal/mol).

## Required input files

### 1. `input`

```text
ncv   umbr_n   T(K)   iwrap
cvmin cvmax dcv
```

| Field   | Meaning |
|---------|---------|
| `ncv`   | Number of CVs. Must be `1` for `wham_1D`. |
| `umbr_n`| Number of umbrella windows. |
| `T(K)`  | Temperature in Kelvin. |
| `iwrap` | `1` = treat CV as 2π-periodic (−π…π wrapping); `0` = no wrapping (distances etc.). |
| `cvmin` | Lower edge of the CV grid. |
| `cvmax` | Upper edge of the CV grid. |
| `dcv`   | Bin width. Number of bins = `nint((cvmax-cvmin)/dcv)+1`. |

Examples:

```text
# non-periodic distance-like CV
1    7    300.0    0
0.00  10.00  0.05

# dihedral / 2pi-periodic CV
1    16   300.0    1
-3.14  3.14  0.10
```

### 2. `whaminput`

```text
tolerance
histogram_file_1
cv0_1   k_1(kJ)   nmax_1
histogram_file_2
cv0_2   k_2(kJ)   nmax_2
...
```

| Field | Meaning |
|-------|---------|
| `tolerance` | SCF convergence threshold on Σ\|Δln aᵢ\| × kT. |
| `histogram_file` | Path to the biased probability file for that window. |
| `cv0` | Umbrella center. |
| `k` | Harmonic force constant in **kJ** (converted internally to kcal). |
| `nmax` | Number of samples / weight for that window. |

There must be exactly `umbr_n` (file, parameters) pairs.

### 3. Histogram files (`PROB/...`)

One file per umbrella window. For `wham_1D` each line is:

```text
cv   probability
```

There must be exactly `nbin` lines, ordered from `cvmin` to `cvmax` with step `dcv`.

## Periodicity (`iwrap`)

When `iwrap = 1`:

- umbrella centers outside (−π, π] are wrapped into that range
- CV displacements `(cv − cv0)` used in the bias Boltzmann factor are wrapped by ±2π

When `iwrap = 0`, no wrapping is applied. Use this for non-periodic CVs such as
distances or unbounded coordinates.
