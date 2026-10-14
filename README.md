# 1D WHAM

Weighted Histogram Analysis Method for stitching umbrella-sampling histograms
along a single collective variable (CV). The umbrella bias is

```text
(1/2) k (cv - cv0)^2
```

## Folder structure

```text
wham1d/
├── Makefile
├── wham_1D.f90                   # main program
├── input.example                 # histograms + bootstrap errors
├── input.example.periodic        # 2π-periodic CV (iwrap=1)
├── input.example.timeseries      # build P(s) from time series
├── whaminput.example
├── whaminput.example.timeseries
├── README.md
└── test/                         # synthetic double-well example
    ├── input
    ├── whaminput
    ├── true_fes.dat
    └── PROB/
```

## Build

Requires `mpif90` (Open MPI / MPICH).

```bash
make          # builds ./wham1d
```

Or directly:

```bash
mpif90 wham_1D.f90 -o wham1d
```

On macOS with Homebrew gfortran, the Makefile sets `SDKROOT` automatically.

## Run

From a directory that contains `input`, `whaminput`, and the histogram / trajectory files:

```bash
mpirun -np 1 /path/to/wham1d
# or
mpirun -np N /path/to/wham1d
```

### Output: `free_energy.dat`

```text
# cv   F(cv)   dF(cv)   P(cv)   [kcal/mol]
```

| Column | Meaning |
|--------|---------|
| `cv` | Collective variable |
| `F(cv)` | Free energy (kcal/mol), shifted so \(\min F = 0\) |
| `dF(cv)` | Bootstrap std. error (kcal/mol); `0` if `nboot = 0` |
| `P(cv)` | Unbiased probability density |

`dF = 0` at the FE minimum by construction (alignment reference bin).

---

## Input file (`input`)

```text
ncv   umbr_n   T(K)   iwrap   [max_iter]   [itseries]   [nboot]
cvmin cvmax dcv
```

| Field | Meaning |
|-------|---------|
| `ncv` | Number of CVs. Must be **1**. |
| `umbr_n` | Number of umbrella windows. |
| `T(K)` | Temperature in Kelvin. |
| `iwrap` | `1` = 2pi-periodic CV (−pi...pi wrapping); `0` = no wrapping. |
| `max_iter` | Max WHAM SCF iterations (default **2000** if omitted). |
| `itseries` | `0` = read histograms (default); `1` = build histograms from time series. |
| `nboot` | Bootstrap samples for FE error (default **0** = no errors). |
| `cvmin` `cvmax` `dcv` | CV grid. `nbin = nint((cvmax-cvmin)/dcv)+1`. |

Uneven umbrella centers are allowed. If SCF does not reach `tolerance` by `max_iter`, a warning is printed and `free_energy.dat` is still written.

### Examples

Histograms, non-periodic CV, 50 bootstrap samples:

```text
1    7    300.0    0    2000    0    50
-2.00  2.00  0.05
```

Periodic (dihedral-like) CV:

```text
1    16   300.0    1    2000    0    0
-3.14  3.14  0.10
```

Time series -> histograms, with errors:

```text
1    7    300.0    0    2000    1    50
-2.00  2.00  0.05
```

---

## Umbrella list (`whaminput`)

Order matters: **file name first**, then umbrella parameters.

### Histogram mode (`itseries = 0`)

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
| `tolerance` | SCF threshold on \(k_BT \sum_i \|\Delta\ln a_i\|\). |
| `histogram_file` | Path to biased probability file. |
| `cv0` | Umbrella center (same units as the histogram CV). |
| `k` | Harmonic force constant in **kJ** (converted to kcal internally). |
| `nmax` | Number of samples used to build that histogram. |

### Time-series mode (`itseries = 1`)

```text
tolerance
timeseries_file_1
cv0_1   k_1(kJ)
timeseries_file_2
cv0_2   k_2(kJ)
...
```

No `nmax` — it is counted from the file. Each time-series file has:

```text
# optional comment lines
time   cv
...
```

Lines starting with `#` are ignored. The code bins `cv` onto `[cvmin,cvmax]` with spacing `dcv` and builds a normalized density \(P(s)\).

---

## Histogram files

For `itseries = 0`, each file is:

```text
cv   probability
```

Rows are matched to the WHAM grid by CV value
(`ibin = nint((cv-cvmin)/dcv)+1`). Points outside `[cvmin,cvmax]` are ignored.

Each histogram should be a **normalized density**:

```text
\Int P(s) ds  ≈  \sum_bins P_bin * dcv  =  1
```

If `|\Int P ds − 1| > 0.1`, the code warns and renormalizes
`P(s) <- P(s)/ [\Int P ds]`. This is required because WHAM weights windows as
`nmax * P_i(s)`.

---

## Periodicity (`iwrap`)

| `iwrap` | Behavior |
|---------|----------|
| `1` | Wrap umbrella centers and \((cv-\mathrm{cv0})\) with period \(2\pi\). |
| `0` | No wrapping (distances, unbound coordinates, etc.). |

---

## Free-energy errors (`nboot`)

With `nboot > 0`, errors use a **Poisson bootstrap** of histogram counts:

1. Bin counts \(N_i(s) = n_{\max,i}\,P_i(s)\,\Delta s\) (or raw counts from time series).
2. Draw \(N_i^*(s) \sim \mathrm{Poisson}(N_i(s))\).
3. Renormalize to \(P_i^*\), align each PMF to the reference bin (\(\min F\)).
4. \(\mathrm{d}F(s) =\) sample standard deviation over the `nboot` PMFs.

This assumes independent bin counts. For strongly correlated MD, a block bootstrap of trajectories is more rigorous (not implemented; full trajectories can be huge). Prefer decorrelated/strided frames and a sufficient `nboot`.


