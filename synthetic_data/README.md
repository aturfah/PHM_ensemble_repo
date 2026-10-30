# Synthetic Data Examples 

The synthetic datasets used to benchmark PHM-ensemble and 7 baseline clustering methods include

    + Double ring
    + Crescent moons
    + Unequal-variance blobs
    + Two spirals
    + Compound
    + High-dim overlap

Another example, double ring of cliques, with multiresolution structure is used to test PHM-ensemble and Leiden algorithm.	  


## Layout

```
synthetic_data/
├── Makefile                            -- make check / make sim / make analysis / make clean
├── check_packages.R                    -- package-dependency check used by `make check`
├── simulation_scripts/                 -- seven data-generating scripts (01-06 + double_ring_of_cliques)
├── analysis_scripts/                   -- seven analysis scripts (01-06 + double_ring_of_cliques) + methods_lib.R
├── sim_data/                           -- created by `make sim`: all simulated datasets
├── analysis_results.txt                -- created by `make analysis`: the recorded ARI tables
└── double_ring_of_cliques_phm_result.pdf  -- created by `make analysis`: PHM's dendrogram + merge-height heatmap
```

`simulation_scripts/` and `analysis_scripts/` hold only code. Data lives
in `sim_data/`, which `make sim` creates at the top level; each analysis
script reaches across into `../sim_data/` for the one `.rds` file it
needs. `double_ring_of_cliques_analysis.R` additionally writes its PDF
one directory *up* from `analysis_scripts/` -- i.e. to this top level,
alongside the Makefile -- not into `analysis_scripts/` where the script
itself lives.


## Quickest path: the Makefile

From this top-level directory:

```bash
make check     # verify R + required packages are installed
make sim       # generate every simulated dataset into sim_data/
make analysis  # run every analysis script against sim_data/,
               # recording all output in analysis_results.txt
make clean     # remove sim_data/, analysis_results.txt, and
               # double_ring_of_cliques_phm_result.pdf
```

`make sim` creates `sim_data/` if it doesn't exist yet and runs each
`simulation_scripts/*_sim.R` script with `sim_data/` as its working
directory, so every `.rds`/`.txt`/`.png` it writes lands there.

`make analysis` requires `sim_data/` to exist (run `make sim` first if
it doesn't) and runs each `analysis_scripts/*_analysis.R` script in
turn, appending its full console output -- including the PHM ensemble
fit log and the final ARI table -- to `analysis_results.txt` at the top
level, under a `===== <script> =====` header. It keeps going even if
one script errors (the error is captured in the file and flagged on
screen) and exits non-zero at the end if anything failed, so a clean
exit status means every dataset was analyzed successfully.

The PHM_ensemble's `M` (replicates) and `rho` (subsample fraction) can
be overridden on the command line, and apply uniformly to all seven
scripts:

```bash
make analysis M=10 RHO=0.5
```

Defaults are `M=5, RHO=1.0`, matching `fit_phm_ensemble()`'s own
defaults -- also what you get running any single
`analysis_scripts/*_analysis.R` directly with no arguments.

`make clean` removes `sim_data/` and anything in it, any stray
simulated data files left directly in `simulation_scripts/`,
`analysis_results.txt`, and `double_ring_of_cliques_phm_result.pdf`.


## Checking dependencies

```bash
make check
```

verifies R itself and every CRAN/GitHub package the simulation and
analysis scripts depend on (`mclust`, `dbscan`, `igraph`, `FNN`, `RANN`,
`ggplot2`, `patchwork`, and the GitHub-only `PHM` package), and prints
install instructions for anything missing.

## simulation_scripts/

Seven self-contained R scripts, one per synthetic dataset used in the
comparison: six numbered 01-06, plus `double_ring_of_cliques_sim.R`
(deliberately unnumbered -- it isn't part of the original six). Each
script is independent -- no shared library files, no command-line
arguments, just `set.seed(...)` followed by the generator.

### Requirements

- R (any reasonably recent version; developed against R 4.x).
- Base R only, **except** `04_two_spirals_sim.R`, which needs the
  `RANN` package (used only for a diagnostic nearest-neighbor distance
  check, not for generating the data itself):
  ```r
  install.packages("RANN")
  ```
  (`make check`, above, covers this automatically.)

### Running manually

Each script writes its output into whatever directory it's run from.
`make sim` runs them with `sim_data/` as that directory (see above); to
do the same by hand, either run them from inside `sim_data/`:

```bash
mkdir -p sim_data
cd sim_data
Rscript ../simulation_scripts/01_double_rings_sim.R
Rscript ../simulation_scripts/02_crescent_moons_sim.R
Rscript ../simulation_scripts/03_unequal_variance_blobs_sim.R
Rscript ../simulation_scripts/04_two_spirals_sim.R
Rscript ../simulation_scripts/05_compound_sim.R
Rscript ../simulation_scripts/06_high_dim_overlap_sim.R
Rscript ../simulation_scripts/double_ring_of_cliques_sim.R
```

or run them from inside `simulation_scripts/` itself, which writes the
output beside the scripts instead:

```bash
cd simulation_scripts
Rscript 01_double_rings_sim.R
# ...etc
```

There's no dependency between scripts -- order doesn't matter, and each
can be re-run on its own to regenerate just that one dataset.

### Reproducibility

Every script's `set.seed(...)` call is the exact seed used for the data
already analyzed in this project. Re-running a script reproduces that
dataset's `X` and `y` (or, for double_ring_of_cliques, `X`, `y`, and
`ring`) bit-for-bit -- this was checked directly (not just assumed) for
the original six by re-running each script and diffing the result
against the originally-analyzed `.rds` files with `all.equal()`; all
six matched exactly. Do not change the seed if you want to reproduce
the analyzed data; change it only if you deliberately want a fresh draw
from the same generator.

### Output files

Each script writes two data files (plus a quick-look PNG of the ground
truth, for 01-06) into whichever directory it's run from (`sim_data/`
via `make sim`; see "Running manually" above for the alternative):

- `<name>.rds` -- an R list, `readRDS()`-able, with at least `X` (the
  point coordinates) and `y` (the integer ground-truth label, 0-indexed)
  plus a few dataset-specific parameters (centers, noise levels, etc.).
- `<name>.txt` -- the same data as a space-delimited plain-text table
  with a header row, for loading outside R (e.g.
  `pandas.read_csv(path, sep=" ")` in Python, or `readtable` in Julia).

| # | Script | Output stem | n | True K | Columns in the .txt file |
|---|--------|-------------|---|--------|---------------------------|
| 01 | `01_double_rings_sim.R` | `concentric_annuli_3000` | 3000 | 2 | `x1 x2 y` |
| 02 | `02_crescent_moons_sim.R` | `two_moons_3000` | 3000 | 2 | `x1 x2 y` |
| 03 | `03_unequal_variance_blobs_sim.R` | `unequal_variance_blobs_3000` | 3000 | 3 | `x1 x2 y` |
| 04 | `04_two_spirals_sim.R` | `two_spirals_3000` | 3000 | 2 | `x1 x2 y` |
| 05 | `05_compound_sim.R` | `compound_3000` | 3000 | 4 | `x1 x2 y` |
| 06 | `06_high_dim_overlap_sim.R` | `high_dim_overlap_3000` | 3000 | 4 | `x1 x2 ... x60 y` (d=60) |
| -- | `double_ring_of_cliques_sim.R` | `double_ring_of_cliques` | 3000 | 36 (fine) / 2 (coarse) | `x1 x2 blob ring` |

Notes specific to a few scripts:

- **06 (high-dim overlap)**: the data is 60-dimensional (`d=60`), so the
  `.txt` file has 61 columns (`x1`...`x60`, `y`). The PNG is a 2D PCA
  projection for visualization only -- the saved `X` in the `.rds`/`.txt`
  is the full 60-dim data.
- **double_ring_of_cliques**: unlike the other six, this dataset has
  *two* ground truths at once -- 36 individual cliques (`y`/`blob`: 18
  inner-ring + 18 outer-ring blobs) and the 2 coarse rings they belong
  to (`ring`). It's designed so both resolutions are simultaneously
  genuine: within-ring neighbors are always closer than the cross-ring
  radial partner. It has no PNG (no quick-look plot script was written
  for it).
- Every other script also writes a `<name>.png` ground-truth scatter
  plot (base R graphics, white background) purely as a visual sanity
  check -- it's not needed for downstream analysis.

## analysis_scripts/

Six analysis scripts (01-06, mirroring the simulation scripts) plus
`double_ring_of_cliques_analysis.R`, plus the shared library
`methods_lib.R` they all `source()`.

Each of the 01-06 scripts runs PHM (the ensemble fit,
`fit_phm_ensemble(X, M = ens_args$M, rho = ens_args$rho, label = label)`
-- `G` left at its own default and auto-extended as needed) alongside
seven baseline methods (k-means++, Ward and single-linkage hierarchical
clustering, mclust, HDBSCAN, Leiden, OPTICS) on its dataset and prints
the resulting ARI table.

`double_ring_of_cliques_analysis.R` is different in three ways: it
calls `fit_phm_ensemble()` too (same defaults, same `--M=`/`--rho=`
overrides), but only compares PHM against Leiden, not the full
7-baseline suite; it reports ARI at *both* levels of ground truth (fine
= 36 cliques, coarse = 2 rings) as two separate tables, since this is
the one example built to have two simultaneous resolutions of
structure; and it additionally writes
`../double_ring_of_cliques_phm_result.pdf` -- a combined dendrogram +
leaf-ordered merge-height heatmap of the PHM fit, with dashed lines
marking where the tree would be cut for each of the two resolutions.

Every analysis script accepts `--M=<value>` and `--rho=<value>`
command-line overrides (parsed by `parse_ensemble_args()` in
`methods_lib.R`), falling back to `M=5, rho=1.0` if run with no
arguments at all.

### Running manually

Each analysis script reads its dataset from `../sim_data/`, so it must
be run with `analysis_scripts/` as the working directory, and `sim_data/`
must already be populated (via `make sim`, or the manual steps above):

```bash
cd analysis_scripts
Rscript 01_double_rings_analysis.R
Rscript 02_crescent_moons_analysis.R
Rscript 03_unequal_variance_blobs_analysis.R
Rscript 04_two_spirals_analysis.R
Rscript 05_compound_analysis.R
Rscript 06_high_dim_overlap_analysis.R
Rscript double_ring_of_cliques_analysis.R
```

Add `--M=10 --rho=0.5` (or similar) to any of these to override the
ensemble parameters for that one run.

Each of 01-06 prints its ARI table to the console and writes nothing to
disk. `double_ring_of_cliques_analysis.R` also prints its two ARI
tables to the console, but additionally writes
`../double_ring_of_cliques_phm_result.pdf` (see above). `make analysis`
(see above) runs all seven this way automatically and
redirects/collects their console output into `analysis_results.txt`.
