# Negative Muskingum coefficients in large hydrofabrics

Experiments and manuscript draft on negative Muskingum coefficients in large hydrofabrics: how many reaches of the
GEOGLOWS River Forecast System v3 (TDX-Hydro) hydrofabric have them at routing time steps from 30 s to 60 min, what
artifacts they produce, how far downstream those propagate, and how the treatments compare. Routing is river-route v3
from the fork in `../river-route`, on the Columbia River basin.

## Layout

| Folder | What it holds |
|---|---|
| `paper/` | the manuscript and its LaTeX build (see `paper/README.md`) |
| `scripts/` | the Python that prepares the inputs, runs the simulations, and makes the figures and tables; `scripts/stability/` is the package they share |
| `hydrography/` | the hydrofabrics of the census as distributed: HydroRIVERS v1.0 and MERIT-Basins (MERIT Hydro v0.7, Basins v1, bugfix 1) |
| `data/inputs/` | the prepared Columbia network, ERA5 catchment runoff, observation reaches, and event windows |
| `data/results/` | the simulation matrix, one folder per cell |
| `data/river-route-changes/` | the comparison of river-route before and after this study (step 17), one folder per version |
| `vendor/river-route-matrix/` | the river-route that routed the matrix, copied unchanged (not tracked) |
| `data/logs/` | the printed output of runs and analyses |
| `data/verify_engine/` | the river-route Router outputs that `verify_engine.py` compares against |
| `figures/` | the figures, as PNG and PDF |
| `tables/` | the summary tables, as CSV |

## Environment

Every script runs from the project folder with the river-route environment, which has every dependency:

```bash
../river-route/.venv/bin/python scripts/<script>.py
```

Two versions of river-route are used. The simulation matrix (steps 4 to 15) was routed with the fork at commit 0ac029a
carrying one change: negative discharge written as routed instead of clamped to zero. That code is copied, unchanged,
into `vendor/river-route-matrix` (not tracked; `COMMIT` and `uncommitted.diff` record its source), and a script runs
with it when that folder is on `PYTHONPATH`:

```bash
PYTHONPATH=vendor/river-route-matrix ../river-route/.venv/bin/python scripts/<script>.py
```

`../river-route` itself now carries the changes this study found it needs (uncommitted): reaches too long for the step
divided into the fewest substeps with x adjusted by Eq. 16, discharge reported as the trapezoidal mean of the
computed discharges (Eq. 9), negative discharge written as routed, and networks with a reach of zero length, an x above 1/2, or a reach
draining into itself refused. Step 17 compares the two versions. `scripts/verify_engine.py` checks the harness against
`rr.Router` under either version, and under the matrix version also checks the matrix's stabilized treatment. After
editing a river-route kernel, delete its numba cache (`river_route/**/*.nbi`, `*.nbc`): `route_job` inlines the routing
methods, but its cache is keyed only on `router/_numba_kernels.py`.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_prepare_columbia.py` | Subsets the RFS v3 region 7020014250 to the Columbia (29,405 reaches): network, metadata, ERA5 weight table |
| 2 | `scripts/02_aggregate_era5.py` | Aggregates hourly ERA5 runoff, Oct 1999 to Dec 2019, to catchment volumes once (`data/inputs/era5_catchment_runoff`) |
| 3 | `scripts/03_select_observations.py` | Chooses 526 observation reaches every 20 km along 14 tributary main stems and the Columbia |
| — | `scripts/verify_engine.py` | Checks that `stability/engine.py` reproduces `rr.Router` bit for bit |
| 4 | `scripts/04_run_matrix.py` | Runs cells of the simulation matrix (below) |
| 5 | `scripts/05_census.py` | Census of coefficient signs, short-reach topology, and stabilization cost, Columbia and global |
| 6 | `scripts/06_theory_figures.py` | Figures of the theory section |
| 7 | `scripts/07_select_events.py` | Chooses the four ERA5 event windows from the reference run (`data/inputs/events.json`) |
| 8 | `scripts/08_analyze_synthetic.py` | Errors of every synthetic cell against the reference |
| 9 | `scripts/09_single_defect.py` | Isolates the downstream influence of single reaches too short for Δt |
| 10 | `scripts/10_analyze_era5.py` | Errors of every ERA5 cell against the reference: annual statistics, events, stem profiles |
| 11 | `scripts/11_era5_negative_flow_example.py` | Reaches with negative discharge under ERA5 at 1 h, and the worst one's hydrograph |
| 12 | `scripts/12_benchmark_cost.py` | Routing cost of every cell, timed alone on in-memory forcing, and accuracy against cost |
| 13 | `scripts/13_study_area_map.py` | Map of the study main stems, observation reaches, and coefficient signs at 1 h |
| 14 | `scripts/14_treatment_tables.py` | Tables of every treatment and step against peak, peak time, and volume; refreshes the manuscript |
| 15 | `scripts/15_restrict_era5_period.py` | Cuts finished 2000–2019 ERA5 cells to 2002–2011; `--check` measures the shorter spin-up |
| 16 | `scripts/16_hydrofabric_census.py` | Census of every reach of TDX-Hydro as delineated (~/data/TDXHydroGeoParquet), RFS v3, the Columbia, HydroRIVERS, and MERIT-Basins, with k from Eq. 17 and x = 0.2, and the census figures |
| 17 | `scripts/17_river_route_changes.py` | river-route before and after this study, each against its own reference: reporting convention (standard network), its stabilized network, and the substep count with x adjusted, under the burst, the square wave, and ERA5 2002–2011, in `data/river-route-changes/<before or after>/` |

Every run times itself as it goes (`scripts/stability/engine.py`, `Timer`): `meta.json` records the seconds spent preparing forcing,
routing with river-route's kernels, and recording outputs, the routing seconds of every chunk with its simulated hours,
and the routing seconds per simulated year, overall and steady (without the first chunk, which also loads the
kernels), so the cost of a cell is known from any part of a run. Runs that share the machine, as the two matrix workers
do, time each other's contention; `12_benchmark_cost.py` measures the same quantity for every cell alone on one month.

The ERA5 reference must run before step 7, and step 7 before the other ERA5 cells, which keep the event windows.

## The simulation matrix

Each cell is a forcing scenario, a treatment, and a routing time step, and writes to
`data/results/<scenario>/<treatment>__dt<NNNN>s/`.

- Scenarios: `era5-2002-2011` (routed from October 2001, three months of spin-up discarded), `synthetic-burst`,
  `synthetic-square`
- Time steps: 30, 60, 300, 600, 900, 1800, 3600 s
- Treatments:

| Treatment | Network type | What it does to reaches outside the window 2kx ≤ Δt ≤ 2k(1−x) |
|---|---|---|
| `standard` | standard | nothing: every reach routed whole at Δt |
| `substeps` | stabilized | reaches too long split into ⌈2kx/Δt⌉ equal sub-reaches, x unchanged |
| `substeps-xadj` | stabilized | as substeps, with x of every sub-reach lowered to 1/2 − N(1/2 − x) to keep the reach's diffusion |
| `subcycles` | stabilized | reaches too short routed in ⌈Δt/(2k(1−x))⌉ steps of their own |
| `stabilized` | stabilized | substeps and subcycles together: river-route `network_type='stabilized'` |
| `inflate-k` | standard | k of reaches too short raised to Δt/(2(1−x)) |
| `merge` | standard | reaches too short deleted; their upstreams and runoff joined to the reach downstream |
| `reference` | stabilized | Δt = 30 s with every reach's own step at most k/5: the time-converged solution |

At Δt = 30 s no reach is too short, so `subcycles`, `inflate-k`, and `merge` are identical to `standard` and are not
run, nor is `stabilized`, which is `substeps` there. Nor is `stabilized` at 60 s, where it differs from `substeps` in
four subcycled reaches (`CUT_CELLS` in `04_run_matrix.py`). Each scenario therefore has 44 cells and the reference.

Cells routed over ERA5 2000–2019 (`data/results/era5-2000-2019`) can be cut to 2002–2011 with
`15_restrict_era5_period.py`, which gives the same values as a run from October 2001 in every year kept; a cut cell's
`meta.json` names the cell it was cut from.

### Outputs of a cell

- `meta.json`: the cell, the reaches its treatment changed, reaches and subcycles routed, wall time
- `river_ids.npy`: the riverId of each row of the arrays below (merge removes rows)
- ERA5: `annual_stats.npz` (every reach and year: peak and hour, minimum, volume, negative hours and volume,
  oscillation hours, total variation), `observations.npy` with `observation_ids.npy` (observation reaches, every
  hour of 2002–2011), and `events/<name>.npy` (every reach, every hour of each event window)
- synthetic: `series.npy` (every reach, every hour)

All arrays are float32 at full precision, as routed. River-route's default zarr writer, which rounds to 12 mantissa
bits, is never used.

## Package

`scripts/stability/` holds the shared code: `config` (paths, matrix), `hydrofabric` (network table, main stems, distances),
`treatments` (the network of each treatment), `engine` (the routing loop over river-route's kernels), `forcing`
(ERA5 and synthetic chunks), `recorders`, `metrics`, `theory` (analytic properties of the scheme), and `plotting`.
