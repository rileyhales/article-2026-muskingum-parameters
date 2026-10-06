"""
Run cells of the simulation matrix: a forcing scenario, a treatment, and a routing time step. Each cell writes to
``data/results/<scenario>/<treatment>__dt<NNNN>s``:

- ``meta.json``: the cell, the rivers its treatment changed, the reaches and subcycles routed, and the wall time
- ``river_ids.npy``: the riverId of each row of every array below (merge removes rows)
- ERA5 scenario: ``annual_stats.npz`` for every river, ``observations.npy`` (observation river, hour) for the
  observation rivers, and ``events/<name>.npy`` (river, hour) for every river inside each event window
- synthetic scenarios: ``series.npy`` (river, hour) for every river

Examples, with the project environment:
    uv run python scripts/04_run_matrix.py --scenario era5-2002-2011 --treatment reference --dt 60
    uv run python scripts/04_run_matrix.py --scenario synthetic-burst --all
"""

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import river_route as rr

from stability import config, engine, forcing, recorders, treatments

OBSERVATIONS_FILE = config.INPUTS / 'observation_rivers.parquet'
EVENTS_FILE = config.INPUTS / 'events.json'
SCENARIOS = (config.ERA5_SCENARIO, *forcing.SYNTHETIC_PERTURBATIONS)
SHORT_RIVER_TREATMENTS = ('subcycles', 'inflate-k', 'merge')
# with no reach too short, the stabilized network is the substeps network, so it is not run again
SAME_AS_SUBSTEPS_WITHOUT_SHORT = ('stabilized',)
# cells left out to save time: at 60 s the stabilized network differs from substeps only in the 4 reaches too short
CUT_CELLS = (('stabilized', 60),)


def matrix_cells() -> list[tuple[str, int]]:
    """
    Every (treatment, dt) of the matrix. At a dt with no river too short, a short-river treatment is the standard
    network and the stabilized network is the substeps network, so those cells are left out, as are CUT_CELLS.
    """
    network = rr.Network(config.NETWORK_FILE)
    cells = []
    for treatment in config.TREATMENTS:
        steps = (config.REFERENCE[1],) if treatment == 'reference' else config.DT_ROUTING
        for dt in steps:
            _, too_short = network.unstable_mask(float(dt))
            if treatment in (*SHORT_RIVER_TREATMENTS, *SAME_AS_SUBSTEPS_WITHOUT_SHORT) and not too_short.any():
                continue
            if (treatment, dt) in CUT_CELLS:
                continue
            cells.append((treatment, dt))
    if not cells:
        raise ValueError('the matrix must have cells')
    return cells


class AfterSpinup:
    """Hands a recorder only the chunks inside the study period, dropping the spin-up months."""

    def __init__(self, recorder, first_date: str) -> None:
        if recorder is None:
            raise ValueError('a recorder is required')
        self.recorder = recorder
        self.first = np.datetime64(first_date, 's')

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        keep = dates >= self.first
        if keep.any():
            self.recorder.record(dates[keep], discharge[:, keep])


def output_rows(treated: treatments.Treated, river_ids: np.ndarray) -> np.ndarray:
    """The row of routed discharge of each river id, -1 where the treatment removed the river."""
    routed = pd.Series(np.arange(treated.network.original_river_ids.shape[0]), index=treated.network.original_river_ids)
    rows = routed.reindex(river_ids).fillna(-1).to_numpy(dtype=np.int64)
    if rows.shape != river_ids.shape:
        raise ValueError('every river id must have a row or -1')
    return rows


def run_era5(treated: treatments.Treated, dt: int, out: Path, threads: int) -> dict:
    """Route the 20 years of ERA5 runoff after its spin-up, recording statistics, observations, and event windows."""
    table = pd.read_parquet(config.NETWORK_FILE, columns=['riverId'])
    n_out = treated.network.original_river_ids.shape[0]
    observed = pd.read_parquet(OBSERVATIONS_FILE)['riverId'].to_numpy()
    observed_rows = output_rows(treated, observed)
    observed_rows = observed_rows[observed_rows >= 0]
    np.save(out / 'observation_ids.npy', treated.network.original_river_ids[observed_rows])
    first, last = f'{config.ERA5_YEARS[0]}-01-01', f'{config.ERA5_YEARS[-1] + 1}-01-01'
    n_hours = int((np.datetime64(last) - np.datetime64(first)) / np.timedelta64(1, 'h'))
    stats = recorders.AnnualStatsRecorder(n_out, config.ERA5_YEARS)
    series = recorders.SeriesRecorder(out / 'observations.npy', n_out, n_hours, rows=observed_rows)
    keep = [stats, series]
    if EVENTS_FILE.exists():
        windows = {name: tuple(window) for name, window in json.loads(EVENTS_FILE.read_text()).items()}
        keep.append(recorders.WindowRecorder(out / 'events', windows, n_out))
    files = forcing.era5_files(config.ERA5_SPINUP[0], f'{config.ERA5_YEARS[-1]}-12')
    chunks = forcing.era5_chunks(files, table['riverId'].to_numpy())
    recorder = AfterSpinup(recorders.TeeRecorder(*keep), first)
    summary = engine.route(treated, dt, config.DT_RUNOFF, chunks, recorder, threads=threads)
    stats.save(out / 'annual_stats.npz')
    series.close()
    if len(keep) == 3:
        keep[2].close()
    return summary


def run_synthetic(scenario: str, treated: treatments.Treated, dt: int, out: Path, threads: int) -> dict:
    """Route a synthetic forcing from steady baseflow, recording every river's whole series."""
    areas = forcing.catchment_areas()
    n_out = treated.network.original_river_ids.shape[0]
    state = engine.steady_state(treated, dt, config.DT_RUNOFF, forcing.synthetic_baseflow_volumes(areas))
    series = recorders.SeriesRecorder(out / 'series.npy', n_out, forcing.SYNTHETIC['hours'])
    chunks = forcing.synthetic_chunks(scenario, areas)
    summary = engine.route(treated, dt, config.DT_RUNOFF, chunks, series, threads=threads, initial_state=state)
    series.close()
    if not np.isfinite(summary['final_state_sum']):
        raise ValueError('the synthetic run must stay finite')
    return summary


def run_cell(scenario: str, treatment: str, dt: int, threads: int, force: bool, root: Path = config.RESULTS) -> None:
    """Route one cell of the matrix, or of the comparison under ``root``, and write its outputs unless they exist."""
    out = config.run_dir(scenario, treatment, dt, root)
    if (out / 'meta.json').exists() and not force:
        print(f'{out.relative_to(config.ROOT)} exists, skipping', flush=True)
        return
    if scenario not in SCENARIOS:
        raise ValueError(f'unknown scenario {scenario}')
    out.mkdir(parents=True, exist_ok=True)
    table = pd.read_parquet(config.NETWORK_FILE)
    too_long, too_short = rr.Network(table).unstable_mask(float(dt))
    treated = treatments.treat(treatment, dt, table)
    np.save(out / 'river_ids.npy', treated.network.original_river_ids)
    if scenario == config.ERA5_SCENARIO:
        summary = run_era5(treated, dt, out, threads)
    else:
        summary = run_synthetic(scenario, treated, dt, out, threads)
    meta = {
        'scenario': scenario, 'treatment': treatment, 'dt_routing': dt, 'dt_runoff': config.DT_RUNOFF,
        'description': (config.TREATMENTS | config.CHANGE_TREATMENTS)[treatment], 'rivers_too_long': int(too_long.sum()),
        'rivers_too_short': int(too_short.sum()), 'rivers_removed': int((~treated.kept).sum()), **summary,
        'threads': threads, 'river_route_version': rr.__version__, 'river_route_path': str(Path(rr.__file__).parent),
        'finished': datetime.now(UTC).isoformat(timespec='seconds'),
    }
    (out / 'meta.json').write_text(json.dumps(meta, indent=2))
    print(f'{out.relative_to(config.ROOT)}  {summary["wall_seconds"]:.1f} s', flush=True)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--scenario', required=True, choices=SCENARIOS)
    parser.add_argument('--treatment', choices=list(config.TREATMENTS))
    parser.add_argument('--dt', type=int)
    parser.add_argument('--all', action='store_true', help='run every cell of the matrix not yet run')
    parser.add_argument('--part', type=str, default='0/1', help='i/n: run every n-th cell starting at i')
    parser.add_argument('--cells', nargs='+', help='treatment:dt cells to run, such as stabilized:3600')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--force', action='store_true', help='rerun a cell whose outputs exist')
    args = parser.parse_args()
    if args.all:
        part, parts = (int(v) for v in args.part.split('/'))
        for cell_treatment, cell_dt in matrix_cells()[part::parts]:
            run_cell(args.scenario, cell_treatment, cell_dt, args.threads, args.force)
    elif args.cells:
        for cell in args.cells:
            cell_treatment, cell_dt = cell.split(':')
            run_cell(args.scenario, cell_treatment, int(cell_dt), args.threads, args.force)
    elif args.treatment and args.dt:
        run_cell(args.scenario, args.treatment, args.dt, args.threads, args.force)
    else:
        parser.error('give --treatment and --dt, --cells, or --all')
