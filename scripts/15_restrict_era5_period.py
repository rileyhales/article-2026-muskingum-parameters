"""
Cut the ERA5 cells routed over 2000-2019 to the matrix period, 2002-2011, so they need not be routed again. A
continuous run's 2002-2011 is the same simulation as one routed over 2002-2011 alone, except for the state it starts 2002
with: two years of history rather than three months of spin-up from a dry network. ``--check`` measures that difference
on the standard network at 1 h, the cheapest cell.

For each finished cell of data/results/era5-2000-2019 without a counterpart in data/results/era5-2002-2011 this
writes the counterpart: annual statistics of the years 2002-2011, the observation series of those hours, the event
windows (hard links, since every window lies inside the period), and meta.json noting the cell it was cut from.

Run with the project environment:
    uv run python scripts/15_restrict_era5_period.py [--check]
"""

import argparse
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from stability import config, engine, forcing, recorders, treatments

FIRST, END = np.datetime64(f'{config.ERA5_YEARS[0]}-01-01'), np.datetime64(f'{config.ERA5_YEARS[-1] + 1}-01-01')
FULL_FIRST = np.datetime64(f'{config.ERA5_FULL_YEARS[0]}-01-01')


def restrict_cell(source: Path, target: Path) -> None:
    """Write the 2002-2011 counterpart of one finished 2000-2019 cell."""
    if not (source / 'meta.json').exists():
        raise FileNotFoundError(f'{source} has not finished')
    target.mkdir(parents=True, exist_ok=True)
    for name in ('river_ids.npy', 'observation_ids.npy'):
        shutil.copy2(source / name, target / name)
    stats = np.load(source / 'annual_stats.npz')
    keep = np.isin(stats['years'], config.ERA5_YEARS)
    if keep.sum() != len(config.ERA5_YEARS):
        raise ValueError(f'{source} does not hold every year of the period')
    np.savez_compressed(target / 'annual_stats.npz', **{key: stats[key][keep] for key in stats.files})
    first = int((FIRST - FULL_FIRST) / np.timedelta64(1, 'h'))
    end = int((END - FULL_FIRST) / np.timedelta64(1, 'h'))
    observed = np.load(source / 'observations.npy', mmap_mode='r')
    np.save(target / 'observations.npy', np.ascontiguousarray(observed[:, first:end]))
    (target / 'events').mkdir(exist_ok=True)
    for window in sorted((source / 'events').glob('*.npy')):
        linked = target / 'events' / window.name
        if not linked.exists():
            os.link(window, linked)
    meta = json.loads((source / 'meta.json').read_text())
    meta.update({'scenario': config.ERA5_SCENARIO, 'cut_from': str(source.relative_to(config.ROOT)),
                 'routed': f'{config.ERA5_FULL_YEARS[0] - 1}-10 to {config.ERA5_FULL_YEARS[-1]}-12',
                 'kept': f'{config.ERA5_YEARS[0]} to {config.ERA5_YEARS[-1]}'})
    (target / 'meta.json').write_text(json.dumps(meta, indent=2))
    return


class KeepPeriod:
    """Hands a recorder only the chunks of the period, dropping the spin-up."""

    def __init__(self, recorder) -> None:
        self.recorder = recorder

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        keep = dates >= FIRST
        if keep.any():
            self.recorder.record(dates[keep], discharge[:, keep])


def check_spinup() -> pd.DataFrame:
    """The largest change in annual peak and volume between a 2002-2011 run and the cut of the 2000-2019 run."""
    table = pd.read_parquet(config.NETWORK_FILE)
    stats = recorders.AnnualStatsRecorder(len(table), config.ERA5_YEARS)
    files = forcing.era5_files(config.ERA5_SPINUP[0], f'{config.ERA5_YEARS[-1]}-12')
    chunks = forcing.era5_chunks(files, table['riverId'].to_numpy())
    engine.route(treatments.treat('standard', 3600, table), 3600, config.DT_RUNOFF, chunks, KeepPeriod(stats), 8)
    cut = np.load(config.run_dir(config.ERA5_SCENARIO, 'standard', 3600) / 'annual_stats.npz')
    rows = []
    for y, year in enumerate(config.ERA5_YEARS):
        peak, volume = cut['max'][y], cut['volume'][y]
        large = peak >= 0.1
        rows.append({
            'year': year,
            'peak_change_max': float(np.max(np.abs(stats.stats['max'][y] - peak)[large] / peak[large])),
            'volume_change_max': float(np.max(np.abs(stats.stats['volume'][y] - volume)[large] / volume[large])),
            'volume_change_basin': float(abs(stats.stats['volume'][y][-1] - volume[-1]) / volume[-1]),
        })
    return pd.DataFrame(rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--check', action='store_true', help='also measure the effect of the shorter spin-up')
    parser.add_argument('--cells', nargs='+', help='cut only these cells, such as substeps-xadj__dt0060s')
    args = parser.parse_args()
    full = config.RESULTS / config.ERA5_FULL_SCENARIO
    for cell in sorted(full.iterdir()):
        if args.cells and cell.name not in args.cells:
            continue
        target = config.RESULTS / config.ERA5_SCENARIO / cell.name
        if (cell / 'meta.json').exists() and not (target / 'meta.json').exists():
            restrict_cell(cell, target)
            print(f'{cell.name}: cut to {config.ERA5_YEARS[0]}-{config.ERA5_YEARS[-1]}', flush=True)
    if args.check:
        print(check_spinup().to_string(index=False))
