"""
Aggregate the hourly ERA5 runoff of 2000 through 2019, and the months of spin-up before it, onto the Columbia
catchments once, as monthly catchment runoff volume files in the schema river-route's CatchmentRunoff reads, so every
simulation of the matrix reads the same forcing.

Run with the project environment:  uv run python scripts/02_aggregate_era5.py
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import river_route as rr
import xarray as xr

from stability import config

THREADS = 8
ENCODING = {'catchment_runoff': {'zlib': True, 'complevel': 1, 'shuffle': True}}


def era5_files() -> list[Path]:
    """Every monthly ERA5 file of the study period, in time order."""
    months = pd.period_range(config.ERA5_SPINUP[0], f'{config.ERA5_YEARS[-1]}-12', freq='M')
    files = [config.ERA5_ROOT / f'year={m.year}' / f'era5_{m.year}{m.month:02d}.nc' for m in months]
    missing = [f for f in files if not f.exists()]
    if missing:
        raise FileNotFoundError(f'{len(missing)} ERA5 files are missing, the first is {missing[0]}')
    if len(files) < 12 * len(config.ERA5_YEARS):
        raise ValueError('expected one ERA5 file per month of the study period')
    return files


def aggregate_month(grid: rr.GridRunoff, era5_file: Path, pool: ThreadPoolExecutor) -> Path:
    """Write one month of catchment runoff volumes and return the path written."""
    out = config.ERA5_CATCHMENT_DIR / f'catchment_runoff_{era5_file.stem.split("_")[1]}.nc'
    if out.exists():
        return out
    ds = grid.to_dataset(era5_file, thread_pool=pool, threads=THREADS)
    steps = np.diff(ds['time'].to_numpy()).astype('timedelta64[s]').astype(int)
    if not np.all(steps == config.DT_RUNOFF):
        raise ValueError(f'{era5_file} is not hourly')
    if ds['catchment_runoff'].attrs.get('units') != 'm3':
        raise ValueError('catchment runoff must be written as volumes')
    ds.to_netcdf(out.with_suffix('.partial'), encoding=ENCODING)
    out.with_suffix('.partial').rename(out)
    return out


if __name__ == '__main__':
    config.ERA5_CATCHMENT_DIR.mkdir(parents=True, exist_ok=True)
    runoff = rr.GridRunoff(
        config.WEIGHTS_FILE, var_grid_runoff='ro', var_x='longitude', var_y='latitude', var_t='valid_time',
        as_volumes=True,
    )
    network_ids = rr.Network(config.NETWORK_FILE).river_ids
    if not np.array_equal(runoff.river_ids, network_ids):
        raise ValueError('the weight table must list the rivers of the network in network order')
    with xr.open_dataset(era5_files()[0]) as first:
        if first['ro'].attrs.get('units') != 'm':
            raise ValueError('ERA5 runoff is expected in meters')
    with ThreadPoolExecutor(THREADS) as thread_pool:
        for month_file in era5_files():
            written = aggregate_month(runoff, month_file, thread_pool)
            print(written.name, flush=True)
