"""
The forcing scenarios of the matrix, each a sequence of ``engine.Chunk`` of hourly catchment runoff volumes.

- ``era5``: the ERA5 runoff of the Columbia catchments, one chunk per month, from scripts/02_aggregate_era5.py.
- ``synthetic-burst`` and ``synthetic-square``: runoff designed to excite the artifacts of negative coefficients. A
  steady baseflow everywhere is interrupted either by a burst, one hour of intense runoff, or by a square wave, a day
  alternating that intensity with no runoff every hour, the highest frequency an hourly forcing can hold.
"""

from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import config
from .engine import Chunk

__all__ = [
    'era5_files',
    'era5_chunks',
    'catchment_areas',
    'synthetic_depths',
    'synthetic_chunks',
    'synthetic_baseflow_volumes',
    'SYNTHETIC',
    'SYNTHETIC_PERTURBATIONS',
]

# the synthetic forcings, depths of runoff per hour over every catchment. Each starts at steady state with the baseflow,
# is perturbed from hour 24, and runs long enough for the perturbation to drain from the longest path of the basin.
SYNTHETIC = {
    'start': '2001-01-01',
    'baseflow_mm_per_hour': 0.02,  # about 0.5 mm per day
    'intensity_mm_per_hour': 2.0,
    'first_hour': 24,
    'hours': 24 * 45,
}
SYNTHETIC_PERTURBATIONS = {
    'synthetic-burst': 1,  # one hour of intense runoff
    'synthetic-square': 24,  # a day of intense runoff alternating with none every hour, from intense to none
}

def era5_files(start: str, end: str) -> list[Path]:
    """The monthly catchment runoff files of the months from start to end, inclusive."""
    months = pd.period_range(start, end, freq='M')
    files = [config.ERA5_CATCHMENT_DIR / f'catchment_runoff_{m.year}{m.month:02d}.nc' for m in months]
    missing = [f for f in files if not f.exists()]
    if missing:
        raise FileNotFoundError(f'{len(missing)} catchment runoff files are missing, the first is {missing[0]}')
    if not files:
        raise ValueError(f'no months between {start} and {end}')
    return files


def era5_chunks(files: list[Path], river_ids: np.ndarray) -> Iterator[Chunk]:
    """One chunk per catchment runoff file, checked against the rivers of the original network."""
    if not files:
        raise ValueError('no runoff files given')
    for path in files:
        with xr.open_dataset(path) as ds:
            if not np.array_equal(ds['riverId'].to_numpy(), river_ids):
                raise ValueError(f'{path} does not list the rivers of the network in network order')
            if ds['catchment_runoff'].attrs.get('units') != 'm3':
                raise ValueError(f'{path} must hold runoff volumes')
            volumes = np.nan_to_num(ds['catchment_runoff'].to_numpy().astype(np.float32, copy=False))
            dates = ds['time'].to_numpy().astype('datetime64[s]')
        yield Chunk(dates, np.ascontiguousarray(volumes))


def catchment_areas() -> np.ndarray:
    """(n,) area of each Columbia catchment in m2, in network order, from the ERA5 weight table."""
    with xr.open_dataset(config.WEIGHTS_FILE) as ds:
        weights = ds[['riverId', 'area_sqm']].to_dataframe()
    river_ids = pd.read_parquet(config.NETWORK_FILE, columns=['riverId'])['riverId'].to_numpy()
    areas = weights.groupby('riverId')['area_sqm'].sum().reindex(river_ids).to_numpy(dtype=np.float64)
    if np.any(~np.isfinite(areas)) or np.any(areas <= 0):
        raise ValueError('every catchment must have a positive area')
    return areas


def synthetic_depths(scenario: str) -> np.ndarray:
    """(hours,) runoff depth in meters of each hour of a synthetic forcing, the same over every catchment."""
    if scenario not in SYNTHETIC_PERTURBATIONS:
        raise ValueError(f'{scenario} is not a synthetic scenario')
    depths = np.full(SYNTHETIC['hours'], SYNTHETIC['baseflow_mm_per_hour'] / 1000)
    perturbed = SYNTHETIC['first_hour'] + np.arange(SYNTHETIC_PERTURBATIONS[scenario])
    depths[perturbed] = np.where(perturbed % 2 == 0, SYNTHETIC['intensity_mm_per_hour'] / 1000, 0.0)
    if depths.min() < 0 or depths[SYNTHETIC['first_hour']] <= depths[0]:
        raise ValueError('the perturbation must start with intense runoff and never be negative')
    return depths


def synthetic_chunks(scenario: str, areas: np.ndarray, chunk_hours: int = 24 * 5) -> Iterator[Chunk]:
    """A synthetic forcing as runoff volumes, in chunks of chunk_hours."""
    hours = SYNTHETIC['hours']
    if chunk_hours <= 0 or hours % chunk_hours != 0:
        raise ValueError(f'chunk_hours must divide the {hours} hours of the synthetic forcing')
    depths = synthetic_depths(scenario)
    dates = np.datetime64(SYNTHETIC['start'], 's') + np.arange(hours) * np.timedelta64(config.DT_RUNOFF, 's')
    for start in range(0, hours, chunk_hours):
        stop = start + chunk_hours
        volumes = np.outer(areas, depths[start:stop]).astype(np.float32)
        yield Chunk(dates[start:stop], volumes)


def synthetic_baseflow_volumes(areas: np.ndarray) -> np.ndarray:
    """(n,) runoff volume per hour of the baseflow the synthetic forcing starts at steady state with."""
    if areas.ndim != 1:
        raise ValueError('areas must be one per catchment')
    volumes = areas * SYNTHETIC['baseflow_mm_per_hour'] / 1000
    if np.any(volumes <= 0):
        raise ValueError('the baseflow must be positive everywhere')
    return volumes
