"""
Check that stability.engine routes exactly as rr.Router does: one month of ERA5 catchment runoff over the Columbia,
routed both ways for each network type, must give identical discharge. Also reports the routing time of each.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/verify_engine.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import river_route as rr
import xarray as xr

from stability import config, engine, forcing, treatments

MONTH = '2000-01'
CHECK_DIR = config.DATA / 'verify_engine'


class KeepAll:
    """Keeps every chunk of discharge in memory."""

    def __init__(self) -> None:
        self.chunks = []

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        self.chunks.append(discharge.copy())


def route_with_router(network_type: str, dt: int, runoff_file: Path) -> np.ndarray:
    """The discharge rr.Router writes for one runoff file, at full precision."""
    out = CHECK_DIR / f'router_{network_type}_{dt}.nc'
    configs = rr.Configs(
        network_file=config.NETWORK_FILE, forcing='catchment', network_type=network_type, dt_routing=dt,
        runoff_files=[runoff_file], discharge_files=[out], unstable_coefficients='ignore', progress_bar=False,
        log_level='ERROR',
    )
    rr.Router(configs).set_discharge_writer(rr.router.writers.netcdf_writer).route()
    with xr.open_dataset(out) as ds:
        return ds['Q'].to_numpy()


def route_with_engine(network_type: str, dt: int, runoff_file: Path) -> tuple[np.ndarray, dict]:
    """The discharge stability.engine routes for one runoff file."""
    table = pd.read_parquet(config.NETWORK_FILE)
    treatment = 'standard' if network_type == 'standard' else 'stabilized'
    treated = treatments.treat(treatment, dt, table)
    keep = KeepAll()
    chunks = forcing.era5_chunks([runoff_file], table['riverId'].to_numpy())
    summary = engine.route(treated, dt, config.DT_RUNOFF, chunks, keep, threads=1)
    return np.concatenate(keep.chunks, axis=1), summary


if __name__ == '__main__':
    CHECK_DIR.mkdir(parents=True, exist_ok=True)
    month_file = forcing.era5_files(MONTH, MONTH)[0]
    for kind, step in (('standard', 3600), ('standard', 300), ('stabilized', 3600), ('stabilized', 300)):
        expected = route_with_router(kind, step, month_file)
        routed, run = route_with_engine(kind, step, month_file)
        identical = np.array_equal(expected, routed)
        print(f'{kind:10s} dt={step:5d}  identical={identical}  max|diff|={np.abs(expected - routed).max():.3g}  {run}')
        if not identical:
            raise SystemExit('the engine does not reproduce rr.Router')
