"""
Measures of a routed hydrograph and of its departure from the reference, vectorized over rivers: series are
(river, hour) arrays of hourly mean discharge.
"""

import json
from pathlib import Path

import numpy as np

from . import config

__all__ = ['Cell', 'load_cell', 'align', 'peaks', 'nyquist_amplitude', 'dip_below_start', 'oscillation_hours']


class Cell:
    """The outputs of one matrix cell: its meta, river ids, and the path of each array."""

    def __init__(self, scenario: str, treatment: str, dt: int) -> None:
        self.directory = config.run_dir(scenario, treatment, dt)
        if not (self.directory / 'meta.json').exists():
            raise FileNotFoundError(f'{self.directory} has not been run')
        self.meta = json.loads((self.directory / 'meta.json').read_text())
        self.river_ids = np.load(self.directory / 'river_ids.npy')
        self.treatment, self.dt, self.scenario = treatment, dt, scenario
        if self.meta['dt_routing'] != dt:
            raise ValueError(f'{self.directory} holds dt {self.meta["dt_routing"]}')

    def array(self, name: str, mmap: bool = True) -> np.ndarray:
        """A .npy array of the cell by its name relative to the cell directory, without the extension."""
        path = self.directory / f'{name}.npy'
        if not path.exists():
            raise FileNotFoundError(path)
        return np.load(path, mmap_mode='r' if mmap else None)


def load_cell(scenario: str, treatment: str, dt: int) -> Cell:
    """The outputs of one cell of the matrix."""
    return Cell(scenario, treatment, dt)


def align(cell_ids: np.ndarray, reference_ids: np.ndarray) -> np.ndarray:
    """The row of each reference river in the cell's arrays, -1 where the cell has no row for it."""
    if np.unique(cell_ids).shape[0] != cell_ids.shape[0]:
        raise ValueError('the cell must list each river once')
    order = np.argsort(cell_ids)
    position = np.searchsorted(cell_ids, reference_ids, sorter=order).clip(0, cell_ids.shape[0] - 1)
    rows = order[position]
    rows = np.where(cell_ids[rows] == reference_ids, rows, -1)
    if rows.shape != reference_ids.shape:
        raise ValueError('every reference river must get a row or -1')
    return rows


def peaks(series: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The largest hourly discharge of each river and the hour it occurs."""
    if series.ndim != 2 or series.shape[1] == 0:
        raise ValueError('series must be (river, hour)')
    peak_hour = np.argmax(series, axis=1)
    peak = np.take_along_axis(np.asarray(series), peak_hour[:, None], axis=1)[:, 0]
    return peak.astype(np.float64), peak_hour


def nyquist_amplitude(series: np.ndarray) -> np.ndarray:
    """
    The largest amplitude of an oscillation of period two hours in each series: the largest |Q[t+1] - 2Q[t] + Q[t-1]|/4,
    which is A for a series alternating Q +- A and nearly zero for a smooth one.
    """
    if series.ndim != 2 or series.shape[1] < 3:
        raise ValueError('series must be (river, hour) with at least three hours')
    values = np.asarray(series, dtype=np.float64)
    second = values[:, 2:] - 2 * values[:, 1:-1] + values[:, :-2]
    amplitude = np.abs(second).max(axis=1) / 4
    if np.any(amplitude < 0):
        raise ValueError('amplitudes cannot be negative')
    return amplitude


def dip_below_start(series: np.ndarray, peak_hour: np.ndarray) -> np.ndarray:
    """How far each series falls below its first hour before its peak: the dip of Muskingum ahead of a rise."""
    if series.shape[0] != peak_hour.shape[0]:
        raise ValueError('one peak hour per river')
    values = np.asarray(series, dtype=np.float64)
    before = np.arange(values.shape[1])[None, :] <= peak_hour[:, None]
    lowest = np.where(before, values, np.inf).min(axis=1)
    dip = np.maximum(values[:, 0] - lowest, 0)
    if not np.all(np.isfinite(dip)):
        raise ValueError('every series must have an hour before its peak')
    return dip


def oscillation_hours(series: np.ndarray) -> np.ndarray:
    """The hours of each series that close four hourly changes alternating in direction."""
    from .recorders import changes_and_oscillation

    if series.ndim != 2:
        raise ValueError('series must be (river, hour)')
    _, oscillating = changes_and_oscillation(np.asarray(series), series.shape[1])
    counts = oscillating.sum(axis=1)
    if counts.shape[0] != series.shape[0]:
        raise ValueError('one count per river')
    return counts


def save_table(frame, name: str) -> Path:
    """Write an analysis table as CSV to the tables directory."""
    if not name:
        raise ValueError('a table needs a name')
    config.TABLES.mkdir(parents=True, exist_ok=True)
    path = config.TABLES / f'{name}.csv'
    frame.to_csv(path, index=False, float_format='%.6g')
    if not path.exists():
        raise OSError(f'{path} was not written')
    return path
