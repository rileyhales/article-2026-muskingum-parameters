"""
Recorders receive each chunk of routed discharge from ``engine.route`` and keep what the analysis reads.

- ``SeriesRecorder`` keeps the whole hourly series of chosen rivers (every river by default) in a .npy file.
- ``AnnualStatsRecorder`` keeps, for every river and calendar year, the statistics of the hourly series that measure
  peaks and the artifacts of negative coefficients (negative discharge, oscillation), carrying the last hours of each
  chunk into the next so a statistic that spans a chunk boundary is counted once.
"""

from pathlib import Path

import numpy as np

__all__ = [
    'SeriesRecorder',
    'AnnualStatsRecorder',
    'WindowRecorder',
    'TeeRecorder',
    'STATISTICS',
    'changes_and_oscillation',
]

# hourly changes smaller than this fraction of the discharge are float32 rounding, not a change of direction
RELATIVE_TOLERANCE = 1e-4
STATISTICS = (
    'max', 'argmax', 'min', 'volume', 'negative_hours', 'negative_volume', 'oscillation_hours', 'total_variation'
)


class SeriesRecorder:
    """Writes the hourly discharge of the chosen rows into a (river, hour) float32 .npy file as it is routed."""

    def __init__(self, path: Path, n_rows: int, n_hours: int, rows: np.ndarray | None = None) -> None:
        if n_rows <= 0 or n_hours <= 0:
            raise ValueError('the series needs at least one river and one hour')
        if rows is not None and (rows.min() < 0 or rows.max() >= n_rows):
            raise ValueError('rows must index the routed rivers')
        self.rows = rows
        n_kept = n_rows if rows is None else rows.shape[0]
        self.series = np.lib.format.open_memmap(path, mode='w+', dtype=np.float32, shape=(n_kept, n_hours))
        self.filled = 0
        return

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        n = dates.shape[0]
        if self.filled + n > self.series.shape[1]:
            raise ValueError('more hours were routed than the series holds')
        if discharge.shape[1] != n:
            raise ValueError('the discharge must have one column per date')
        self.series[:, self.filled : self.filled + n] = discharge if self.rows is None else discharge[self.rows]
        self.filled += n
        return

    def close(self) -> None:
        if self.filled != self.series.shape[1]:
            raise ValueError(f'{self.filled} of {self.series.shape[1]} hours were routed')
        self.series.flush()
        return


def changes_and_oscillation(extended: np.ndarray, n_new: int) -> tuple[np.ndarray, np.ndarray]:
    """
    The absolute hourly change into each of the last n_new hours of ``extended`` (river, hour), and whether each of
    those hours closes an oscillation: four hourly changes in a row that alternate in direction, up-down-up-down or
    down-up-down-up. A natural peak changes direction once, and even the dip of Muskingum ahead of a sharp rise (down,
    up, then down after the peak) changes direction twice; only an oscillation at the period of two steps alternates
    three times running.
    """
    if extended.ndim != 2 or n_new > extended.shape[1]:
        raise ValueError('extended must be (river, hour) and hold the new hours')
    if n_new <= 0:
        raise ValueError('there must be new hours to evaluate')
    change = np.diff(extended.astype(np.float64), axis=1)
    scale = np.maximum(np.abs(extended[:, 1:]), 1e-3)
    change = np.where(np.abs(change) > RELATIVE_TOLERANCE * scale, change, 0.0)
    n_changes = change.shape[1]
    first = n_changes - n_new  # change index of the change into the first new hour
    total = np.zeros((extended.shape[0], n_new))
    oscillating = np.zeros((extended.shape[0], n_new), dtype=bool)
    usable = slice(max(first, 0), n_changes)
    total[:, usable.start - first :] = np.abs(change[:, usable])
    last_four = np.arange(max(first, 3), n_changes)
    if last_four.shape[0]:
        d1, d2 = change[:, last_four - 3], change[:, last_four - 2]
        d3, d4 = change[:, last_four - 1], change[:, last_four]
        oscillating[:, last_four - first] = (d1 * d2 < 0) & (d2 * d3 < 0) & (d3 * d4 < 0)
    return total, oscillating


class AnnualStatsRecorder:
    """Per river and year: peak and its hour, minimum, volume, negative hours and volume, oscillation, variation."""

    def __init__(self, n_rivers: int, years: tuple[int, ...]) -> None:
        if n_rivers <= 0 or not years:
            raise ValueError('statistics need rivers and years')
        if list(years) != sorted(years):
            raise ValueError('years must be in order')
        self.years = np.asarray(years)
        shape = (len(years), n_rivers)
        self.stats = {name: np.zeros(shape, dtype=np.float64) for name in STATISTICS}
        self.stats['max'][:] = -np.inf
        self.stats['min'][:] = np.inf
        self.carry = np.zeros((n_rivers, 0), dtype=np.float32)
        self.hours = 0
        return

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        year_of_hour = dates.astype('datetime64[Y]').astype(int) + 1970
        if not np.all(np.isin(year_of_hour, self.years)):
            raise ValueError('every hour must fall in a recorded year')
        if discharge.shape[1] != dates.shape[0]:
            raise ValueError('the discharge must have one column per date')
        changes, oscillating = changes_and_oscillation(np.concatenate([self.carry, discharge], axis=1), dates.shape[0])
        for year in np.unique(year_of_hour):
            cols = np.flatnonzero(year_of_hour == year)
            self._accumulate(int(np.searchsorted(self.years, year)), cols, discharge, changes, oscillating)
        self.carry = discharge[:, -4:].copy()
        self.hours += dates.shape[0]
        return

    def _accumulate(self, y: int, cols: np.ndarray, q: np.ndarray, changes: np.ndarray, osc: np.ndarray) -> None:
        """Fold the hours ``cols`` of one chunk into year index y."""
        if cols.shape[0] == 0 or y >= self.years.shape[0]:
            raise ValueError('a year must have hours to accumulate')
        block = q[:, cols]
        block_max = block.max(axis=1)
        newer = block_max > self.stats['max'][y]
        self.stats['argmax'][y] = np.where(newer, self.hours + cols[block.argmax(axis=1)], self.stats['argmax'][y])
        self.stats['max'][y] = np.maximum(self.stats['max'][y], block_max)
        self.stats['min'][y] = np.minimum(self.stats['min'][y], block.min(axis=1))
        self.stats['volume'][y] += block.astype(np.float64).sum(axis=1) * 3600
        self.stats['negative_hours'][y] += np.count_nonzero(block < 0, axis=1)
        self.stats['negative_volume'][y] -= np.minimum(block, 0).astype(np.float64).sum(axis=1) * 3600
        self.stats['oscillation_hours'][y] += np.count_nonzero(osc[:, cols], axis=1)
        self.stats['total_variation'][y] += changes[:, cols].sum(axis=1)
        return

    def save(self, path: Path) -> None:
        if not np.all(np.isfinite(self.stats['max'])):
            raise ValueError('every year must have been routed')
        np.savez_compressed(path, years=self.years, **self.stats)
        return


class TeeRecorder:
    """Hands each chunk to several recorders."""

    def __init__(self, *recorders) -> None:
        if not recorders:
            raise ValueError('a tee needs at least one recorder')
        self.recorders = recorders
        return

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        for recorder in self.recorders:
            recorder.record(dates, discharge)
        return


class WindowRecorder:
    """Writes the hourly discharge of every routed river inside each named window of dates to ``<dir>/<name>.npy``."""

    def __init__(self, directory: Path, windows: dict[str, tuple[str, str]], n_rows: int) -> None:
        if not windows:
            raise ValueError('a window recorder needs windows')
        if n_rows <= 0:
            raise ValueError('a window recorder needs rivers')
        directory.mkdir(parents=True, exist_ok=True)
        self.windows = {}
        for name, (start, end) in windows.items():
            first, last = np.datetime64(start, 's'), np.datetime64(end, 's')
            hours = int((last - first) / np.timedelta64(3600, 's'))
            if hours <= 0:
                raise ValueError(f'window {name} must end after it starts')
            series = np.lib.format.open_memmap(directory / f'{name}.npy', mode='w+', dtype=np.float32,
                                               shape=(n_rows, hours))
            self.windows[name] = (first, hours, series)
        return

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        if discharge.shape[1] != dates.shape[0]:
            raise ValueError('the discharge must have one column per date')
        if dates.shape[0] == 0:
            raise ValueError('a chunk must hold hours')
        for first, hours, series in self.windows.values():
            position = ((dates - first) / np.timedelta64(3600, 's')).astype(np.int64)
            inside = (position >= 0) & (position < hours)
            if inside.any():
                series[:, position[inside]] = discharge[:, inside]
        return

    def close(self) -> None:
        for _, _, series in self.windows.values():
            series.flush()
        return
