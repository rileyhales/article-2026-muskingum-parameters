"""
Choose the event windows every ERA5 cell of the matrix keeps the whole network's hourly series for, from the
reference run's observation series: the largest peak of the period at the Columbia outlet, at the mouths of the Snake
and the Willamette, and the flashiest peak at the mouth of a smaller tributary, the one whose peak most exceeds its
median flow. Each window opens WINDOW_BEFORE days before its peak and closes WINDOW_AFTER days after it. The windows of
data/inputs/events.json were chosen from the 2000-2019 reference and all fall within 2002-2011.

Run with the project environment:  uv run python scripts/07_select_events.py
"""

import json

import numpy as np
import pandas as pd

from stability import config, hydrofabric

WINDOW_BEFORE = 30
WINDOW_AFTER = 15
EVENTS_FILE = config.INPUTS / 'events.json'
FLASHY_CANDIDATES = ('John Day', 'Deschutes', 'Yakima', 'Cowlitz', 'Okanogan', 'Spokane', 'Kettle', 'Clearwater')


def series_of(mouth_id: int, ids: np.ndarray, series: np.ndarray) -> np.ndarray:
    """The reference hourly series at a tributary mouth."""
    rows = np.flatnonzero(ids == mouth_id)
    if rows.shape[0] != 1:
        raise ValueError(f'river {mouth_id} must be one observation river')
    values = np.asarray(series[rows[0]], dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError('the reference series must be finite')
    return values


def window_around(peak_hour: int) -> tuple[str, str]:
    """The window of dates around a peak, whole days, clipped to the study period."""
    start = np.datetime64(f'{config.ERA5_YEARS[0]}-01-01T00:00:00') + np.timedelta64(int(peak_hour), 'h')
    first = max(start.astype('datetime64[D]') - np.timedelta64(WINDOW_BEFORE, 'D'),
                np.datetime64(f'{config.ERA5_YEARS[0]}-01-01'))
    last = min(start.astype('datetime64[D]') + np.timedelta64(WINDOW_AFTER, 'D'),
               np.datetime64(f'{config.ERA5_YEARS[-1] + 1}-01-01'))
    if last <= first:
        raise ValueError('a window must have hours')
    return str(first), str(last)


if __name__ == '__main__':
    reference = config.run_dir(config.ERA5_SCENARIO, *config.REFERENCE)
    ids = np.load(reference / 'observation_ids.npy')
    series = np.load(reference / 'observations.npy', mmap_mode='r')
    events = {}
    for name in ('Columbia', 'Snake', 'Willamette'):
        flow = series_of(hydrofabric.TRIBUTARIES[name], ids, series)
        events[f'{name.lower()}-peak'] = window_around(int(np.argmax(flow)))
    flashiness = {}
    for name in FLASHY_CANDIDATES:
        flow = series_of(hydrofabric.TRIBUTARIES[name], ids, series)
        flashiness[name] = (flow.max() / np.median(flow), int(np.argmax(flow)))
    names = list(flashiness)
    flashiest = names[int(np.argmax([flashiness[n][0] for n in names]))]
    events[f'{flashiest.lower().replace(" ", "-")}-flash'] = window_around(flashiness[flashiest][1])
    print(pd.DataFrame(flashiness, index=['peak/median', 'hour']).T.to_string())
    EVENTS_FILE.write_text(json.dumps(events, indent=2))
    print(json.dumps(events, indent=2))
