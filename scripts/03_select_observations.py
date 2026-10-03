"""
Choose the observation rivers whose whole hourly series every 20-year simulation keeps: points every SPACING_KM along
the main stem of each major tributary of the Columbia, and the mouth of each. Every other river keeps annual
statistics, and every river keeps its series inside the event windows.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/03_select_observations.py
"""


import numpy as np
import pandas as pd

from stability import config, hydrofabric

SPACING_KM = 20.0
OBSERVATIONS_FILE = config.INPUTS / 'observation_rivers.parquet'


def stem_points(table: pd.DataFrame, name: str, mouth_id: int) -> pd.DataFrame:
    """Rows every SPACING_KM up the main stem of a tributary, with the channel distance of each above the mouth."""
    rows = np.flatnonzero(table['riverId'].to_numpy() == mouth_id)
    if rows.shape[0] != 1:
        raise ValueError(f'the mouth of {name} must be one river of the network')
    stem = hydrofabric.main_stem(table, int(rows[0]))
    outlet_km = table['km_to_outlet'].to_numpy()[stem] - table['Length'].to_numpy()[stem] / 1000
    above_mouth = outlet_km - outlet_km[-1]
    marks = np.arange(0, above_mouth.max(), SPACING_KM)
    # the first river, walking up from the mouth, whose outlet is at or above each mark
    picks = {int(stem[np.flatnonzero(above_mouth >= mark)[-1]]) for mark in marks} | {int(stem[0]), int(stem[-1])}
    if not picks:
        raise ValueError(f'{name} has no observation points')
    points = table.loc[sorted(picks), ['riverId', 'area_km2', 'muskingumK', 'Length', 'lat', 'lon']].copy()
    points['row'] = points.index.to_numpy()
    points['stem'] = name
    distance_of_row = pd.Series(above_mouth, index=stem)
    points['km_above_mouth'] = distance_of_row.reindex(points['row'].to_numpy()).to_numpy()
    return points


if __name__ == '__main__':
    columbia = hydrofabric.load()
    observations = pd.concat(
        [stem_points(columbia, name, mouth) for name, mouth in hydrofabric.TRIBUTARIES.items()], ignore_index=True
    )
    observations = observations.drop_duplicates('riverId').sort_values('row').reset_index(drop=True)
    observations.to_parquet(OBSERVATIONS_FILE, index=False)
    print(observations.groupby('stem').size().to_string())
    print(f'{len(observations)} observation rivers written to {OBSERVATIONS_FILE}')
