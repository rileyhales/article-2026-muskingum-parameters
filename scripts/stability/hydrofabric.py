"""
The Columbia hydrofabric as a table: topology, Muskingum parameters, stream metadata, and the derived quantities the
analysis measures along the network, such as the distance and travel time from each river to the outlet, the main stem
of a tributary, and the rivers downstream of a river at chosen distances.
"""

import numpy as np
import pandas as pd

from . import config

__all__ = ['load', 'main_stem', 'downstream_rows', 'rows_at_distances', 'TRIBUTARIES', 'tributary_mouths']

# the major tributaries of the Columbia, each named by the riverId of its mouth reach in the hydrofabric. The named ones
# come from the River Forecast System river names; the others were identified by drainage area and mouth location.
TRIBUTARIES = {
    'Columbia': config.COLUMBIA_OUTLET,  # the max-area path, which follows the Kootenay above Castlegar
    'Upper Columbia': 720201501,  # the Columbia above its confluence with the Kootenay at Castlegar
    'Kootenay': 720187612,
    'Pend Oreille': 720213439,
    'Kettle': 720149983,
    'Spokane': 720195705,
    'Okanogan': 720169881,
    'Yakima': 720164131,
    'Snake': 720221669,
    'Clearwater': 720197840,
    'Salmon': 720174090,
    'John Day': 720229661,
    'Deschutes': 720146339,
    'Willamette': 720184033,
    'Cowlitz': 720221682,
}


def load() -> pd.DataFrame:
    """The Columbia network table joined to its stream metadata, with derived columns, in network (DFS) order."""
    network = pd.read_parquet(config.NETWORK_FILE)
    metadata = pd.read_parquet(config.METADATA_FILE)
    if not np.array_equal(network['riverId'].to_numpy(), metadata['riverId'].to_numpy()):
        raise ValueError('the network and metadata must list the same rivers in the same order')
    table = network.join(metadata.drop(columns=['riverId', 'nextRiverId']))
    index = pd.Series(np.arange(len(table)), index=table['riverId'].to_numpy())
    has_down = table['nextRiverId'].to_numpy() != -1
    down = np.full(len(table), -1, dtype=np.int64)
    down[has_down] = index.reindex(table['nextRiverId'].to_numpy()[has_down]).to_numpy()
    if np.any(down[has_down] <= np.flatnonzero(has_down)):
        raise ValueError('the table must be in topological order')
    table['down'] = down
    length_km = table['Length'].to_numpy() / 1000
    k = table['muskingumK'].to_numpy(dtype=np.float64)
    to_outlet_km = np.zeros(len(table))
    to_outlet_s = np.zeros(len(table))
    for i in range(len(table) - 1, -1, -1):  # downstream rivers come later in the table, so they are already done
        d = down[i]
        to_outlet_km[i] = length_km[i] + (to_outlet_km[d] if d >= 0 else 0)
        to_outlet_s[i] = k[i] + (to_outlet_s[d] if d >= 0 else 0)
    table['km_to_outlet'] = to_outlet_km  # from the upstream end of the river to the basin outlet
    table['travel_s_to_outlet'] = to_outlet_s
    table['area_km2'] = table['DSContArea'].to_numpy() / 1e6
    return table


def main_stem(table: pd.DataFrame, mouth: int) -> np.ndarray:
    """Rows of the main stem ending at row ``mouth``, from its headwater down: always the upstream of largest area."""
    if not 0 <= mouth < len(table):
        raise ValueError(f'row {mouth} is not in the table')
    area = table['area_km2'].to_numpy()
    down = table['down'].to_numpy()
    upstreams = pd.Series(np.arange(len(table)))[down >= 0].groupby(down[down >= 0]).agg(list)
    stem = [mouth]
    for _ in range(len(table)):  # a path can visit each river at most once
        ups = upstreams.get(stem[-1])
        if ups is None:
            break
        stem.append(int(ups[int(np.argmax(area[ups]))]))
    if len(stem) > len(table):
        raise ValueError('the main stem must end at a headwater')
    return np.asarray(stem[::-1])


def downstream_rows(table: pd.DataFrame, start: int) -> np.ndarray:
    """Rows from ``start`` down to the basin outlet, in flow order."""
    down = table['down'].to_numpy()
    if not 0 <= start < len(table):
        raise ValueError(f'row {start} is not in the table')
    path = [start]
    for _ in range(len(table)):
        if down[path[-1]] < 0:
            break
        path.append(int(down[path[-1]]))
    return np.asarray(path)


def rows_at_distances(table: pd.DataFrame, start: int, distances_km: tuple[float, ...]) -> dict[float, int]:
    """
    For each distance, the river at that distance of channel downstream of the outlet of ``start``: the first river
    whose outlet lies at least that far below it. Distances past the basin outlet are left out.
    """
    if any(d < 0 for d in distances_km):
        raise ValueError('distances must not be negative')
    path = downstream_rows(table, start)
    outlet_km = table['km_to_outlet'].to_numpy() - table['Length'].to_numpy() / 1000  # outlet of each river
    below = outlet_km[path[0]] - outlet_km[path]  # channel distance from the outlet of start to the outlet of each
    found = {}
    for distance in distances_km:
        at = np.flatnonzero(below >= distance)
        if at.shape[0]:
            found[distance] = int(path[at[0]])
    if found and min(found.values()) < 0:
        raise ValueError('rows must be in the table')
    return found


def tributary_mouths(table: pd.DataFrame, stem: np.ndarray, min_area_km2: float) -> pd.DataFrame:
    """The rivers draining into the main stem ``stem`` whose drainage area is at least min_area_km2."""
    if min_area_km2 <= 0:
        raise ValueError('the area threshold must be positive')
    on_stem = np.zeros(len(table), dtype=bool)
    on_stem[stem] = True
    down = table['down'].to_numpy()
    joins = (~on_stem) & (down >= 0) & on_stem[np.clip(down, 0, None)]
    mouths = table[joins & (table['area_km2'].to_numpy() >= min_area_km2)]
    if mouths.empty:
        raise ValueError('no tributary is that large')
    return mouths.sort_values('area_km2', ascending=False)
