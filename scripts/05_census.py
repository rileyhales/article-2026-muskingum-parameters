"""
Census of the hydrofabric: how travel times were assigned, how many rivers have a negative c1 (too long) or c3 (too
short) at each routing time step, where the short rivers sit in the topology, and what stabilizing them costs. Covers
the Columbia basin and every river of the global River Forecast System v3 hydrofabric. Their figures are drawn with
the other hydrofabrics by 16_hydrofabric_census.py.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/05_census.py
"""


import numpy as np
import pandas as pd
import river_route as rr

from stability import census, config, hydrofabric

GLOBAL_METADATA = config.RFS_ROOT / 'hydrography' / 'global' / 'metadata.parquet'


def sign_table(name: str, k: np.ndarray, x: np.ndarray) -> pd.DataFrame:
    """Counts and shares of rivers by coefficient sign at each routing time step."""
    rows = []
    for dt in config.DT_ROUTING:
        too_long, too_short = census.coefficient_signs(k, x, dt)
        rows.append({
            'domain': name, 'dt_s': dt, 'rivers': k.shape[0], 'too_long': int(too_long.sum()),
            'too_short': int(too_short.sum()), 'valid': int((~too_long & ~too_short).sum()),
        })
    table = pd.DataFrame(rows)
    for column in ('too_long', 'too_short', 'valid'):
        table[f'{column}_pct'] = 100 * table[column] / table['rivers']
    if not np.all(table[['too_long', 'too_short', 'valid']].sum(axis=1) == table['rivers']):
        raise ValueError('every river must fall in exactly one class')
    return table


def stabilization_cost(network: rr.Network) -> pd.DataFrame:
    """The sub-reaches and subcycles river-route's stabilized network routes at each dt, and the work they cost."""
    rows = []
    for dt in config.DT_ROUTING:
        substeps, subcycles, resolvable = network.conditioning(float(dt))
        rows.append({
            'dt_s': dt, 'reaches': int(substeps.sum()), 'substepped_rivers': int((substeps > 1).sum()),
            'max_substeps': int(substeps.max()), 'subcycled_rivers': int((subcycles > 1).sum()),
            'max_subcycles': int(subcycles.max()), 'unresolvable': int((~resolvable).sum()),
            'work_vs_standard': float((substeps * subcycles).sum() / network.size),
            'work_vs_standard_at_3600': float((substeps * subcycles).sum() / network.size * 3600 / dt),
        })
    if not rows:
        raise ValueError('no time steps')
    return pd.DataFrame(rows)


def short_river_topology(table: pd.DataFrame) -> pd.DataFrame:
    """Where the rivers too short for each dt sit: between two confluences, below one, above one, or alone."""
    place = census.positions(table['down'].to_numpy())
    rows = []
    k = table['muskingumK'].to_numpy(dtype=np.float64)
    x = table['muskingumX'].to_numpy(dtype=np.float64)
    for dt in config.DT_ROUTING:
        _, short = census.coefficient_signs(k, x, dt)
        counts = np.bincount(place[short], minlength=len(census.POSITIONS))
        rows.append({
            'dt_s': dt, 'too_short': int(short.sum()),
            **{position: int(count) for position, count in zip(census.POSITIONS, counts, strict=True)},
            'median_length_m': float(np.median(table['Length'].to_numpy()[short])) if short.any() else np.nan,
        })
    result = pd.DataFrame(rows)
    if np.any(result['too_short'] < result['between_confluences']):
        raise ValueError('a class cannot hold more rivers than are short')
    return result


def by_order(table: pd.DataFrame, k_column: str, length_column: str, order_column: str) -> pd.DataFrame:
    """Length, velocity, and travel time of the rivers of each Strahler order."""
    grouped = table.groupby(order_column)
    summary = grouped.agg(
        rivers=(k_column, 'size'), velocity_m_s=('velocity_factor', 'first'),
        length_median_m=(length_column, 'median'), length_p01_m=(length_column, lambda s: s.quantile(0.01)),
        k_median_s=(k_column, 'median'), k_min_s=(k_column, 'min'), k_max_s=(k_column, 'max'),
    )
    if summary['rivers'].sum() != len(table):
        raise ValueError('every river must have an order')
    return summary.reset_index()


if __name__ == '__main__':
    config.TABLES.mkdir(parents=True, exist_ok=True)
    basin = hydrofabric.load()
    world = pd.read_parquet(GLOBAL_METADATA, columns=['strahlerOrder', 'Length', 'musk_k', 'musk_x', 'velocity_factor'])
    signs = pd.concat([
        sign_table('Columbia', basin['muskingumK'].to_numpy(np.float64), basin['muskingumX'].to_numpy(np.float64)),
        sign_table('Global', world['musk_k'].to_numpy(np.float64), world['musk_x'].to_numpy(np.float64)),
    ], ignore_index=True)
    signs.to_csv(config.TABLES / 'census_coefficient_signs.csv', index=False)
    cost = stabilization_cost(rr.Network(config.NETWORK_FILE))
    cost.to_csv(config.TABLES / 'census_stabilization_cost.csv', index=False)
    topology = short_river_topology(basin)
    topology.to_csv(config.TABLES / 'census_short_river_topology.csv', index=False)
    orders = by_order(basin, 'muskingumK', 'Length', 'strahlerOrder')
    orders.to_csv(config.TABLES / 'census_columbia_by_order.csv', index=False)
    world_orders = by_order(world, 'musk_k', 'Length', 'strahlerOrder')
    world_orders.to_csv(config.TABLES / 'census_global_by_order.csv', index=False)
    for name, frame in (('signs', signs), ('cost', cost), ('topology', topology), ('Columbia orders', orders),
                        ('global orders', world_orders)):
        print(f'\n{name}\n{frame.to_string(index=False)}')
    print('\nglobal x values:', np.unique(world['musk_x'].to_numpy()))
