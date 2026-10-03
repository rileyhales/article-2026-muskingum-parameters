"""
Census of the hydrofabric: how travel times were assigned, how many rivers have a negative c1 (too long) or c3 (too
short) at each routing time step, where the short rivers sit in the topology, and what stabilizing them costs. Covers
the Columbia basin and every river of the global River Forecast System v3 hydrofabric.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/05_census.py
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import river_route as rr

from stability import config, hydrofabric, plotting

GLOBAL_METADATA = config.RFS_ROOT / 'hydrography' / 'global' / 'metadata.parquet'


def coefficient_signs(k: np.ndarray, x: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Masks of the rivers too long (c1 < 0, dt < 2kx) and too short (c3 < 0, dt > 2k(1-x)) for dt."""
    if dt <= 0:
        raise ValueError('dt must be positive')
    if k.shape != x.shape:
        raise ValueError('k and x must be one per river')
    return 2 * k * x > dt, 2 * k * (1 - x) < dt


def sign_table(name: str, k: np.ndarray, x: np.ndarray) -> pd.DataFrame:
    """Counts and shares of rivers by coefficient sign at each routing time step."""
    rows = []
    for dt in config.DT_ROUTING:
        too_long, too_short = coefficient_signs(k, x, dt)
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
    down = table['down'].to_numpy()
    n_up = np.bincount(down[down >= 0], minlength=len(table))
    below_confluence = n_up >= 2
    above_confluence = np.zeros(len(table), dtype=bool)
    above_confluence[down >= 0] = n_up[down[down >= 0]] >= 2
    rows = []
    k = table['muskingumK'].to_numpy(dtype=np.float64)
    x = table['muskingumX'].to_numpy(dtype=np.float64)
    for dt in config.DT_ROUTING:
        _, short = coefficient_signs(k, x, dt)
        rows.append({
            'dt_s': dt, 'too_short': int(short.sum()),
            'between_confluences': int((short & below_confluence & above_confluence).sum()),
            'below_confluence_only': int((short & below_confluence & ~above_confluence).sum()),
            'above_confluence_only': int((short & ~below_confluence & above_confluence & (n_up > 0)).sum()),
            'headwater': int((short & (n_up == 0)).sum()),
            'chain': int((short & (n_up == 1) & ~above_confluence).sum()),
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


def plot_signs(signs: pd.DataFrame) -> None:
    """Share of rivers too long, valid, and too short at each dt, for the Columbia and the global hydrofabric."""
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    classes = (('too_long_pct', 'c1 < 0 (too long)', '#eb6834'), ('valid_pct', 'All coefficients ≥ 0', '#2a78d6'),
               ('too_short_pct', 'c3 < 0 (too short)', '#4a3aa7'))
    for axis, (domain, rows) in zip(axes, signs.groupby('domain', sort=False), strict=True):
        steps = rows['dt_s'].to_numpy()
        for column, label, color in classes:
            axis.plot(steps, rows[column], marker='o', markersize=4, color=color, label=label)
        axis.set_xscale('log')
        axis.set_xticks(steps, [f'{s // 60} min' if s >= 60 else f'{s} s' for s in steps], rotation=45)
        axis.set_title(f'{domain} ({rows["rivers"].iloc[0]:,} rivers)')
        axis.set_xlabel('Routing time step')
    axes[0].set_ylabel('Share of rivers (%)')
    axes[0].legend(loc='center left')
    plotting.save(figure, 'census_coefficient_signs')
    return


def plot_travel_times(columbia: pd.DataFrame, global_k: np.ndarray) -> None:
    """Distribution of k with the window of k each time step keeps all coefficients non-negative."""
    figure, axis = plt.subplots(figsize=(7.2, 3.0))
    for k, label, color in ((global_k, 'Global hydrofabric', '#86b6ef'),
                            (columbia['muskingumK'].to_numpy(), 'Columbia', '#2a78d6')):
        ordered = np.sort(k)
        axis.plot(ordered, np.linspace(0, 1, ordered.shape[0]), color=color, label=label)
    x = 0.2
    for dt, shade in ((60, '#f0efec'), (900, '#e4e3df'), (3600, '#d6d5d0')):
        low, high = dt / (2 * (1 - x)), dt / (2 * x)
        axis.axvspan(low, high, color=shade, zorder=0, linewidth=0)
        axis.text(np.sqrt(low * high), 0.5, f'Δt = {dt // 60} min', ha='center', va='center', color=plotting.MUTED,
                  fontsize=8, rotation=90, transform=axis.get_xaxis_transform())
    axis.set_xscale('log')
    axis.set_xlabel('Muskingum k, the travel time of a river (s)')
    axis.set_ylabel('Cumulative share of rivers')
    axis.legend(loc='upper left')
    plotting.save(figure, 'census_travel_time_distribution')
    return


if __name__ == '__main__':
    plotting.apply_style()
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
    plot_signs(signs)
    plot_travel_times(basin, world['musk_k'].to_numpy(np.float64))
    for name, frame in (('signs', signs), ('cost', cost), ('topology', topology), ('Columbia orders', orders),
                        ('global orders', world_orders)):
        print(f'\n{name}\n{frame.to_string(index=False)}')
    print('\nglobal x values:', np.unique(world['musk_x'].to_numpy()))
