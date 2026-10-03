"""
How far downstream does one river with a negative c3 disturb the discharge, and how far downstream does splitting one
river with a negative c1 into substeps change it? Every river too short for dt is a seed of the first kind; a fixed
random sample of the rivers too long for dt are seeds of the second.

For short seeds the network is routed with every short river subcycled, and again with a batch of seeds left untreated. Routing is
linear, so the difference between the two at a river downstream of a seed is exactly what that seed contributes, as
long as no other seed of its batch lies upstream of that river: seeds are batched so that none lies in the watershed
of another's farthest downstream point. The difference is measured at fixed channel distances below each seed.

For long seeds the network is routed standard, and again with a batch of seeds split into substeps, x held or x
adjusted, and every other river whole: the difference is what treating that one river changes downstream.

Writes tables/single_defect.csv, one row per kind, scenario, dt, seed, and distance, and figures
single_defect_distance and single_defect_long_distance.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/09_single_defect.py
"""

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import river_route as rr

from stability import (
    config,
    engine,
    forcing,
    hydrofabric,
    metrics,
    plotting,
    treatments,
)

DISTANCES_KM = (0.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0)
STEPS = (3600, 1800)
LONG_STEPS = (900, 3600)
LONG_KINDS = {'long-substeps': 'substeps', 'long-substeps-xadj': 'substeps-xadj'}
LONG_SEEDS = 400  # long reaches are too many to all be seeds; a fixed random sample of this many is drawn
LONG_SEED = 20261002
SPINUP_DAYS = 60
EVENTS = ('willamette-peak', 'columbia-peak')


def seed_paths(table: pd.DataFrame, seeds: np.ndarray) -> dict[int, dict[float, int]]:
    """For every seed row, the row at each distance below it that the basin reaches."""
    paths = {int(s): hydrofabric.rows_at_distances(table, int(s), DISTANCES_KM) for s in seeds}
    if any(0.0 not in path for path in paths.values()):
        raise ValueError('every seed is its own point at distance zero')
    return paths


def batch_seeds(table: pd.DataFrame, paths: dict[int, dict[float, int]]) -> list[np.ndarray]:
    """Group seeds so that no seed of a batch lies in the watershed of another seed's farthest point."""
    upstream = table['upstreamCount'].to_numpy()
    reach = {seed: (max(path.values()) - upstream[max(path.values())], max(path.values())) for seed, path in
             paths.items()}
    batches: list[list[int]] = []
    for seed in sorted(paths):
        lo, hi = reach[seed]
        placed = False
        for batch in batches:
            clash = any(lo <= other <= hi or reach[other][0] <= seed <= reach[other][1] for other in batch)
            if not clash:
                batch.append(seed)
                placed = True
                break
        if not placed:
            batches.append([seed])
    if sum(len(b) for b in batches) != len(paths):
        raise ValueError('every seed must be in exactly one batch')
    return [np.asarray(b) for b in batches]


class KeepRows:
    """Keeps the whole series of chosen rows of every chunk."""

    def __init__(self, rows: np.ndarray) -> None:
        if rows.ndim != 1 or rows.shape[0] == 0:
            raise ValueError('choose at least one row')
        self.rows, self.chunks, self.dates = rows, [], []

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        self.chunks.append(discharge[self.rows].copy())
        self.dates.append(dates.copy())

    def series(self) -> tuple[np.ndarray, np.ndarray]:
        if not self.chunks:
            raise ValueError('nothing was routed')
        return np.concatenate(self.dates), np.concatenate(self.chunks, axis=1)


def forcing_of(scenario: str, table: pd.DataFrame, areas: np.ndarray) -> tuple[list, np.ndarray | None, tuple]:
    """The chunks of a scenario, the baseflow volumes to start from (None from a dry start), and the window kept."""
    if scenario in forcing.SYNTHETIC_PERTURBATIONS:
        chunks = list(forcing.synthetic_chunks(scenario, areas))
        return chunks, forcing.synthetic_baseflow_volumes(areas), (chunks[0].dates[0], chunks[-1].dates[-1])
    windows = json.loads((config.INPUTS / 'events.json').read_text())
    first, last = (np.datetime64(d, 's') for d in windows[scenario])
    start = (first - np.timedelta64(SPINUP_DAYS, 'D')).astype('datetime64[M]')
    files = forcing.era5_files(str(start), str((last - np.timedelta64(1, 's')).astype('datetime64[M]')))
    return list(forcing.era5_chunks(files, table['riverId'].to_numpy())), None, (first, last)


def route_rows(treated: treatments.Treated, dt: int, chunks: list, base: np.ndarray | None, rows: np.ndarray):
    """The series of the chosen rows routed over a treated network, from steady baseflow or a dry start."""
    keep = KeepRows(rows)
    state = None if base is None else engine.steady_state(treated, dt, config.DT_RUNOFF, base)
    engine.route(treated, dt, config.DT_RUNOFF, iter(chunks), keep, threads=4, initial_state=state)
    return keep.series()


def seed_network(kind: str, dt: int, table: pd.DataFrame, batch: np.ndarray | None) -> treatments.Treated:
    """
    The network a kind of experiment routes: for short seeds every short reach subcycled, except the batch when given;
    for long seeds the standard network, with only the batch split into substeps (x held or x adjusted) when given.
    """
    every_river = np.ones(len(table), dtype=bool)
    mask = np.zeros(len(table), dtype=bool)
    if batch is not None:
        mask[batch] = True
    if kind == 'short':
        return treatments.Treated(treatments.ConditionedNetwork(table, 'subcycles', dt, mask), 'stabilized', None,
                                  every_river)
    if kind not in LONG_KINDS:
        raise ValueError(f'unknown kind {kind}')
    if batch is None:
        return treatments.Treated(rr.Network(table), 'standard', None, every_river)
    return treatments.Treated(treatments.ConditionedNetwork(table, LONG_KINDS[kind], dt, ~mask), 'stabilized', None,
                              every_river)


def seeds_of(kind: str, dt: int, table: pd.DataFrame) -> np.ndarray:
    """Rows of the seeds: every reach too short for dt, or a fixed random sample of the reaches too long for it."""
    too_long, too_short = rr.Network(table).unstable_mask(float(dt))
    has_downstream = table['down'].to_numpy() >= 0
    if kind == 'short':
        return np.flatnonzero(too_short & has_downstream)
    candidates = np.flatnonzero(too_long & has_downstream)
    if candidates.shape[0] == 0:
        raise ValueError('no reach is too long for dt')
    picked = np.random.default_rng(LONG_SEED).choice(candidates, min(LONG_SEEDS, candidates.shape[0]), replace=False)
    return np.sort(picked)


def difference_row(good: np.ndarray, bad: np.ndarray, seed_scale: float) -> dict:
    """The measures of the difference ``bad - good`` at one reach, scaled by its range and by its seed's range."""
    scale = max(good.max() - good.min(), 1e-6)
    alternating = float(metrics.nyquist_amplitude((bad - good)[None, :])[0])
    if not np.isfinite(alternating):
        raise ValueError('the difference must be finite')
    return {
        'max_difference': float(np.abs(bad - good).max() / scale),
        'max_difference_of_seed': float(np.abs(bad - good).max() / seed_scale),
        'alternating': alternating / scale, 'alternating_of_seed': alternating / seed_scale,
        'peak_difference': float((bad.max() - good.max()) / scale),
        'peak_relative': float((bad.max() - good.max()) / max(good.max(), 1e-6)),
        'timing_difference_h': int(np.argmax(bad) - np.argmax(good)),
        'negative_hours': int(np.count_nonzero(bad < 0)), 'treated_negative_hours': int(np.count_nonzero(good < 0)),
        'most_negative': float(min(bad.min(), 0) / scale),
    }


def measure(kind: str, scenario: str, dt: int, table: pd.DataFrame, areas: np.ndarray) -> pd.DataFrame:
    """The difference each seed makes at each distance below it, for one kind of seed, scenario, and dt."""
    paths = seed_paths(table, seeds_of(kind, dt, table))
    chunks, base, (first, last) = forcing_of(scenario, table, areas)
    every_row = np.unique([r for path in paths.values() for r in path.values()])
    dates, baseline = route_rows(seed_network(kind, dt, table, None), dt, chunks, base, every_row)
    window = (dates >= first) & (dates < last)
    column_of = {int(r): i for i, r in enumerate(every_row)}
    rows = []
    for batch in batch_seeds(table, paths):
        _, perturbed = route_rows(seed_network(kind, dt, table, batch), dt, chunks, base, every_row)
        for seed in batch:
            own = baseline[column_of[int(seed)], window].astype(np.float64)
            seed_scale = max(own.max() - own.min(), 1e-6)
            for distance, row in paths[int(seed)].items():
                good = baseline[column_of[row], window].astype(np.float64)
                bad = perturbed[column_of[row], window].astype(np.float64)
                rows.append({
                    'kind': kind, 'scenario': scenario, 'dt': dt, 'seed': int(table['riverId'].iloc[seed]),
                    'distance_km': distance, 'river': int(table['riverId'].iloc[row]),
                    'area_km2': float(table['area_km2'].iloc[row]), 'seed_k': float(table['muskingumK'].iloc[seed]),
                    **difference_row(good, bad, seed_scale),
                })
    return pd.DataFrame(rows)


def plot(results: pd.DataFrame) -> None:
    """Median and 90th percentile of the alternating and largest differences against distance below the seed."""
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True, layout='constrained')
    x = np.arange(len(DISTANCES_KM))
    ramp = dict(zip(results['scenario'].unique(), plotting.DT_RAMP[::-1], strict=False))
    for axis, column, title in ((axes[0], 'alternating', 'Alternating part'), (axes[1], 'max_difference', 'Largest')):
        for scenario, color in ramp.items():
            data = results[(results['scenario'] == scenario) & (results['dt'] == 3600)]
            q50 = data.groupby('distance_km')[column].median().reindex(DISTANCES_KM)
            q90 = data.groupby('distance_km')[column].quantile(0.9).reindex(DISTANCES_KM)
            axis.plot(x, 100 * q50, color=color, marker='o', markersize=3, label=f'{scenario}, median')
            axis.plot(x, 100 * q90, color=color, linestyle='--', linewidth=1, label=f'{scenario}, 90th percentile')
        axis.set_yscale('log')
        axis.set_xticks(x, [f'{d:g}' for d in DISTANCES_KM])
        axis.set_xlabel('Channel distance below the seed (km)')
        axis.set_title(f'{title} difference, Δt = 1 h', fontsize=9)
    axes[0].set_ylabel('% of the range of discharge')
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='outside lower center', ncol=2, fontsize=7)
    plotting.save(figure, 'single_defect_distance')
    return


def plot_long(results: pd.DataFrame) -> None:
    """
    Largest difference below one seed, as a share of the seed's own range of discharge, which removes dilution: a short
    reach left untreated against one long reach split into substeps, x held or x adjusted, at 1 h.
    """
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True, layout='constrained')
    x = np.arange(len(DISTANCES_KM))
    kinds = (('short', plotting.TREATMENT_COLORS['standard'], 'Short reach left untreated'),
             ('long-substeps', plotting.TREATMENT_COLORS['substeps'], 'Long reach in substeps'),
             ('long-substeps-xadj', plotting.TREATMENT_COLORS['substeps-xadj'], 'Long reach in substeps, x adjusted'))
    for axis, scenario in zip(axes, ('synthetic-burst', 'columbia-peak'), strict=True):
        for kind, color, label in kinds:
            data = results[(results['kind'] == kind) & (results['scenario'] == scenario) & (results['dt'] == 3600)]
            grouped = data.groupby('distance_km')['max_difference_of_seed']
            axis.plot(x, 100 * grouped.median().reindex(DISTANCES_KM), color=color, marker='o', markersize=3,
                      label=f'{label}, median')
            axis.plot(x, 100 * grouped.quantile(0.9).reindex(DISTANCES_KM), color=color, linestyle='--',
                      linewidth=1, label=f'{label}, 90th percentile')
        axis.set_yscale('log')
        axis.set_xticks(x, [f'{d:g}' for d in DISTANCES_KM])
        axis.set_xlabel('Channel distance below the seed (km)')
        axis.set_title(f'{scenario}, Δt = 1 h', fontsize=9)
    axes[0].set_ylabel("Largest difference\n(% of the seed's range)")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='outside lower center', ncol=2, fontsize=7)
    plotting.save(figure, 'single_defect_long_distance')
    return


if __name__ == '__main__':
    plotting.apply_style()
    columbia = hydrofabric.load()
    catchment_area = forcing.catchment_areas()
    frames = []
    for name in (*forcing.SYNTHETIC_PERTURBATIONS, *EVENTS):
        for step in STEPS:
            frames.append(measure('short', name, step, columbia, catchment_area))
            print('short', name, step, len(frames[-1]), flush=True)
        for kind in LONG_KINDS:
            for step in LONG_STEPS:
                frames.append(measure(kind, name, step, columbia, catchment_area))
                print(kind, name, step, len(frames[-1]), flush=True)
    single = pd.concat(frames, ignore_index=True)
    metrics.save_table(single, 'single_defect')
    plot(single[single['kind'] == 'short'])
    plot_long(single)
    columns = ['alternating', 'max_difference', 'max_difference_of_seed', 'peak_relative']
    summary = single.groupby(['kind', 'scenario', 'dt', 'distance_km'])[columns]
    print(summary.median().unstack('distance_km').round(5).to_string())
    print(summary.quantile(0.9).unstack('distance_km').round(5).to_string())
