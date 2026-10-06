"""
Histogram of the routing time steps at which the reaches of each hydrofabric have non-negative Muskingum coefficients.
Every reach of every published hydrofabric of step 16 is given k = L / v with one celerity v and one x, so its window
of steps is 2kx <= dt <= 2k(1 - x). A bin of steps counts every reach whose window holds any step of the bin, as a share
of the reaches of the hydrofabric. At v = 1 m/s and x = 0.25, one figure is drawn per hydrofabric and one for the
reaches of all of them together; the reaches of all of them together are drawn again at every pair of v = 1 ± 0.5
m/s and x = 0.25 ± 0.15, in one figure of nine panels.

``analyze`` loads the hydrofabrics, writes tables/dt_window_histogram.csv with every hydrofabric at every pair, then
draws. ``draw`` draws from that table alone.

Run with the project environment:  uv run python scripts/18_dt_window_histogram.py
"""

import argparse
import importlib
import itertools

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from stability import config, plotting

hydrofabric_census = importlib.import_module('16_hydrofabric_census')

VELOCITY = 1.0  # m/s: the flood wave celerity of every reach in the figures of each hydrofabric, so k = L / v
X = 0.25  # the Muskingum x of every reach in those figures
VELOCITIES = (0.5, 1.0, 1.5)  # m/s: 1 ± 0.5, the celerities of the figures of all hydrofabrics together
XS = (0.10, 0.25, 0.40)  # 0.25 ± 0.15, and their values of x
DECADES = (0, 7)  # the bins span 10^0 to 10^7 s
BINS_PER_DECADE = 10  # evenly spaced in log
COMBINED = 'Combined'  # the reaches of every hydrofabric together
TABLE = 'dt_window_histogram'
MARKED_STEPS = (300, 900, 1800, 3600, 7200)  # s: the routing time steps every histogram marks
X_LIMITS = (1e1, 1e6)  # s: where the axis of every histogram starts and stops; the table holds the bins outside it too
BAR_COLOR = '#2a78d6'  # the blue of the first categorical slot of the style, for every histogram
# inches: the grid runs into the margins of the page as far as a float may, so the marked steps can be labeled above
# panels a third as wide as a figure, in a smaller size
GRID_WIDTH = 7.0
GRID_STEP_FONT_SIZE = 8


def bin_edges() -> np.ndarray:
    """The edges of the bins of routing time steps, s."""
    bins = (DECADES[1] - DECADES[0]) * BINS_PER_DECADE
    if bins < 1:
        raise ValueError('the bins must span at least one decade')
    edges = np.logspace(DECADES[0], DECADES[1], bins + 1)
    if np.any(np.diff(edges) <= 0):
        raise ValueError('bin edges must rise')
    return edges


def window(length_m: np.ndarray, velocity: float, x: float) -> tuple[np.ndarray, np.ndarray]:
    """The shortest and longest step with non-negative coefficients for every reach, 2kx and 2k(1 - x), s."""
    if not np.all(np.isfinite(length_m)) or np.any(length_m < 0):
        raise ValueError('lengths must be finite and non-negative')
    if velocity <= 0 or not 0 <= x <= 0.5:
        raise ValueError('the celerity must be positive and x between 0 and 0.5')
    k = length_m / velocity
    return 2 * k * x, 2 * k * (1 - x)


def valid_counts(length_m: np.ndarray, velocity: float, x: float, edges: np.ndarray) -> np.ndarray:
    """The reaches whose window holds any step of each bin [a, b), those with 2kx < b and 2k(1 - x) >= a."""
    shortest, longest = window(length_m, velocity, x)
    if np.any(longest > edges[-1]):
        raise ValueError(f'a window reaches past the last bin, {edges[-1]:g} s')
    # every reach with 2k(1 - x) < a also has 2kx < b, so the second count is part of the first
    counts = np.searchsorted(np.sort(shortest), edges[1:]) - np.searchsorted(np.sort(longest), edges[:-1])
    if np.any(counts < 0) or np.any(counts > length_m.shape[0]):
        raise ValueError('a bin cannot hold fewer than none or more than every reach')
    return counts


def hydrofabric_table(name: str, edges: np.ndarray) -> pd.DataFrame:
    """The reaches of one hydrofabric valid in each bin, at every pair of celerity and x."""
    length = hydrofabric_census.HYDROFABRICS[name]()['length_m'].to_numpy()
    parts, below = [], []
    for velocity, x in itertools.product(VELOCITIES, XS):
        parts.append(pd.DataFrame({
            'hydrofabric': name, 'velocity_m_s': velocity, 'x': x, 'dt_lo_s': edges[:-1], 'dt_hi_s': edges[1:],
            'reaches_valid': valid_counts(length, velocity, x, edges), 'reaches_total': length.shape[0],
        }))
        below.append(int(np.sum(window(length, velocity, x)[1] < edges[0])))
    print(f'{name}: {length.shape[0]:,} reaches, {min(below):,} to {max(below):,} valid only below {edges[0]:g} s')
    return pd.concat(parts, ignore_index=True)


def analyze() -> None:
    """Count the reaches valid in each bin of every published hydrofabric and of all of them, and write the table."""
    edges = bin_edges()
    table = pd.concat(
        [hydrofabric_table(name, edges) for name in hydrofabric_census.FIGURE_HYDROFABRICS], ignore_index=True,
    )
    keys = ['velocity_m_s', 'x', 'dt_lo_s', 'dt_hi_s']
    combined = table.groupby(keys, sort=False)[['reaches_valid', 'reaches_total']].sum().reset_index()
    if len(combined) != len(VELOCITIES) * len(XS) * (edges.shape[0] - 1):
        raise ValueError('every hydrofabric must have the same bins and pairs')
    table = pd.concat([table, combined.assign(hydrofabric=COMBINED)[table.columns]], ignore_index=True)
    config.TABLES.mkdir(parents=True, exist_ok=True)
    table.to_csv(config.TABLES / f'{TABLE}.csv', index=False)
    peaks = table.loc[table.groupby(['hydrofabric', 'velocity_m_s', 'x'], sort=False)['reaches_valid'].idxmax()]
    peaks = peaks.assign(valid_pct=100 * peaks['reaches_valid'] / peaks['reaches_total'])
    base = (peaks['velocity_m_s'] == VELOCITY) & (peaks['x'] == X)
    print(peaks[base | (peaks['hydrofabric'] == COMBINED)].to_string(index=False))
    return


def bars(axis: plt.Axes, rows: pd.DataFrame, top: float) -> None:
    """Draw the histogram of the rows on an axis, with a line at each marked step."""
    if rows.empty:
        raise ValueError('no bins to draw')
    lo, hi = rows['dt_lo_s'].to_numpy(), rows['dt_hi_s'].to_numpy()
    if np.any(hi <= lo):
        raise ValueError('every bin must have a positive width')
    share = 100 * rows['reaches_valid'].to_numpy() / rows['reaches_total'].to_numpy()
    if np.any(share > top):
        raise ValueError(f'a bar rises past {top:g}%')
    axis.bar(lo, share, width=hi - lo, align='edge', color=BAR_COLOR, edgecolor=plotting.SURFACE, linewidth=0.5)
    axis.set_axisbelow(True)
    axis.set_xscale('log')
    axis.set_xlim(*X_LIMITS)
    axis.set_ylim(0, top)
    for dt in MARKED_STEPS:
        axis.axvline(dt, color=plotting.MUTED, linewidth=1, linestyle='--', zorder=1.5)
    return


def label_steps(axis: plt.Axes, size: float) -> None:
    """Label the marked steps above an axis, each over its line."""
    if axis.get_xscale() != 'log' or size <= 0:
        raise ValueError('the steps are labeled above a log axis, in a positive size')
    # the marked steps are too close together on a log axis to be labeled inside it, so each is labeled above its line
    marks = axis.secondary_xaxis('top')
    marks.set_ticks(MARKED_STEPS, labels=[f'{dt // 60} m' for dt in MARKED_STEPS], rotation=90,
                    color=plotting.MUTED)
    marks.xaxis.set_minor_locator(ticker.NullLocator())
    marks.tick_params(length=0, pad=2, labelsize=size)  # the lines are the ticks
    return


def plot(rows: pd.DataFrame, title: str | None, name: str, top: float) -> None:
    """One histogram, with the marked steps labeled above it and the title if any, saved as the figure of the name."""
    figure, axis = plt.subplots()
    bars(axis, rows, top)
    label_steps(axis, plotting.LEGEND_FONT_SIZE)
    if title is not None:
        axis.set_title(title)
    axis.set_xlabel('Routing time step Δt (s)')
    axis.set_ylabel('Share of reaches (%)')
    plotting.save(figure, name)
    return


def plot_grid(rows: pd.DataFrame, top: float) -> None:
    """The histograms of every pair, a row per celerity and a column per x, on shared axes."""
    figure, axes = plt.subplots(len(VELOCITIES), len(XS), sharex=True, sharey=True, figsize=(GRID_WIDTH, 7.5))
    for row, velocity in enumerate(VELOCITIES):
        for column, x in enumerate(XS):
            pair = rows[(rows['velocity_m_s'] == velocity) & (rows['x'] == x)]
            if pair.empty:
                raise ValueError(f'no bins at v = {velocity:g} m/s, x = {x:g}')
            bars(axes[row, column], pair, top)
        axes[row, -1].annotate(f'v = {velocity:.2f} m s⁻¹', (1, 0.5), xycoords='axes fraction', xytext=(6, 0),
                               textcoords='offset points', rotation=270, ha='left', va='center')
    for column, x in enumerate(XS):
        label_steps(axes[0, column], GRID_STEP_FONT_SIZE)
        axes[0, column].set_title(f'x = {x:.2f}')
    axes[len(VELOCITIES) // 2, 0].set_ylabel('Share of reaches (%)')
    axes[-1, len(XS) // 2].set_xlabel('Routing time step Δt (s)')
    plotting.save(figure, f'{TABLE}_combined_grid')
    return


def share_top(rows: pd.DataFrame) -> float:
    """The top of a scale of shares that holds every bar of the rows, percent, to the next 10."""
    if rows.empty:
        raise ValueError('no rows')
    top = 10 * np.ceil(10 * (rows['reaches_valid'] / rows['reaches_total']).max())
    if not 0 < top <= 100:
        raise ValueError('shares must lie between 0 and 100%')
    return float(top)


def draw() -> None:
    """Draw every hydrofabric at the one pair, on one scale, and all of them together at every pair, on another."""
    table = pd.read_csv(config.TABLES / f'{TABLE}.csv')
    names = (*hydrofabric_census.FIGURE_HYDROFABRICS, COMBINED)
    pairs = set(itertools.product(VELOCITIES, XS))
    if set(table['hydrofabric']) != set(names) or set(zip(table['velocity_m_s'], table['x'])) != pairs:
        raise ValueError(f'the table must hold {", ".join(names)} at every pair of celerity and x')
    base = table[(table['velocity_m_s'] == VELOCITY) & (table['x'] == X)]
    top = share_top(base)
    for name in names:
        title = None if name == COMBINED else name  # the combined figure goes in the paper, whose caption names it
        plot(base[base['hydrofabric'] == name], title, f'{TABLE}_{name.lower().replace(" ", "-")}', top)
    combined = table[table['hydrofabric'] == COMBINED]
    plot_grid(combined, share_top(combined))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', nargs='?', default='analyze', choices=('analyze', 'draw'))
    args = parser.parse_args()
    plotting.apply_style()
    if args.action == 'analyze':
        analyze()
    draw()
