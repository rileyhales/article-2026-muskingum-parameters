"""
Do negative coefficients make discharge oscillate between positive and negative values under real forcing? This finds
the reach whose ERA5 discharge at 1 h over the matrix period carried the largest share of its volume as negative
discharge, routes the
year of its lowest value again at 1 h standard, with each short-reach treatment, and as the reference, and plots the
days around that value. It also tabulates every reach whose negative discharge exceeds 0.1% of its volume.

``analyze`` routes the example and writes tables/era5_negative_reaches.csv, tables/era5_negative_flow_example.csv, and
tables/era5_negative_flow_series.csv (the routed year of the example reach under every cell), then draws. ``draw``
draws figure era5_negative_flow_example from those tables alone, so nothing is routed again to redraw it.

Run with the river-route environment:
    ../river-route/.venv/bin/python scripts/11_era5_negative_flow_example.py [analyze|draw]
"""

import argparse


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
    theory,
    treatments,
)

DT = 3600
SHOWN_HOURS = 48
CELLS = (('reference', config.REFERENCE[1]), ('standard', DT), ('inflate-k', DT), ('subcycles', DT), ('merge', DT))


def negative_reaches(table: pd.DataFrame) -> pd.DataFrame:
    """Every reach whose negative discharge at DT over the period is more than 0.1% of its volume, worst first."""
    cell = metrics.load_cell(config.ERA5_SCENARIO, 'standard', DT)
    stats = np.load(cell.directory / 'annual_stats.npz')
    share = stats['negative_volume'].sum(axis=0) / np.maximum(stats['volume'].sum(axis=0), 1.0)
    _, too_short = rr.Network(table).unstable_mask(float(DT))
    found = pd.DataFrame({
        'river': cell.river_ids, 'negative_share': share, 'too_short': too_short,
        'negative_hours': stats['negative_hours'].sum(axis=0).astype(int), 'lowest_m3s': stats['min'].min(axis=0),
        'worst_year': stats['years'][np.argmin(stats['min'], axis=0)], 'k_s': table['muskingumK'].to_numpy(),
        'length_m': table['Length'].to_numpy(), 'order': table['strahlerOrder'].to_numpy(),
        'area_km2': table['area_km2'].to_numpy(),
    })
    found = found[found['negative_share'] > 1e-3].sort_values('negative_share', ascending=False)
    if found.empty:
        raise ValueError('no reach carries more than 0.1% of its volume as negative discharge')
    return found.reset_index(drop=True)


class KeepRiver:
    """Keeps the hourly series of one river id from every chunk."""

    def __init__(self, river_ids: np.ndarray, river: int) -> None:
        rows = np.flatnonzero(river_ids == river)
        if rows.shape[0] != 1:
            raise ValueError(f'river {river} must be routed exactly once')
        self.row, self.values, self.dates = int(rows[0]), [], []

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        self.values.append(discharge[self.row].copy())
        self.dates.append(dates.copy())


def route_year(treatment: str, dt: int, table: pd.DataFrame, river: int, year: int) -> tuple | None:
    """
    One river's hourly dates and discharge through ``year``, routed from a dry network three months before it, or None
    when the treatment removes the river.
    """
    treated = treatments.treat(treatment, dt, table)
    if river not in treated.network.original_river_ids:
        return None
    keep = KeepRiver(treated.network.original_river_ids, river)
    files = forcing.era5_files(f'{year - 1}-10', f'{year}-12')
    engine.route(treated, dt, config.DT_RUNOFF, forcing.era5_chunks(files, table['riverId'].to_numpy()), keep, 8)
    dates, values = np.concatenate(keep.dates), np.concatenate(keep.values)
    if dates.shape != values.shape:
        raise ValueError('one value per date')
    in_year = dates >= np.datetime64(f'{year}-01-01')
    return dates[in_year], values[in_year]


def plot(example: pd.Series, series: pd.DataFrame) -> None:
    """The days around the lowest discharge of the example reach at DT, under every cell, from the tables."""
    standard = series[(series['treatment'] == 'standard') & (series['dt'] == DT)].sort_values('date')
    if standard.empty:
        raise ValueError(f'the standard network at {DT} s must be drawn')
    lowest = standard['date'].iloc[int(np.argmin(standard['discharge'].to_numpy()))]
    first, last = lowest - np.timedelta64(SHOWN_HOURS // 2, 'h'), lowest + np.timedelta64(SHOWN_HOURS // 2, 'h')
    window = series[(series['date'] >= first) & (series['date'] < last)]
    start = window['date'].min()
    figure, axis = plt.subplots(figsize=(plotting.WIDTH, 4.0))
    for treatment, _ in CELLS:
        data = window[window['treatment'] == treatment].sort_values('date')
        if data.empty:
            continue  # a treatment that removed the reach, or a step left out of the figures
        on_top = treatment == 'reference'
        axis.plot((data['date'] - start) / np.timedelta64(1, 'h'), data['discharge'],
                  color=plotting.TREATMENT_COLORS[treatment], linewidth=1.2, linestyle='--' if on_top else '-',
                  zorder=5 if on_top else 2, label=plotting.TREATMENT_LABELS[treatment] + ('' if on_top else ', 1 h'))
    axis.axhline(0, color=plotting.MUTED, linewidth=0.6)
    axis.set_title(f'Reach {int(example["river"])}: L = {example["length_m"]:.0f} m, k = {example["k_s"]:.0f} s, '
                   f'c3 = {example["c3"]:.2f}, c4 = {example["c4"]:.2f}')
    axis.set_xlabel(f'Hours from {str(start)[:13]}')
    axis.set_ylabel('Discharge (m³/s)')
    axis.legend(loc='upper right')
    plotting.save(figure, 'era5_negative_flow_example')
    return


def analyze() -> None:
    """Find the example reach, route its worst year under every cell, and write the tables the figure is drawn from."""
    columbia = hydrofabric.load()
    reaches = negative_reaches(columbia)
    metrics.save_table(reaches, 'era5_negative_reaches')
    worst = reaches.iloc[0]
    river, year = int(worst['river']), int(worst['worst_year'])
    routed = {cell: route_year(*cell, columbia, river, year) for cell in CELLS}
    routed = {cell: series for cell, series in routed.items() if series is not None}
    series = pd.concat([pd.DataFrame({'treatment': t, 'dt': dt, 'date': dates, 'discharge': values})
                        for (t, dt), (dates, values) in routed.items()], ignore_index=True)
    metrics.save_table(series, 'era5_negative_flow_series')
    _, _, c3, c4 = (float(c) for c in theory.coefficients(DT / worst['k_s'], 0.2))
    summary = {f'{t}_{dt}': {'min': float(v.min()), 'max': float(v.max()), 'negative_hours': int((v < 0).sum())}
               for (t, dt), (_, v) in routed.items()}
    example = pd.DataFrame([{**worst.to_dict(), 'year': year, 'c3': c3, 'c4': c4,
                             **{f'{cell}_{key}': value for cell, values in summary.items()
                                for key, value in values.items()}}])
    metrics.save_table(example, 'era5_negative_flow_example')
    print(reaches.head(15).to_string())
    print(example.T.to_string())
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', nargs='?', default='analyze', choices=('analyze', 'draw'))
    args = parser.parse_args()
    plotting.apply_style()
    if args.action == 'analyze':
        analyze()
    flows = metrics.shown(metrics.load_table('era5_negative_flow_series'))
    plot(metrics.load_table('era5_negative_flow_example').iloc[0], flows.assign(date=pd.to_datetime(flows['date'])))
