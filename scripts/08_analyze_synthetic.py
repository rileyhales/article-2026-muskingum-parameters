"""
Analyze the synthetic scenarios: for every cell, each river's peak, peak hour, dip, and negative discharge, and the
alternating (period of two hours) and largest parts of its departure from the reference. Errors are normalized by
the reference rise of each river, its peak less its steady baseflow, so a headwater and the outlet weigh the same.

``analyze`` writes tables/synthetic_summary.csv, tables/synthetic_by_class.csv, and tables/synthetic_examples.csv (the
hydrographs of the example rivers) from the routed cells, then draws. ``draw`` draws the synthetic figures from those
tables alone, showing the steps of config.DT_ROUTING, so a step left out needs nothing analyzed or routed again.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/08_analyze_synthetic.py [analyze|draw]
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import river_route as rr

from stability import config, forcing, metrics, plotting

SCENARIOS = tuple(forcing.SYNTHETIC_PERTURBATIONS)
OSCILLATION_THRESHOLD = 0.01  # an oscillation counts when its amplitude passes 1% of the river's rise
SHOWN = ('standard', 'substeps', 'substeps-xadj', 'subcycles', 'stabilized', 'inflate-k', 'merge')
# the panels of the example hydrographs: scenario, example river, the cells drawn with the reference, and title
EXAMPLE_PANELS = (
    ('synthetic-square', 'short', (('standard', 3600), ('subcycles', 3600), ('inflate-k', 3600), ('merge', 3600)),
     'River too short for 1 h'),
    ('synthetic-burst', 'dip', (('standard', 60), ('substeps', 60), ('substeps-xadj', 60)),
     'Long river ahead of a sharp rise'),
    ('synthetic-burst', 'outlet', (('standard', 3600), ('substeps', 60), ('substeps-xadj', 60), ('stabilized', 3600)),
     'Columbia at the Pacific'),
)


def cells_of(scenario: str, root=config.RESULTS, shown=SHOWN) -> list[tuple[str, int]]:
    """Every (treatment, dt) that has been run for a scenario under ``root``, the reference first."""
    found = [(config.REFERENCE[0], config.REFERENCE[1])]
    for treatment in shown:
        for dt in config.DT_ROUTING:
            if (config.run_dir(scenario, treatment, dt, root) / 'meta.json').exists():
                found.append((treatment, dt))
    if len(found) < 2:
        raise FileNotFoundError(f'no cells of {scenario} have been run')
    return found


def river_measures(series: np.ndarray) -> dict[str, np.ndarray]:
    """Per river measures of a series: peak, peak hour, oscillation amplitude, dip, zero and oscillation hours."""
    peak, peak_hour = metrics.peaks(series)
    measures = {
        'peak': peak, 'peak_hour': peak_hour, 'nyquist': metrics.nyquist_amplitude(series),
        'dip': metrics.dip_below_start(series, peak_hour), 'negative_hours': (np.asarray(series) < 0).sum(axis=1),
        'negative_volume': -np.minimum(np.asarray(series, dtype=np.float64), 0).sum(axis=1),
        'oscillation_hours': metrics.oscillation_hours(series), 'start': np.asarray(series[:, 0], dtype=np.float64),
    }
    if not np.all(np.isfinite(measures['peak'])):
        raise ValueError('every peak must be finite')
    return measures


def compare(
    scenario: str, network: rr.Network, area_km2: np.ndarray, root=config.RESULTS, shown=SHOWN
) -> pd.DataFrame:
    """One row per river and cell under ``root``: its errors against the reference and its class at the cell's dt."""
    reference = metrics.load_cell(scenario, *config.REFERENCE, root)
    reference_series = np.asarray(reference.array('series'), dtype=np.float64)
    ref = river_measures(reference_series)
    rise = ref['peak'] - ref['start']
    if np.any(rise <= 0):
        raise ValueError('every river must rise in the reference')
    frames = []
    for treatment, dt in cells_of(scenario, root, shown):
        cell = metrics.load_cell(scenario, treatment, dt, root)
        rows = metrics.align(cell.river_ids, reference.river_ids)
        kept = rows >= 0
        series = np.asarray(cell.array('series'))[rows[kept]]
        measured = river_measures(series)
        error = series - reference_series[kept]
        too_long, too_short = network.unstable_mask(float(dt))
        river_class = np.where(too_long, 'too long', np.where(too_short, 'too short', 'valid'))[kept]
        frames.append(pd.DataFrame({
            'scenario': scenario, 'treatment': treatment, 'dt': dt, 'river': reference.river_ids[kept],
            'class': river_class, 'area_km2': area_km2[kept],
            'peak_error': (measured['peak'] - ref['peak'][kept]) / rise[kept],
            'timing_error_h': measured['peak_hour'] - ref['peak_hour'][kept],
            'oscillation': metrics.nyquist_amplitude(error) / rise[kept],
            'max_error': np.abs(error).max(axis=1) / rise[kept], 'dip': measured['dip'] / rise[kept],
            'volume_error': error.sum(axis=1) / reference_series[kept].sum(axis=1),
            'negative_hours': measured['negative_hours'],
            'negative_share': measured['negative_volume'] / np.maximum(series.sum(axis=1), 1e-6),
            'oscillation_hours': measured['oscillation_hours'],
        }))
    return pd.concat(frames, ignore_index=True)


def summarize(rivers: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Distribution of the errors over the rivers of each group."""
    if rivers.empty:
        raise ValueError('nothing to summarize')
    grouped = rivers.groupby(keys, sort=False)
    summary = grouped.agg(
        rivers=('peak_error', 'size'),
        peak_error_median=('peak_error', 'median'),
        peak_error_p05=('peak_error', lambda s: s.quantile(0.05)),
        peak_error_p95=('peak_error', lambda s: s.quantile(0.95)),
        peak_error_abs_mean=('peak_error', lambda s: s.abs().mean()),
        timing_error_share=('timing_error_h', lambda s: (s != 0).mean()),
        timing_error_median_h=('timing_error_h', 'median'),
        timing_error_mean_h=('timing_error_h', 'mean'),
        timing_error_abs_mean_h=('timing_error_h', lambda s: s.abs().mean()),
        timing_error_abs_p95_h=('timing_error_h', lambda s: s.abs().quantile(0.95)),
        volume_error_abs_mean=('volume_error', lambda s: s.abs().mean()),
        volume_error_abs_max=('volume_error', lambda s: s.abs().max()),
        oscillating_share=('oscillation', lambda s: (s > OSCILLATION_THRESHOLD).mean()),
        oscillation_p99=('oscillation', lambda s: s.quantile(0.99)),
        oscillation_max=('oscillation', 'max'),
        max_error_median=('max_error', 'median'),
        max_error_p95=('max_error', lambda s: s.quantile(0.95)),
        dip_share=('dip', lambda s: (s > 0.01).mean()),
        dip_max=('dip', 'max'),
        rivers_with_negative_hours=('negative_hours', lambda s: (s > 0).sum()),
        negative_hours=('negative_hours', 'sum'),
        negative_share_max=('negative_share', 'max'),
        rivers_oscillating_hours=('oscillation_hours', lambda s: (s > 0).sum()),
    ).reset_index()
    if summary['rivers'].sum() != len(rivers):
        raise ValueError('every river must be summarized once')
    return summary


def plot_errors(summary: pd.DataFrame) -> None:
    """Median and spread of the peak error, and the share of rivers whose peak hour moved, against dt."""
    figure, axes = plt.subplots(3, 2, figsize=(plotting.WIDTH, 8.5), sharex=True)
    steps = np.asarray(config.DT_ROUTING)
    for column, scenario in enumerate(SCENARIOS):
        rows = summary[summary['scenario'] == scenario]
        for treatment in SHOWN:
            data = rows[rows['treatment'] == treatment].sort_values('dt')
            style = {'color': plotting.TREATMENT_COLORS[treatment], 'marker': 'o', 'markersize': 3,
                     'label': plotting.TREATMENT_LABELS[treatment]}
            positions = plotting.step_positions(data['dt'], steps)
            axes[0, column].plot(positions, 100 * data['peak_error_median'], **style)
            axes[1, column].plot(positions, 100 * data['peak_error_abs_mean'], **style)
            axes[2, column].plot(positions, 100 * data['timing_error_share'], **style)
        axes[0, column].set_title(scenario.replace('synthetic-', 'Synthetic ').title())
        axes[0, column].axhline(0, color=plotting.MUTED, linewidth=0.6)
        plotting.step_ticks(axes[2, column], steps)
        axes[2, column].set_xlabel('Routing time step')
    axes[0, 0].set_ylabel('Median peak error\n(% of rise)')
    axes[1, 0].set_ylabel('Mean |peak error|\n(% of rise)')
    axes[2, 0].set_ylabel('Rivers with peak hour\nmoved (%)')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2)
    plotting.save(figure, 'synthetic_peak_errors')
    return


def plot_artifacts(summary: pd.DataFrame) -> None:
    """Share of rivers oscillating, and with a dip, against dt, for the square wave and the burst."""
    figure, grid = plt.subplots(2, 2, figsize=(plotting.WIDTH, 7.0))
    axes, key = grid.flat[:3], grid.flat[3]  # three panels, and the fourth cell holds the legend
    steps = np.asarray(config.DT_ROUTING)
    panels = (('synthetic-square', 'oscillating_share', '(a) Alternating, square', 'Rivers > 1% of rise (%)'),
              ('synthetic-burst', 'dip_share', '(b) Dipping, burst', 'Rivers > 1% of rise (%)'),
              ('synthetic-burst', 'rivers_with_negative_hours', '(c) Negative, burst', 'Rivers'))
    for axis, (scenario, column, title, label) in zip(axes, panels, strict=True):
        rows = summary[summary['scenario'] == scenario]
        scale = 1 if column == 'rivers_with_negative_hours' else 100
        for treatment in SHOWN:
            data = rows[rows['treatment'] == treatment].sort_values('dt')
            axis.plot(plotting.step_positions(data['dt'], steps), scale * data[column], marker='o', markersize=3,
                      color=plotting.TREATMENT_COLORS[treatment], label=plotting.TREATMENT_LABELS[treatment])
        reference = rows[rows['treatment'] == 'reference'][column].iloc[0] * scale
        axis.axhline(reference, color=plotting.INK, linewidth=1, linestyle='--', label='Reference')
        plotting.step_ticks(axis, steps)
        axis.set_title(title)
        axis.set_ylabel(label)
    key.axis('off')
    key.legend(*axes[0].get_legend_handles_labels(), loc='center')
    plotting.save(figure, 'synthetic_artifacts')
    return


def example_rivers(rivers: pd.DataFrame) -> dict[str, int]:
    """A river too short at one hour that oscillates most, and a long river that dips most in the reference."""
    square = rivers[(rivers['scenario'] == 'synthetic-square') & (rivers['treatment'] == 'standard')
                    & (rivers['dt'] == 3600) & (rivers['class'] == 'too short') & (rivers['area_km2'] > 100)]
    burst = rivers[(rivers['scenario'] == 'synthetic-burst') & (rivers['treatment'] == 'reference')]
    if square.empty or burst.empty:
        raise ValueError('the example rivers must exist')
    return {
        'short': int(square.sort_values('oscillation')['river'].iloc[-1]),
        'dip': int(burst.sort_values('dip')['river'].iloc[-1]),
        'outlet': config.COLUMBIA_OUTLET,
    }


def example_series(examples: dict[str, int]) -> pd.DataFrame:
    """The whole hourly series of each example river under the reference and the cells of its panel."""
    frames = []
    for scenario, key, cells, _ in EXAMPLE_PANELS:
        for treatment, dt in (config.REFERENCE, *cells):
            cell = metrics.load_cell(scenario, treatment, dt)
            row = np.flatnonzero(cell.river_ids == examples[key])
            if row.shape[0] != 1:
                continue  # the treatment removed the river
            flow = np.asarray(cell.array('series')[row[0]], dtype=np.float64)
            frames.append(pd.DataFrame({
                'panel': key, 'scenario': scenario, 'river': examples[key], 'treatment': treatment, 'dt': dt,
                'hour': np.arange(flow.shape[0]), 'discharge': flow,
            }))
    if not frames:
        raise ValueError('no example series')
    return pd.concat(frames, ignore_index=True)


def plot_examples(series: pd.DataFrame) -> None:
    """Hydrographs of the example rivers under the treatments that matter to each, from the examples table."""
    figure, axes = plt.subplots(1, 3, figsize=(plotting.WIDTH, 4.0))
    for axis, (_, key, cells, title) in zip(axes, EXAMPLE_PANELS, strict=True):
        hours = (18, 18 + 48) if key != 'outlet' else (0, forcing.SYNTHETIC['hours'])
        panel = series[(series['panel'] == key) & series['hour'].between(hours[0], hours[1] - 1)]
        for treatment, dt in (('reference', None), *cells):
            data = panel[panel['treatment'] == treatment]
            if dt is not None:
                data = data[data['dt'] == dt]
            if data.empty:
                continue  # a step left out of the figures, or a river the treatment removed
            label = plotting.TREATMENT_LABELS[treatment] + ('' if dt is None else f', {dt} s')
            on_top = treatment == 'reference'
            axis.plot(data['hour'] / (24 if key == 'outlet' else 1), data['discharge'], label=label, linewidth=1.2,
                      color=plotting.TREATMENT_COLORS[treatment], zorder=5 if on_top else 2,
                      linestyle='--' if on_top else '-')
        axis.set_title(title)
        axis.set_xlabel('Day' if key == 'outlet' else 'Hour')
        axis.legend(loc='upper right')
    axes[0].set_ylabel('Discharge (m³/s)')
    plotting.save(figure, 'synthetic_examples')
    return


def analyze() -> None:
    """Measure every routed cell against the reference and write the tables the figures are drawn from."""
    network = rr.Network(config.NETWORK_FILE)
    metadata = pd.read_parquet(config.METADATA_FILE, columns=['riverId', 'DSContArea'])
    area = metadata['DSContArea'].to_numpy() / 1e6
    rivers = pd.concat([compare(s, network, area) for s in SCENARIOS], ignore_index=True)
    summary = summarize(rivers, ['scenario', 'treatment', 'dt'])
    metrics.save_table(summary, 'synthetic_summary')
    metrics.save_table(summarize(rivers, ['scenario', 'treatment', 'dt', 'class']), 'synthetic_by_class')
    examples = example_rivers(rivers)
    metrics.save_table(example_series(examples), 'synthetic_examples')
    print(examples)
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(summary.round(4).to_string(index=False))
    return


def draw() -> None:
    """Draw the synthetic figures from the analysis tables, at the steps of config.DT_ROUTING."""
    summary = metrics.shown(metrics.load_table('synthetic_summary'))
    plot_errors(summary)
    plot_artifacts(summary)
    if metrics.has_table('synthetic_examples'):
        plot_examples(metrics.shown(metrics.load_table('synthetic_examples')))
    else:
        print('synthetic_examples not drawn: tables/synthetic_examples.csv is written by analyze')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', nargs='?', default='analyze', choices=('analyze', 'draw'))
    args = parser.parse_args()
    plotting.apply_style()
    if args.action == 'analyze':
        analyze()
    draw()
