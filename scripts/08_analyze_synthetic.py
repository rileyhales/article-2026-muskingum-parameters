"""
Analyze the synthetic scenarios: for every cell, each river's peak, peak hour, dip, and negative discharge, and the
alternating (period of two hours) and largest parts of its departure from the reference. Errors are normalized by
the reference rise of each river, its peak less its steady baseflow, so a headwater and the outlet weigh the same.

Writes tables/synthetic_summary.csv, tables/synthetic_by_class.csv, and the synthetic figures.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/08_analyze_synthetic.py
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import river_route as rr

from stability import config, forcing, metrics, plotting

SCENARIOS = tuple(forcing.SYNTHETIC_PERTURBATIONS)
OSCILLATION_THRESHOLD = 0.01  # an oscillation counts when its amplitude passes 1% of the river's rise
SHOWN = ('standard', 'substeps', 'substeps-xadj', 'subcycles', 'stabilized', 'inflate-k', 'merge')


def cells_of(scenario: str) -> list[tuple[str, int]]:
    """Every (treatment, dt) that has been run for a scenario, the reference first."""
    found = [(config.REFERENCE[0], config.REFERENCE[1])]
    for treatment in SHOWN:
        for dt in config.DT_ROUTING:
            if (config.run_dir(scenario, treatment, dt) / 'meta.json').exists():
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


def compare(scenario: str, network: rr.Network, area_km2: np.ndarray) -> pd.DataFrame:
    """One row per river and cell: its errors against the reference and its class at the cell's dt."""
    reference = metrics.load_cell(scenario, *config.REFERENCE)
    reference_series = np.asarray(reference.array('series'), dtype=np.float64)
    ref = river_measures(reference_series)
    rise = ref['peak'] - ref['start']
    if np.any(rise <= 0):
        raise ValueError('every river must rise in the reference')
    frames = []
    for treatment, dt in cells_of(scenario):
        cell = metrics.load_cell(scenario, treatment, dt)
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
    figure, axes = plt.subplots(3, 2, figsize=(7.2, 7.6), sharex=True, layout='constrained')
    steps = np.asarray(config.DT_ROUTING)
    for column, scenario in enumerate(SCENARIOS):
        rows = summary[summary['scenario'] == scenario]
        for treatment in SHOWN:
            data = rows[rows['treatment'] == treatment].sort_values('dt')
            style = {'color': plotting.TREATMENT_COLORS[treatment], 'marker': 'o', 'markersize': 3,
                     'label': plotting.TREATMENT_LABELS[treatment]}
            axes[0, column].plot(data['dt'], 100 * data['peak_error_median'], **style)
            axes[1, column].plot(data['dt'], 100 * data['peak_error_abs_mean'], **style)
            axes[2, column].plot(data['dt'], 100 * data['timing_error_share'], **style)
        axes[0, column].set_title(scenario.replace('synthetic-', 'Synthetic ').title())
        axes[0, column].axhline(0, color=plotting.MUTED, linewidth=0.6)
        axes[2, column].set_xscale('log')
        axes[2, column].set_xticks(steps, [f'{s // 60}m' if s >= 60 else f'{s}s' for s in steps], rotation=45)
        axes[2, column].set_xlabel('Routing time step')
    axes[0, 0].set_ylabel('Median peak error\n(% of rise)')
    axes[1, 0].set_ylabel('Mean |peak error|\n(% of rise)')
    axes[2, 0].set_ylabel('Rivers with peak hour\nmoved (%)')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='outside lower center', ncol=4)
    plotting.save(figure, 'synthetic_peak_errors')
    return


def plot_artifacts(summary: pd.DataFrame) -> None:
    """Share of rivers oscillating, and with a dip, against dt, for the square wave and the burst."""
    figure, axes = plt.subplots(1, 3, figsize=(7.2, 3.6), layout='constrained')
    steps = np.asarray(config.DT_ROUTING)
    panels = (('synthetic-square', 'oscillating_share', 'Rivers with alternating error\n> 1% of rise, square (%)'),
              ('synthetic-burst', 'dip_share', 'Rivers dipping > 1%\nof rise, burst (%)'),
              ('synthetic-burst', 'rivers_with_negative_hours', 'Rivers with negative\ndischarge, burst'))
    for axis, (scenario, column, label) in zip(axes, panels, strict=True):
        rows = summary[summary['scenario'] == scenario]
        scale = 1 if column == 'rivers_with_negative_hours' else 100
        for treatment in SHOWN:
            data = rows[rows['treatment'] == treatment].sort_values('dt')
            axis.plot(data['dt'], scale * data[column], color=plotting.TREATMENT_COLORS[treatment], marker='o',
                      markersize=3, label=plotting.TREATMENT_LABELS[treatment])
        reference = rows[rows['treatment'] == 'reference'][column].iloc[0] * scale
        axis.axhline(reference, color=plotting.INK, linewidth=1, linestyle='--', label='Reference')
        axis.set_xscale('log')
        axis.set_xticks(steps, [f'{s // 60}m' if s >= 60 else f'{s}s' for s in steps], rotation=45)
        axis.set_ylabel(label)
        axis.set_xlabel('Routing time step')
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='outside lower center', ncol=4)
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


def plot_examples(examples: dict[str, int]) -> None:
    """Hydrographs of the example rivers under the treatments that matter to each."""
    panels = (
        ('synthetic-square', 'short', (('reference', 30), ('standard', 3600), ('subcycles', 3600),
                                       ('inflate-k', 3600), ('merge', 3600)), 'River too short for 1 h'),
        ('synthetic-burst', 'dip', (('reference', 30), ('standard', 60), ('substeps', 60), ('substeps-xadj', 60)),
         'Long river ahead of a sharp rise'),
        ('synthetic-burst', 'outlet', (('reference', 30), ('standard', 3600), ('substeps', 60),
                                       ('substeps-xadj', 60), ('stabilized', 3600)), 'Columbia at the Pacific'),
    )
    figure, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), gridspec_kw={'wspace': 0.4})
    for axis, (scenario, key, cells, title) in zip(axes, panels, strict=True):
        hours = slice(18, 18 + 48) if key != 'outlet' else slice(0, forcing.SYNTHETIC['hours'])
        for treatment, dt in cells:
            cell = metrics.load_cell(scenario, treatment, dt)
            row = np.flatnonzero(cell.river_ids == examples[key])
            if row.shape[0] != 1:
                continue
            flow = np.asarray(cell.array('series')[row[0], hours])
            time = np.arange(hours.start, hours.stop) / (24 if key == 'outlet' else 1)
            label = plotting.TREATMENT_LABELS[treatment] + ('' if treatment == 'reference' else f', {dt} s')
            on_top = treatment == 'reference'
            axis.plot(time, flow, color=plotting.TREATMENT_COLORS[treatment], label=label, linewidth=1.2,
                      zorder=5 if on_top else 2, linestyle='--' if on_top else '-')
        axis.set_title(title, fontsize=9)
        axis.set_xlabel('Day' if key == 'outlet' else 'Hour')
        axis.legend(loc='upper right', fontsize=6)
    axes[0].set_ylabel('Discharge (m³/s)')
    plotting.save(figure, 'synthetic_examples')
    return


if __name__ == '__main__':
    plotting.apply_style()
    network = rr.Network(config.NETWORK_FILE)
    metadata = pd.read_parquet(config.METADATA_FILE, columns=['riverId', 'DSContArea'])
    area = metadata['DSContArea'].to_numpy() / 1e6
    rivers = pd.concat([compare(s, network, area) for s in SCENARIOS], ignore_index=True)
    summary = summarize(rivers, ['scenario', 'treatment', 'dt'])
    metrics.save_table(summary, 'synthetic_summary')
    metrics.save_table(summarize(rivers, ['scenario', 'treatment', 'dt', 'class']), 'synthetic_by_class')
    plot_errors(summary)
    plot_artifacts(summary)
    examples = example_rivers(rivers)
    plot_examples(examples)
    print(examples)
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(summary.round(4).to_string(index=False))
