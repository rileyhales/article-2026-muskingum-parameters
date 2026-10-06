"""
Analyze the ERA5 cells of the matrix period against the reference:

- annual statistics of every river: annual peak and its hour, volume, hours and volume of negative discharge,
  oscillation hours, and total variation, as errors against the reference, summarized per cell and per class of river
  at the cell's dt
- event windows of every river: hours of negative discharge
- event windows of every river: peak, peak hour, largest and alternating departures from the reference
- profiles of the event peak error along the main stems of the major tributaries

``analyze`` writes tables era5_annual_summary, era5_annual_by_class, era5_events_summary, era5_events_by_class,
era5_profiles, era5_cost, era5_negative, and era5_reference_artifacts from the routed cells, then draws. ``draw`` draws
the era5 figures from those tables alone, at the steps of config.DT_ROUTING, so a step left out needs nothing analyzed
or routed again.

Run with the project environment:  uv run python scripts/10_analyze_era5.py [analyze|draw]
"""

import argparse
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import river_route as rr

from stability import config, metrics, plotting

SCENARIO = config.ERA5_SCENARIO
TREATMENTS = ('standard', 'substeps', 'substeps-xadj', 'subcycles', 'stabilized', 'inflate-k', 'merge')
SAME_EVENT_HOURS = 72  # an annual peak more than this far from the reference peak is a different event
MIN_PEAK = 0.1  # m3/s; rivers whose reference peak is smaller are left out of relative errors
# the cells, the stems, and the event of the figure of peak error along the main stems
PROFILE_CELLS = (('standard', 3600), ('substeps', 60), ('substeps-xadj', 60), ('subcycles', 3600), ('substeps', 3600))
PROFILE_STEMS = ('Columbia', 'Snake', 'Willamette')
PROFILE_EVENT = 'columbia-peak'


def finished_cells() -> list[tuple[str, int]]:
    """Every (treatment, dt) of the ERA5 scenario whose run has finished."""
    cells = [(t, dt) for t in TREATMENTS for dt in config.DT_ROUTING
             if (config.run_dir(SCENARIO, t, dt) / 'meta.json').exists()]
    if not cells:
        raise FileNotFoundError('no ERA5 cells have finished')
    return cells


def river_classes(network: rr.Network, dt: int) -> np.ndarray:
    """The class of every river at dt: too long, too short, or valid."""
    too_long, too_short = network.unstable_mask(float(dt))
    classes = np.where(too_long, 'too long', np.where(too_short, 'too short', 'valid'))
    if classes.shape[0] != network.size:
        raise ValueError('one class per river')
    return classes


def annual_errors(cell: metrics.Cell, reference: dict, reference_ids: np.ndarray, classes: np.ndarray) -> pd.DataFrame:
    """One row per river and year: the errors of a cell's annual statistics against the reference."""
    stats = np.load(cell.directory / 'annual_stats.npz')
    rows = metrics.align(cell.river_ids, reference_ids)
    kept = rows >= 0
    take = rows[kept]
    n_years = stats['max'].shape[0]
    ref_max = reference['max'][:, kept]
    usable = ref_max >= MIN_PEAK
    frame = pd.DataFrame({
        'treatment': cell.treatment, 'dt': cell.dt,
        'year': np.repeat(stats['years'], kept.sum()),
        'river': np.tile(reference_ids[kept], n_years),
        'class': np.tile(classes[kept], n_years),
        'peak_error': np.where(usable, (stats['max'][:, take] - ref_max) / np.maximum(ref_max, MIN_PEAK), np.nan)
        .ravel(),
        'timing_error_h': (stats['argmax'][:, take] - reference['argmax'][:, kept]).ravel(),
        'volume_error': ((stats['volume'][:, take] - reference['volume'][:, kept])
                         / np.maximum(reference['volume'][:, kept], 1.0)).ravel(),
        'negative_hours': stats['negative_hours'][:, take].ravel(),
        'negative_share': (stats['negative_volume'][:, take] / np.maximum(stats['volume'][:, take], 1.0)).ravel(),
        'excess_negative_hours': (stats['negative_hours'][:, take] - reference['negative_hours'][:, kept]).ravel(),
        'oscillation_hours': stats['oscillation_hours'][:, take].ravel(),
        'excess_oscillation_hours': (stats['oscillation_hours'][:, take]
                                     - reference['oscillation_hours'][:, kept]).ravel(),
        'variation_ratio': (stats['total_variation'][:, take]
                            / np.maximum(reference['total_variation'][:, kept], 1e-6)).ravel(),
    })
    if len(frame) != n_years * kept.sum():
        raise ValueError('one row per river and year')
    return frame


def summarize_annual(frame: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Distribution of the annual errors over the river-years of each group."""
    same = frame['timing_error_h'].abs() <= SAME_EVENT_HOURS
    frame = frame.assign(same_event=same, timing_same=frame['timing_error_h'].where(same))
    summary = frame.groupby(keys, sort=False).agg(
        river_years=('peak_error', 'size'),
        peak_error_median=('peak_error', 'median'),
        peak_error_abs_mean=('peak_error', lambda s: s.abs().mean()),
        peak_error_p05=('peak_error', lambda s: s.quantile(0.05)),
        peak_error_p95=('peak_error', lambda s: s.quantile(0.95)),
        same_event_share=('same_event', 'mean'),
        timing_moved_share=('timing_same', lambda s: (s.dropna() != 0).mean()),
        timing_mean_h=('timing_same', lambda s: s.dropna().mean()),
        timing_abs_mean_h=('timing_same', lambda s: s.dropna().abs().mean()),
        timing_abs_p95_h=('timing_same', lambda s: s.dropna().abs().quantile(0.95)),
        volume_error_abs_mean=('volume_error', lambda s: s.abs().mean()),
        volume_error_abs_max=('volume_error', lambda s: s.abs().max()),
        volume_error_mean=('volume_error', 'mean'),
        river_years_negative=('negative_hours', lambda s: (s > 0).sum()),
        negative_hours=('negative_hours', 'sum'),
        excess_negative_hours=('excess_negative_hours', 'sum'),
        negative_share_max=('negative_share', 'max'),
        river_years_oscillating=('oscillation_hours', lambda s: (s > 0).sum()),
        oscillation_hours=('oscillation_hours', 'sum'),
        excess_oscillation_hours=('excess_oscillation_hours', 'sum'),
        variation_ratio_median=('variation_ratio', 'median'),
        variation_ratio_p99=('variation_ratio', lambda s: s.quantile(0.99)),
    ).reset_index()
    if summary['river_years'].sum() != len(frame):
        raise ValueError('every river-year must be summarized once')
    return summary


def event_errors(cell: metrics.Cell, reference: metrics.Cell, classes: np.ndarray, event: str) -> pd.DataFrame:
    """One row per river: the peak, timing, largest, and alternating errors of an event window."""
    rows = metrics.align(cell.river_ids, reference.river_ids)
    kept = rows >= 0
    good = np.asarray(reference.array(f'events/{event}'), dtype=np.float64)[kept]
    bad = np.asarray(cell.array(f'events/{event}'), dtype=np.float64)[rows[kept]]
    scale = np.maximum(good.max(axis=1) - good.min(axis=1), 1e-3)
    ref_peak, ref_hour = metrics.peaks(good)
    peak, hour = metrics.peaks(bad)
    difference = bad - good
    frame = pd.DataFrame({
        'treatment': cell.treatment, 'dt': cell.dt, 'event': event, 'river': reference.river_ids[kept],
        'class': classes[kept], 'peak_error': np.where(ref_peak >= MIN_PEAK, (peak - ref_peak) / ref_peak, np.nan),
        'timing_error_h': hour - ref_hour, 'max_error': np.abs(difference).max(axis=1) / scale,
        'alternating': metrics.nyquist_amplitude(difference) / scale,
        'negative_hours': (bad < 0).sum(axis=1),
        'reference_negative_hours': (good < 0).sum(axis=1),
        'most_negative': np.minimum(bad.min(axis=1), 0) / np.maximum(ref_peak, MIN_PEAK),
    })
    if frame['max_error'].isna().any():
        raise ValueError('every error must be finite')
    return frame


def summarize_events(frame: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Distribution of the event errors over the rivers of each group."""
    summary = frame.groupby(keys, sort=False).agg(
        rivers=('peak_error', 'size'),
        peak_error_median=('peak_error', 'median'),
        peak_error_abs_mean=('peak_error', lambda s: s.abs().mean()),
        peak_error_p95=('peak_error', lambda s: s.abs().quantile(0.95)),
        timing_moved_share=('timing_error_h', lambda s: (s != 0).mean()),
        timing_mean_h=('timing_error_h', 'mean'),
        timing_abs_mean_h=('timing_error_h', lambda s: s.abs().mean()),
        timing_abs_p95_h=('timing_error_h', lambda s: s.abs().quantile(0.95)),
        max_error_median=('max_error', 'median'),
        max_error_p95=('max_error', lambda s: s.quantile(0.95)),
        alternating_p99=('alternating', lambda s: s.quantile(0.99)),
        alternating_share=('alternating', lambda s: (s > 0.01).mean()),
        rivers_negative=('negative_hours', lambda s: (s > 0).sum()),
        negative_hours=('negative_hours', 'sum'),
        reference_negative_hours=('reference_negative_hours', 'sum'),
        most_negative=('most_negative', 'min'),
    ).reset_index()
    if summary['rivers'].sum() != len(frame):
        raise ValueError('every river must be summarized once')
    return summary


def stem_profiles(cells: list[tuple[str, int]], reference: metrics.Cell, event: str) -> pd.DataFrame:
    """The event peak error at every observation river along the main stems, for each cell."""
    observations = pd.read_parquet(config.INPUTS / 'observation_rivers.parquet')
    position = metrics.align(reference.river_ids, observations['riverId'].to_numpy())
    good = np.asarray(reference.array(f'events/{event}'), dtype=np.float64)[position]
    ref_peak, ref_hour = metrics.peaks(good)
    frames = []
    for treatment, dt in cells:
        cell = metrics.load_cell(SCENARIO, treatment, dt)
        rows = metrics.align(cell.river_ids, observations['riverId'].to_numpy())
        kept = rows >= 0
        peak, hour = metrics.peaks(np.asarray(cell.array(f'events/{event}'), dtype=np.float64)[rows[kept]])
        frames.append(observations[kept].assign(
            treatment=treatment, dt=dt, event=event, peak_error=(peak - ref_peak[kept]) / ref_peak[kept],
            timing_error_h=hour - ref_hour[kept]))
    return pd.concat(frames, ignore_index=True)


def negative_table(cells: list[tuple[str, int]], network: rr.Network) -> pd.DataFrame:
    """
    Negative discharge of every cell over the period: its share of the total volume, the reaches whose negative
    discharge exceeds 0.1% and 1% of their volume, the largest such share, and the lowest discharge.
    """
    rows = []
    for treatment, dt in [config.REFERENCE, *cells]:
        stats = np.load(metrics.load_cell(SCENARIO, treatment, dt).directory / 'annual_stats.npz')
        negative, volume = stats['negative_volume'].sum(axis=0), stats['volume'].sum(axis=0)
        share = negative / np.maximum(volume, 1.0)
        rows.append({
            'treatment': treatment, 'dt': dt, 'negative_volume_share': float(negative.sum() / volume.sum()),
            'reaches_over_0.1pct': int((share > 1e-3).sum()), 'reaches_over_1pct': int((share > 1e-2).sum()),
            'largest_share': float(share.max()), 'lowest_m3s': float(stats['min'].min()),
            'too_short_rivers': int(network.unstable_mask(float(dt))[1].sum()),
        })
    table = pd.DataFrame(rows)
    if (table['largest_share'] < 0).any():
        raise ValueError('a share of negative volume cannot be negative')
    return table


def cost_table(cells: list[tuple[str, int]]) -> pd.DataFrame:
    """Reaches, subcycles, reach-steps per simulated hour, and wall time of every cell, from its meta.json."""
    rows = []
    for treatment, dt in [config.REFERENCE, *cells]:
        meta = metrics.load_cell(SCENARIO, treatment, dt).meta
        rows.append({
            'treatment': treatment, 'dt': dt, 'reaches': meta['reaches'], 'subcycled_rivers': meta['subcycled_rivers'],
            'reach_steps_per_hour': meta['reach_steps_per_routing_step'] * config.DT_RUNOFF // dt,
            'wall_seconds': meta['wall_seconds'], 'rivers_removed': meta['rivers_removed'],
        })
    table = pd.DataFrame(rows)
    standard_hour = table[(table['treatment'] == 'standard') & (table['dt'] == 3600)]['reach_steps_per_hour']
    if standard_hour.shape[0] == 1:
        table['work_vs_standard_1h'] = table['reach_steps_per_hour'] / float(standard_hour.iloc[0])
    if table.empty:
        raise ValueError('no cells')
    return table


def plot_annual(summary: pd.DataFrame) -> None:
    """Annual peak error, peak timing, negative hours, and oscillation hours against dt, per treatment."""
    figure, axes = plt.subplots(2, 2, figsize=(plotting.WIDTH, 7.0), sharex=True)
    steps = np.asarray(config.DT_ROUTING)
    panels = ((axes[0, 0], 'peak_error_abs_mean', 100, 'Mean |annual peak error| (%)'),
              (axes[0, 1], 'timing_moved_share', 100, 'Annual peaks whose hour\nmoved (%)'),
              (axes[1, 0], 'excess_negative_hours', 1, 'Hours of negative discharge\nbeyond the reference'),
              (axes[1, 1], 'excess_oscillation_hours', 1, 'Oscillation hours\nbeyond the reference'))
    for axis, column, scale, label in panels:
        for treatment in TREATMENTS:
            data = summary[summary['treatment'] == treatment].sort_values('dt')
            if data.empty:
                continue
            axis.plot(plotting.step_positions(data['dt'], steps), scale * data[column], marker='o', markersize=3,
                      color=plotting.TREATMENT_COLORS[treatment], label=plotting.TREATMENT_LABELS[treatment])
        if column in ('excess_negative_hours', 'excess_oscillation_hours'):
            axis.set_yscale('symlog', linthresh=10)
        axis.set_ylabel(label)
    for axis in axes[1]:
        plotting.step_ticks(axis, steps)
        axis.set_xlabel('Routing time step')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2)
    plotting.save(figure, 'era5_annual_errors')
    return


def plot_profiles(profiles: pd.DataFrame, stems: tuple[str, ...], cells: tuple[tuple[str, int], ...]) -> None:
    """Event peak error along main stems, distance above the mouth on x."""
    figure, axes = plt.subplots(1, len(stems), figsize=(plotting.WIDTH, 4.0), sharey=True)
    for axis, stem in zip(np.atleast_1d(axes), stems, strict=True):
        for treatment, dt in cells:
            data = profiles[(profiles['stem'] == stem) & (profiles['treatment'] == treatment)
                            & (profiles['dt'] == dt)].sort_values('km_above_mouth')
            axis.plot(data['km_above_mouth'], 100 * data['peak_error'], color=plotting.TREATMENT_COLORS[treatment],
                      label=f'{plotting.TREATMENT_LABELS[treatment]}, {dt} s', marker='o', markersize=2)
        axis.axhline(0, color=plotting.MUTED, linewidth=0.6)
        axis.invert_xaxis()
        axis.set_title(stem)
    np.atleast_1d(axes)[0].set_ylabel('Event peak error (%)')
    np.atleast_1d(axes)[len(stems) // 2].set_xlabel('km above the mouth')
    handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2)
    plotting.save(figure, 'era5_stem_profiles')
    return


def analyze() -> None:
    """Measure every routed cell against the reference and write the tables the figures are drawn from."""
    network = rr.Network(config.NETWORK_FILE)
    reference_cell = metrics.load_cell(SCENARIO, *config.REFERENCE)
    reference_stats = dict(np.load(reference_cell.directory / 'annual_stats.npz'))
    cells = finished_cells()
    annual = pd.concat([
        annual_errors(metrics.load_cell(SCENARIO, t, dt), reference_stats, reference_cell.river_ids,
                      river_classes(network, dt)) for t, dt in cells
    ], ignore_index=True)
    reference_row = {'treatment': 'reference', 'dt': config.REFERENCE[1],
                     'negative_hours': int(reference_stats['negative_hours'].sum()),
                     'negative_volume_m3': float(reference_stats['negative_volume'].sum()),
                     'oscillation_hours': int(reference_stats['oscillation_hours'].sum()),
                     'rivers_negative': int((reference_stats['negative_hours'].sum(axis=0) > 0).sum())}
    metrics.save_table(pd.DataFrame([reference_row]), 'era5_reference_artifacts')
    metrics.save_table(cost_table(cells), 'era5_cost')
    metrics.save_table(negative_table(cells, network), 'era5_negative')
    annual_summary = summarize_annual(annual, ['treatment', 'dt'])
    metrics.save_table(annual_summary, 'era5_annual_summary')
    metrics.save_table(summarize_annual(annual, ['treatment', 'dt', 'class']), 'era5_annual_by_class')
    events = json.loads((config.INPUTS / 'events.json').read_text())
    event_frames = [event_errors(metrics.load_cell(SCENARIO, t, dt), reference_cell, river_classes(network, dt), e)
                    for t, dt in cells for e in events]
    event_rows = pd.concat(event_frames, ignore_index=True)
    metrics.save_table(summarize_events(event_rows, ['event', 'treatment', 'dt']), 'era5_events_summary')
    metrics.save_table(summarize_events(event_rows, ['event', 'treatment', 'dt', 'class']), 'era5_events_by_class')
    profiles = pd.concat([stem_profiles(cells, reference_cell, e) for e in events], ignore_index=True)
    metrics.save_table(profiles, 'era5_profiles')
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(annual_summary.round(4).to_string(index=False))
    return


def draw() -> None:
    """Draw the era5 figures from the analysis tables, at the steps of config.DT_ROUTING."""
    plot_annual(metrics.shown(metrics.load_table('era5_annual_summary')))
    profiles = metrics.shown(metrics.load_table('era5_profiles'))
    profiles = profiles[profiles['event'] == PROFILE_EVENT]
    measured = set(zip(profiles['treatment'], profiles['dt']))
    plot_profiles(profiles, PROFILE_STEMS, tuple(cell for cell in PROFILE_CELLS if cell in measured))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', nargs='?', default='analyze', choices=('analyze', 'draw'))
    args = parser.parse_args()
    plotting.apply_style()
    if args.action == 'analyze':
        analyze()
    draw()
