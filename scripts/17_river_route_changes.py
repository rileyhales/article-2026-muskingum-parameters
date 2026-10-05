"""
Compare river-route before and after the changes this study made to it, each version against its own time-converged
reference, so every difference is attributable to one change:

- reporting: the standard network, whose discharge before is the plain mean of the discharges computed at the
  routing steps of each hour and after is their trapezoidal mean (Eq. 9)
- river-route: its own stabilized network, which before divides a river too long for dt into ceil(2kx/dt) substeps of
  its own x (the stabilized treatment of the matrix) and after into the fewest substeps with x adjusted by Eq. 16
- substep count: after only, substeps with x adjusted, ceil(2kx/dt) of them (substeps-xadj) against the fewest
  (substeps-fewest), without subcycles

Cells are routed for the synthetic burst, the square wave, and ERA5 2002-2011 at every step of the matrix, and write to
``data/river-route-changes/<code>/<scenario>/<treatment>__dt<NNNN>s`` as the matrix cells do (scripts/04_run_matrix.py).
The before code is the river-route that routed the matrix, copied into vendor/river-route-matrix; a run checks that the
river-route it imported is the version it was asked for.

``analyze`` writes tables/river_route_changes_<scenario>.csv with every measure of every cell, then draws. ``draw``
writes tables/river_route_changes.csv with the measures of the manuscript, which it also writes between the
river_route_changes markers of the manuscript, and the figure river_route_changes, from those tables alone, at the steps
of config.DT_ROUTING, so a step left out needs nothing analyzed or routed again.

Run with the river-route environment (python is ../river-route/.venv/bin/python):
    PYTHONPATH=vendor/river-route-matrix python scripts/17_river_route_changes.py run --code before --scenario era5-2002-2011
    python scripts/17_river_route_changes.py run --code after --scenario synthetic-burst
    python scripts/17_river_route_changes.py analyze
    python scripts/17_river_route_changes.py draw
"""

import argparse
import importlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import river_route as rr

from stability import config, engine, forcing, metrics, plotting

run_matrix = importlib.import_module('04_run_matrix')
synthetic = importlib.import_module('08_analyze_synthetic')
era5 = importlib.import_module('10_analyze_era5')

SCENARIOS = (*forcing.SYNTHETIC_PERTURBATIONS, config.ERA5_SCENARIO)
TREATMENTS = {
    'before': ('standard', 'river-route'),
    'after': ('standard', 'river-route', 'substeps-xadj', 'substeps-fewest'),
}
# the rows of the tables and figure: each (version, treatment) and its label
ROWS = (
    ('before', 'standard', 'Standard, before'),
    ('after', 'standard', 'Standard, after'),
    ('before', 'river-route', 'river-route, before'),
    ('after', 'river-route', 'river-route, after'),
    ('after', 'substeps-xadj', 'Substeps, x adjusted, ⌈2kx/Δt⌉'),
    ('after', 'substeps-fewest', 'Substeps, x adjusted, fewest'),
)
STANDARD_WORK_PER_HOUR = 29_405  # reach-steps per simulated hour of the standard network at 1 h
MINUS = '−'


def cells(code: str) -> list[tuple[str, int]]:
    """Every (treatment, dt) a version routes, the reference first."""
    if code not in TREATMENTS:
        raise ValueError(f'unknown code {code}')
    found = [config.REFERENCE] + [(t, dt) for t in TREATMENTS[code] for dt in config.DT_ROUTING]
    if len(found) != 1 + len(TREATMENTS[code]) * len(config.DT_ROUTING):
        raise ValueError('one cell per treatment and step')
    return found


def run(code: str, scenario: str, threads: int, force: bool) -> None:
    """Route every cell of a version for a scenario, with the river-route of that version."""
    imported = engine.river_route_code()
    if imported != code:
        raise SystemExit(
            f'river_route is the {imported} version, imported from {rr.__file__}; the before version needs '
            f'PYTHONPATH={config.CODES["before"].relative_to(config.ROOT)}'
        )
    for treatment, dt in cells(code):
        run_matrix.run_cell(scenario, treatment, dt, threads, force, root=config.CHANGES / code)
    return


def work_per_hour(code: str, scenario: str) -> pd.DataFrame:
    """Reach-steps per simulated hour of every cell, relative to the standard network at 1 h."""
    rows = []
    for treatment, dt in cells(code)[1:]:
        meta = metrics.load_cell(scenario, treatment, dt, config.CHANGES / code).meta
        per_hour = meta['reach_steps_per_routing_step'] * config.DT_RUNOFF / dt
        rows.append({'code': code, 'treatment': treatment, 'dt': dt, 'work': per_hour / STANDARD_WORK_PER_HOUR,
                     'reaches': meta['reaches'], 'subcycled_rivers': meta['subcycled_rivers']})
    if not rows:
        raise ValueError('no cells')
    return pd.DataFrame(rows)


def summarize_synthetic(code: str, scenario: str, network: rr.Network, area: np.ndarray) -> pd.DataFrame:
    """The synthetic summary of scripts/08 for every cell of a version: peaks, timing, volume, artifacts."""
    rivers = synthetic.compare(scenario, network, area, root=config.CHANGES / code, shown=TREATMENTS[code])
    summary = synthetic.summarize(rivers[rivers['treatment'] != 'reference'], ['treatment', 'dt'])
    short = synthetic.summarize(rivers[rivers['class'] == 'too short'], ['treatment', 'dt'])
    short = short[['treatment', 'dt', 'peak_error_abs_mean', 'rivers_with_negative_hours']].rename(
        columns={'peak_error_abs_mean': 'short_peak_error_abs_mean', 'rivers_with_negative_hours': 'short_negative'})
    summary = summary.merge(short, on=['treatment', 'dt'], how='left')
    return summary.rename(columns={'timing_error_mean_h': 'timing_mean_h', 'timing_error_abs_mean_h': 'timing_abs_mean_h'})


def summarize_era5(code: str, network: rr.Network) -> pd.DataFrame:
    """The annual summary of scripts/10 for every cell of a version: annual peaks, timing, volume, negatives."""
    root = config.CHANGES / code
    reference = metrics.load_cell(config.ERA5_SCENARIO, *config.REFERENCE, root)
    stats = dict(np.load(reference.directory / 'annual_stats.npz'))
    annual = pd.concat([
        era5.annual_errors(metrics.load_cell(config.ERA5_SCENARIO, t, dt, root), stats, reference.river_ids,
                           era5.river_classes(network, dt)) for t, dt in cells(code)[1:]
    ], ignore_index=True)
    summary = era5.summarize_annual(annual, ['treatment', 'dt'])
    short = era5.summarize_annual(annual[annual['class'] == 'too short'], ['treatment', 'dt'])
    short = short[['treatment', 'dt', 'peak_error_abs_mean', 'excess_negative_hours']].rename(
        columns={'peak_error_abs_mean': 'short_peak_error_abs_mean', 'excess_negative_hours': 'short_negative'})
    return summary.merge(short, on=['treatment', 'dt'], how='left')


def summarize(scenario: str, network: rr.Network, area: np.ndarray) -> pd.DataFrame:
    """Both versions of a scenario, with the work of every cell."""
    frames = []
    for code in TREATMENTS:
        if scenario == config.ERA5_SCENARIO:
            summary = summarize_era5(code, network)
        else:
            summary = summarize_synthetic(code, scenario, network, area)
        frames.append(summary.assign(code=code).merge(work_per_hour(code, scenario), on=['code', 'treatment', 'dt']))
    table = pd.concat(frames, ignore_index=True)
    if table.empty:
        raise ValueError(f'nothing to summarize for {scenario}')
    return table


MEASURES = (  # (scenario, row label, column of its summary, scale, decimals or None for work, signed)
    ('synthetic-burst', 'Burst peak, bias (% of rise)', 'peak_error_median', 100, 1, True),
    ('synthetic-burst', 'Burst peak time, bias (h)', 'timing_mean_h', 1, 2, True),
    (config.ERA5_SCENARIO, 'Annual peak, bias (%)', 'peak_error_median', 100, 2, True),
    (config.ERA5_SCENARIO, 'Annual peak, mean abs. error (%)', 'peak_error_abs_mean', 100, 2, False),
    (config.ERA5_SCENARIO, 'Annual peak time, bias (h)', 'timing_mean_h', 1, 2, True),
    ('synthetic-burst', 'Work, vs standard at 1 h', 'work', 1, None, False),
)


def work_text(value: float) -> str:
    """Work to three significant figures at most, with thousands separated."""
    if value <= 0:
        raise ValueError('work is positive')
    decimals = 0 if value >= 100 else 1 if value >= 10 else 2
    return f'{value:,.{decimals}f}'


def build_table(summaries: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """A row per version, treatment, and measure, and a column per routing time step."""
    tables = importlib.import_module('14_treatment_tables')
    indexed = {scenario: summary.set_index(['code', 'treatment', 'dt']) for scenario, summary in summaries.items()}
    rows = []
    for code, treatment, label in ROWS:
        for index, (scenario, name, column, scale, decimals, signed) in enumerate(MEASURES):
            row = {'Version and treatment': label if index == 0 else '', 'Measure': name}
            for dt, step in tables.step_labels().items():
                if (code, treatment, dt) not in indexed[scenario].index:
                    row[step] = '…'
                    continue
                value = indexed[scenario].loc[(code, treatment, dt), column]
                row[step] = work_text(value) if decimals is None else tables.number_text(scale * value, decimals, signed)
            rows.append(row)
    table = pd.DataFrame(rows)
    if len(table) != len(ROWS) * len(MEASURES):
        raise ValueError('one row per version, treatment, and measure')
    return table


def plot(summaries: dict[str, pd.DataFrame]) -> None:
    """Peak time bias and median peak error against dt for the burst and ERA5, and the work of every cell."""
    figure, axes = plt.subplots(3, 2, figsize=(plotting.WIDTH, 9.0))
    steps = np.asarray(config.DT_ROUTING)
    panels = (('peak time', 'timing_mean_h', 1, 'Mean peak time\nbias (h)'),
              ('peak', 'peak_error_median', 100, 'Median peak\nbias (%)'))
    for column, scenario in enumerate(('synthetic-burst', config.ERA5_SCENARIO)):
        summary = summaries[scenario]
        for row, (_, measure, scale, label) in enumerate(panels):
            axis = axes[row, column]
            for code, treatment, name in ROWS:
                if treatment == 'substeps-fewest' or (row == 0 and treatment == 'substeps-xadj'):
                    continue
                data = summary[(summary['code'] == code) & (summary['treatment'] == treatment)].sort_values('dt')
                axis.plot(plotting.step_positions(data['dt'], steps), scale * data[measure], label=name,
                          color=plotting.TREATMENT_COLORS[treatment], marker=plotting.TREATMENT_MARKERS[treatment],
                          markersize=5, linestyle='--' if code == 'before' else '-',
                          markerfacecolor='none' if code == 'before' else plotting.TREATMENT_COLORS[treatment])
            axis.axhline(0, color=plotting.MUTED, linewidth=0.6)
            plotting.step_ticks(axis, steps)
            axis.set_ylabel(label if column == 0 else '')
        axes[0, column].set_title('Synthetic burst' if column == 0 else 'ERA5 annual peaks')
    work = summaries['synthetic-burst']
    for code, treatment, name in ROWS[1:]:
        if code == 'before' and treatment == 'standard':
            continue
        data = work[(work['code'] == code) & (work['treatment'] == treatment)].sort_values('dt')
        axes[2, 0].plot(plotting.step_positions(data['dt'], steps), data['work'], label=name,
                        color=plotting.TREATMENT_COLORS[treatment], marker=plotting.TREATMENT_MARKERS[treatment],
                        markersize=5, linestyle='--' if code == 'before' else '-',
                        markerfacecolor='none' if code == 'before' else plotting.TREATMENT_COLORS[treatment])
    axes[2, 0].set_yscale('log')
    axes[2, 0].set_ylabel('Work per hour, vs\nstandard at 1 h')
    plotting.step_ticks(axes[2, 0], steps)
    axes[2, 1].axis('off')
    handles, labels = axes[2, 0].get_legend_handles_labels()
    standard = axes[0, 0].get_legend_handles_labels()
    axes[2, 1].legend(standard[0][:1] + handles, standard[1][:1] + labels, loc='center')
    plotting.save(figure, 'river_route_changes')
    return


def analyze() -> None:
    """Summarize every cell of both versions and write the tables the manuscript table and figure are drawn from."""
    network = rr.Network(config.NETWORK_FILE)
    area = pd.read_parquet(config.METADATA_FILE, columns=['DSContArea'])['DSContArea'].to_numpy() / 1e6
    shown = ['code', 'treatment', 'dt', 'peak_error_median', 'peak_error_abs_mean', 'timing_mean_h',
             'timing_abs_mean_h', 'volume_error_abs_mean', 'short_peak_error_abs_mean', 'short_negative', 'work']
    for scenario in SCENARIOS:
        summary = summarize(scenario, network, area)
        metrics.save_table(summary, f'river_route_changes_{scenario}')
        with pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.max_rows', 200):
            print(scenario)
            print(summary[shown].round(4).to_string(index=False))
    return


def draw() -> None:
    """The manuscript table and the figure of both versions, from the analysis tables, at config.DT_ROUTING."""
    plotting.apply_style()
    tables = importlib.import_module('14_treatment_tables')
    summaries = {s: metrics.shown(metrics.load_table(f'river_route_changes_{s}')) for s in SCENARIOS}
    table = build_table(summaries)
    metrics.save_table(table, 'river_route_changes')
    refreshed = tables.refresh_manuscript(config.MANUSCRIPT, 'river_route_changes', tables.to_markdown(table))
    print('tables/river_route_changes.csv written' + (', manuscript refreshed' if refreshed else ''))
    plot(summaries)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', choices=('run', 'analyze', 'draw'))
    parser.add_argument('--code', choices=tuple(TREATMENTS))
    parser.add_argument('--scenario', choices=SCENARIOS)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--force', action='store_true', help='rerun a cell whose outputs exist')
    args = parser.parse_args()
    if args.action == 'run':
        if not args.code or not args.scenario:
            parser.error('run needs --code and --scenario')
        run(args.code, args.scenario, args.threads, args.force)
    else:
        if args.action == 'analyze':
            analyze()
        draw()
