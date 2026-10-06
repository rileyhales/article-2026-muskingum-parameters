"""
Measure what every cell of the matrix costs to route, without running its whole period. The engine times its stages
as it goes (stability.engine.Timer), so every run of the matrix records its routing seconds per simulated year in its
meta.json; those runs share the machine, so this script measures the same quantity alone: one month of ERA5 catchment
runoff is read into memory once, and each treated network is routed over it with a recorder that discards everything,
at one thread and at eight, keeping the fastest of REPEATS runs. The cost is reported as routing seconds per simulated
year and as reach-steps per simulated hour, the work the treatment implies.

Run after the matrix, with nothing else running:
    uv run python scripts/12_benchmark_cost.py [analyze|draw]

``analyze`` times every cell and writes tables/cost_benchmark.csv, then draws. ``draw`` draws figures cost_benchmark
and cost_tradeoff (which reads the synthetic and ERA5 summary tables of scripts 08 and 10) from the tables alone, at
the steps of config.DT_ROUTING, so a step left out needs nothing timed or routed again.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from stability import config, engine, forcing, metrics, plotting, treatments

MONTH = '2011-06'  # routed as two copies, so the steady rate excludes the first, which loads the kernels
REPEATS = 3
THREADS = (1, 8)
HOURS_PER_YEAR = 8760


class Discard:
    """A recorder that keeps nothing, so only routing is timed."""

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None:
        if discharge.shape[1] != dates.shape[0]:
            raise ValueError('the discharge must have one column per date')


def cells() -> list[tuple[str, int]]:
    """Every (treatment, dt) of the ERA5 matrix that has been run, the reference first."""
    found = [config.REFERENCE]
    for treatment in config.TREATMENTS:
        for dt in config.DT_ROUTING:
            if treatment != 'reference' and (config.run_dir(config.ERA5_SCENARIO, treatment, dt) / 'meta.json').exists():
                found.append((treatment, dt))
    if len(found) < 2:
        raise FileNotFoundError('run the matrix before benchmarking it')
    return found


def time_cell(treatment: str, dt: int, table: pd.DataFrame, chunks: list, threads: int) -> dict:
    """
    The fastest of REPEATS routings of one month over one cell's network, timed by the engine's own timer, which counts
    only river-route's kernels, and the work the treatment implies.
    """
    treated = treatments.treat(treatment, dt, table)
    engine.route(treated, dt, config.DT_RUNOFF, chunks, Discard(), threads=threads)  # compile and warm caches
    runs = [engine.route(treated, dt, config.DT_RUNOFF, chunks, Discard(), threads=threads) for _ in range(REPEATS)]
    fastest = min(run['steady_routing_seconds_per_simulated_year'] for run in runs)
    if not np.isfinite(fastest) or fastest <= 0:
        raise ValueError('a cell must be timed')
    return {
        'treatment': treatment, 'dt': dt, 'threads': threads, 'reaches': runs[0]['reaches'],
        'subcycled_rivers': runs[0]['subcycled_rivers'],
        'reach_steps_per_hour': runs[0]['reach_steps_per_routing_step'] * config.DT_RUNOFF // dt,
        'seconds_per_simulated_year': fastest,
    }


def plot(costs: pd.DataFrame) -> None:
    """Seconds to route a simulated year against dt, per treatment, at eight threads."""
    figure, axis = plt.subplots(figsize=(plotting.WIDTH, 4.0))
    steps = np.asarray(config.DT_ROUTING)
    data = costs[costs['threads'] == max(THREADS)]
    for treatment in config.TREATMENTS:
        rows = data[data['treatment'] == treatment].sort_values('dt')
        if rows.empty or treatment == 'reference':
            continue
        axis.plot(plotting.step_positions(rows['dt'], steps), rows['seconds_per_simulated_year'], marker='o',
                  markersize=3, color=plotting.TREATMENT_COLORS[treatment], label=plotting.TREATMENT_LABELS[treatment])
    reference = data[data['treatment'] == 'reference']['seconds_per_simulated_year'].iloc[0]
    axis.axhline(reference, color=plotting.INK, linestyle='--', linewidth=1, label='Reference')
    axis.set_yscale('log')
    plotting.step_ticks(axis, steps)
    axis.set_xlabel('Routing time step')
    axis.set_ylabel(f'Seconds per simulated year\n({max(THREADS)} threads, Columbia)')
    axis.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
    plotting.save(figure, 'cost_benchmark')
    return


def plot_tradeoff(costs: pd.DataFrame) -> None:
    """Accuracy against cost: mean absolute peak error of every cell against the seconds to route a simulated year."""
    errors = (('synthetic_summary', 'synthetic-burst', 'Synthetic burst: mean |peak error| (% of rise)'),
              ('era5_annual_summary', None, 'ERA5: mean |annual peak error| (%)'))
    figure, axes = plt.subplots(1, 2, figsize=(plotting.WIDTH, 4.0))
    timing = costs[costs['threads'] == max(THREADS)][['treatment', 'dt', 'seconds_per_simulated_year']]
    for axis, (table, scenario, label) in zip(axes, errors, strict=True):
        summary = metrics.shown(metrics.load_table(table))
        if scenario is not None:
            summary = summary[summary['scenario'] == scenario]
        joined = summary.merge(timing, on=['treatment', 'dt'])
        for treatment in config.TREATMENTS:
            rows = joined[joined['treatment'] == treatment].sort_values('dt')
            if rows.empty or treatment == 'reference':
                continue
            axis.plot(rows['seconds_per_simulated_year'], 100 * rows['peak_error_abs_mean'], marker='o',
                      markersize=3, color=plotting.TREATMENT_COLORS[treatment],
                      label=plotting.TREATMENT_LABELS[treatment])
        axis.set_xscale('log')
        axis.set_yscale('symlog', linthresh=0.01)
        axis.set_xlabel(f'Seconds per simulated year ({max(THREADS)} threads)')
        axis.set_ylabel(label)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2)
    plotting.save(figure, 'cost_tradeoff')
    return


def analyze() -> None:
    """Time every cell of the matrix and write the table the figures are drawn from."""
    network_table = pd.read_parquet(config.NETWORK_FILE)
    month = next(forcing.era5_chunks(forcing.era5_files(MONTH, MONTH), network_table['riverId'].to_numpy()))
    month = [month, month]
    rows = []
    for cell in cells():
        for count in THREADS:
            rows.append(time_cell(*cell, network_table, month, count))
            print(rows[-1], flush=True)
    benchmark = pd.DataFrame(rows)
    standard_hour = benchmark[(benchmark['treatment'] == 'standard') & (benchmark['dt'] == 3600)]
    for count in THREADS:
        base = standard_hour[standard_hour['threads'] == count]['seconds_per_simulated_year'].iloc[0]
        chosen = benchmark['threads'] == count
        benchmark.loc[chosen, 'time_vs_standard_1h'] = benchmark.loc[chosen, 'seconds_per_simulated_year'] / base
    metrics.save_table(benchmark, 'cost_benchmark')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', nargs='?', default='analyze', choices=('analyze', 'draw'))
    args = parser.parse_args()
    plotting.apply_style()
    if args.action == 'analyze':
        analyze()
    costs = metrics.shown(metrics.load_table('cost_benchmark'))
    plot(costs)
    plot_tradeoff(costs)
