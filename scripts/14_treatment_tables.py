"""
Tables of the effect of every treatment at every routing time step on peak flow, peak timing, and volume, for the
manuscript: one from the ERA5 period of the matrix (annual peaks of every reach) and one from the synthetic burst (the peak of
every reach). Each cell gives the bias and the mean absolute value of a measure against the reference.

Reads tables/synthetic_summary.csv and tables/era5_annual_summary.csv (scripts 08 and 10), replaces the tables between
the treatment_effects markers in paper-claude/manuscript.md, and writes tables/treatment_effects_<scenario>.csv, with a
column per step of config.DT_ROUTING, so a step left out needs nothing analyzed or routed again. Cells not yet run are
marked with an ellipsis. The stabilized network at 1 min, which differs from the substeps network only in the four
reaches too short at that step, was cut from the matrix and is marked b, as the table captions explain.

Run with the project environment:  uv run python scripts/14_treatment_tables.py
"""

from pathlib import Path

import pandas as pd

from stability import config, metrics, plotting

ROWS = ('standard', 'substeps', 'substeps-xadj', 'subcycles', 'stabilized', 'inflate-k', 'merge')
CUT_CELLS = (('stabilized', 60),)  # cut from the ERA5 matrix, scripts/04_run_matrix.py
# each measure as (row label, column of the summary, scale, decimals, signed)
MEASURES = {
    'synthetic-burst': (
        ('Peak, bias (% of rise)', 'peak_error_median', 100, 1, True),
        ('Peak, mean abs. error (% of rise)', 'peak_error_abs_mean', 100, 1, False),
        ('Peak time, bias (h)', 'timing_error_mean_h', 1, 2, True),
        ('Peak time, mean abs. error (h)', 'timing_error_abs_mean_h', 1, 2, False),
        ('Volume, mean abs. error (%)', 'volume_error_abs_mean', 100, 3, False),
    ),
    config.ERA5_SCENARIO: (
        ('Annual peak, bias (%)', 'peak_error_median', 100, 2, True),
        ('Annual peak, mean abs. error (%)', 'peak_error_abs_mean', 100, 2, False),
        ('Peak time, bias (h)', 'timing_mean_h', 1, 2, True),
        ('Peak time, mean abs. error (h)', 'timing_abs_mean_h', 1, 2, False),
        ('Annual volume, mean abs. error (%)', 'volume_error_abs_mean', 100, 3, False),
    ),
}
MINUS = '\u2212'


def step_labels() -> dict[int, str]:
    """The column of each routing time step of config.DT_ROUTING, in minutes."""
    if any(dt % 60 for dt in config.DT_ROUTING):
        raise ValueError('the tables give every routing time step in whole minutes')
    return {dt: f'{dt // 60} min' for dt in config.DT_ROUTING}


def summary_of(scenario: str) -> pd.DataFrame:
    """The per-cell summary of a scenario, from the analysis tables, at the steps of config.DT_ROUTING."""
    if scenario == config.ERA5_SCENARIO:
        table = metrics.shown(metrics.load_table('era5_annual_summary'))
    else:
        table = metrics.shown(metrics.load_table('synthetic_summary'))
        table = table[table['scenario'] == scenario]
    if table.empty:
        raise ValueError(f'no summary for {scenario}')
    return table.set_index(['treatment', 'dt'])


def number_text(value: float, decimals: int, signed: bool) -> str:
    """A number with a true minus sign, a plus sign when signed, and no sign on a zero that rounds away."""
    text = f'{abs(value):.{decimals}f}'
    if float(text) == 0:
        return text
    if value < 0:
        return MINUS + text
    return ('+' if signed else '') + text


def cell_text(summary: pd.DataFrame, treatment: str, dt: int, measure: tuple) -> str:
    """One cell: the value, or a mark for a cell cut from the matrix (b) or not yet routed (…)."""
    _, column, scale, decimals, signed = measure
    if (treatment, dt) not in summary.index:
        if (treatment, dt) in CUT_CELLS:
            return 'b'
        return '…'
    return number_text(scale * summary.loc[(treatment, dt)][column], decimals, signed)


def build(scenario: str) -> pd.DataFrame:
    """The table of a scenario: a row per treatment and measure, a column per routing time step."""
    summary = summary_of(scenario)
    rows = []
    for treatment in ROWS:
        for index, measure in enumerate(MEASURES[scenario]):
            row = {'Treatment': plotting.TREATMENT_LABELS[treatment] if index == 0 else '', 'Measure': measure[0]}
            row.update({label: cell_text(summary, treatment, dt, measure) for dt, label in step_labels().items()})
            rows.append(row)
    table = pd.DataFrame(rows)
    if table.shape != (len(ROWS) * len(MEASURES[scenario]), 2 + len(config.DT_ROUTING)):
        raise ValueError('the table must have a row per treatment and measure and a column per time step')
    return table


def to_markdown(table: pd.DataFrame) -> str:
    """A pipe table of the cells, its columns padded to one width as the manuscript writes tables."""
    if table.empty:
        raise ValueError('an empty table')
    widths = [max(len(str(value)) for value in (column, *table[column])) for column in table.columns]

    def line(values) -> str:
        return '| ' + ' | '.join(str(value).ljust(width) for value, width in zip(values, widths, strict=True)) + ' |'

    rule = '|' + '|'.join('-' * (width + 2) for width in widths) + '|'
    return '\n'.join([line(table.columns), rule, *(line(row) for row in table.itertuples(index=False))]) + '\n'


def refresh_manuscript(manuscript: Path, name: str, markdown: str) -> bool:
    """Replace the table between the begin and end markers of ``name`` in the manuscript; False if it has none."""
    text = manuscript.read_text()
    begin, end = f'<!-- begin {name} -->', f'<!-- end {name} -->'
    if text.count(begin) != 1 or text.count(end) != 1:
        return False
    head, rest = text.split(begin)
    _, tail = rest.split(end)
    manuscript.write_text(f'{head}{begin}\n\n{markdown}\n{end}{tail}')
    return True


if __name__ == '__main__':
    for name in MEASURES:
        effects = build(name)
        effects.to_csv(config.TABLES / f'treatment_effects_{name}.csv', index=False)
        refreshed = refresh_manuscript(config.MANUSCRIPT, f'treatment_effects_{name}', to_markdown(effects))
        print(f'{name}: written to tables/' + (', refreshed in the manuscript' if refreshed else ''))
