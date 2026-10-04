"""
Tables of the effect of every treatment at every routing time step on peak flow, peak timing, and volume, for the
manuscript: one from the ERA5 period of the matrix (annual peaks of every reach) and one from the synthetic burst (the peak of
every reach). Each cell gives the bias and the mean absolute value of a measure against the reference.

Reads tables/synthetic_summary.csv and tables/era5_annual_summary.csv (scripts 08 and 10), replaces the tables between
the treatment_effects markers in paper-claude/manuscript.md, and writes
tables/treatment_effects_<scenario>.csv. Cells not yet run are marked with an ellipsis. At 30 s no reach is
too short, so the short-reach treatments are the standard network and the stabilized network is the substeps network;
those cells are marked a rather than run. The stabilized network at 1 min, which differs from the substeps network
only in the four reaches too short at that step, was cut from the matrix and is marked b.
The table captions explain both marks.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/14_treatment_tables.py
"""

from pathlib import Path

import pandas as pd

from stability import config, plotting

ROWS = ('standard', 'substeps', 'substeps-xadj', 'subcycles', 'stabilized', 'inflate-k', 'merge')
SAME_AS_STANDARD = ('subcycles', 'inflate-k', 'merge')  # not run at a dt with no reach too short
CUT_CELLS = (('stabilized', 60),)  # cut from the ERA5 matrix, scripts/04_run_matrix.py
STEP_LABELS = {30: '30 s', 60: '1 min', 300: '5 min', 600: '10 min', 900: '15 min', 1800: '30 min', 3600: '60 min'}
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


def summary_of(scenario: str) -> pd.DataFrame:
    """The per-cell summary of a scenario, from the analysis tables."""
    if scenario == config.ERA5_SCENARIO:
        table = pd.read_csv(config.TABLES / 'era5_annual_summary.csv')
    else:
        table = pd.read_csv(config.TABLES / 'synthetic_summary.csv')
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
    """One cell: the value, or a mark for a cell that is another treatment's (a), cut (b), or not yet routed (…)."""
    _, column, scale, decimals, signed = measure
    if (treatment, dt) not in summary.index:
        if dt == min(STEP_LABELS) and treatment in (*SAME_AS_STANDARD, 'stabilized'):
            return 'a'
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
            row.update({label: cell_text(summary, treatment, dt, measure) for dt, label in STEP_LABELS.items()})
            rows.append(row)
    table = pd.DataFrame(rows)
    if table.shape != (len(ROWS) * len(MEASURES[scenario]), 2 + len(STEP_LABELS)):
        raise ValueError('the table must have a row per treatment and measure and a column per time step')
    return table


def to_markdown(table: pd.DataFrame) -> str:
    """A pipe table of the cells, as the manuscript writes tables."""
    header = '| ' + ' | '.join(table.columns) + ' |'
    rule = '|' + '|'.join('---' for _ in table.columns) + '|'
    body = ['| ' + ' | '.join(str(value) for value in row) + ' |' for row in table.itertuples(index=False)]
    if not body:
        raise ValueError('an empty table')
    return '\n'.join([header, rule, *body]) + '\n'


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
