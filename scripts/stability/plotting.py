"""
The figure style of the manuscript: one color per treatment, held across every figure, and a four-step ordinal ramp
for the routing time steps a hydrograph figure compares. Both palettes passed the categorical and ordinal checks of
the data visualization validator on the light surface. Three treatment hues sit below 3:1 contrast with the surface,
so every figure that uses them carries a legend.
"""

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.layout_engine import TightLayoutEngine

from . import config

__all__ = [
    'TREATMENT_COLORS', 'TREATMENT_LABELS', 'TREATMENT_MARKERS', 'HYDROFABRIC_COLORS', 'HYDROFABRIC_MARKERS', 'DT_RAMP', 'INK', 'MUTED',
    'GRID', 'SURFACE', 'WIDTH', 'FONT_SIZE', 'LEGEND_FONT_SIZE', 'apply_style', 'save', 'step_label',
    'step_positions', 'step_ticks',
]

WIDTH = 6.0  # inches: every figure is as wide as the text block of a letter page, and 3 to 4 inches tall per row
FONT_SIZE = 12  # points, for every text of every figure but legends
LEGEND_FONT_SIZE = 10
SURFACE = '#fcfcfb'
INK = '#0b0b0b'
MUTED = '#52514e'
GRID = '#e4e3df'
TREATMENT_COLORS = {
    'reference': INK,
    'standard': '#2a78d6',
    'substeps': '#eb6834',
    'substeps-xadj': '#1baf7a',
    'subcycles': '#eda100',
    'stabilized': '#e87ba4',
    'inflate-k': '#008300',
    'merge': '#4a3aa7',
    # slot 8, for river-route's own stabilized network in the comparison of its versions; it fails the
    # normal-vision floor against stabilized, so the two are never drawn together
    'river-route': '#e34948',
}
# a marker per treatment for figures whose hues fall in the colorblind floor band, as red and aqua do
TREATMENT_MARKERS = {'standard': 'o', 'river-route': 's', 'substeps-xadj': '^', 'substeps-fewest': 'v'}
TREATMENT_LABELS = {
    'reference': 'Reference',
    'standard': 'Standard',
    'substeps': 'Substeps',
    'substeps-xadj': 'Substeps, x adjusted',
    'subcycles': 'Subcycles',
    'stabilized': 'Stabilized',
    'inflate-k': 'Inflate k',
    'merge': 'Merge',
    'substeps-fewest': 'Substeps, x adjusted, fewest',
    'river-route': 'river-route',
}
# one color per hydrofabric of the census figures, on the categorical slots they have always held
HYDROFABRIC_COLORS = {
    'TDX-Hydro': '#2a78d6',
    'HydroRIVERS': '#1baf7a',
    'MERIT-Basins': '#eda100',
    # slots 6 and 7; slot 2, orange, fails the normal-vision floor against MERIT-Basins
    'NHDPlus V2': '#008300',
    'NHDPlus HR': '#4a3aa7',
}
# and one marker, so a hydrofabric is never told by color alone
HYDROFABRIC_MARKERS = {'TDX-Hydro': 'o', 'HydroRIVERS': '^', 'MERIT-Basins': 'D', 'NHDPlus V2': 'v', 'NHDPlus HR': 's'}
DT_RAMP = ('#86b6ef', '#3987e5', '#1c5cab', '#0d366b')  # light to dark, for up to four time steps


def apply_style() -> None:
    """Set the matplotlib defaults every figure script uses."""
    if not set(config.TREATMENTS) <= set(TREATMENT_COLORS):
        raise ValueError('every treatment of the matrix needs a color')
    if len(set(TREATMENT_COLORS.values())) != len(TREATMENT_COLORS):
        raise ValueError('treatment colors must be distinct')
    mpl.rcParams.update({
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
        'font.family': 'sans-serif', 'font.size': FONT_SIZE, 'axes.titlesize': FONT_SIZE, 'axes.labelsize': FONT_SIZE,
        'xtick.labelsize': FONT_SIZE, 'ytick.labelsize': FONT_SIZE,
        'legend.fontsize': LEGEND_FONT_SIZE, 'legend.title_fontsize': LEGEND_FONT_SIZE,
        'figure.titlesize': FONT_SIZE, 'figure.figsize': (WIDTH, 4.0), 'figure.autolayout': True,
        'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'axes.linewidth': 0.6,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'grid.linestyle': '-',
        'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.labelcolor': INK, 'ytick.labelcolor': INK,
        'lines.linewidth': 1.5, 'lines.solid_capstyle': 'round', 'lines.solid_joinstyle': 'round',
        'legend.frameon': False, 'figure.dpi': 150, 'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })
    return


def save(figure: plt.Figure, name: str) -> None:
    """Write a figure as PNG and PDF to the figures directory and close it."""
    if not name or '/' in name:
        raise ValueError('a figure is saved by a bare name')
    if not isinstance(figure.get_layout_engine(), TightLayoutEngine):
        raise ValueError(f'{name} must use the tight layout of the style')
    config.FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(config.FIGURES / f'{name}.png')
    figure.savefig(config.FIGURES / f'{name}.pdf')
    plt.close(figure)
    if not (config.FIGURES / f'{name}.png').exists():
        raise OSError(f'{name}.png was not written')
    return


def step_label(seconds: int) -> str:
    """A routing time step as an axis labels it: 30 s, 5 min, 1 h."""
    if seconds <= 0:
        raise ValueError('steps are positive')
    if seconds < 60:
        return f'{seconds} s'
    return f'{seconds // 60} min' if seconds < 3600 else f'{seconds // 3600} h'


def step_positions(values, steps) -> np.ndarray:
    """The position of each step on an axis that spaces the steps evenly, in order."""
    values, steps = np.asarray(values), np.asarray(steps)
    positions = np.searchsorted(steps, values)
    found = steps[np.minimum(positions, steps.shape[0] - 1)]
    if np.any(positions >= steps.shape[0]) or not np.array_equal(found, values):
        raise ValueError('every value must be one of the steps')
    return positions


def step_ticks(axis: plt.Axes, steps) -> None:
    """Label an axis of evenly spaced steps with the steps themselves."""
    if np.any(np.diff(np.asarray(steps)) <= 0):
        raise ValueError('steps must increase')
    axis.set_xticks(np.arange(len(steps)), [step_label(int(s)) for s in steps], rotation=90)
    return
