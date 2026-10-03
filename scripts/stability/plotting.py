"""
The figure style of the manuscript: one color per treatment, held across every figure, and a four-step ordinal ramp
for the routing time steps a hydrograph figure compares. Both palettes passed the categorical and ordinal checks of
the data visualization validator on the light surface. Three treatment hues sit below 3:1 contrast with the surface,
so every figure that uses them carries a legend.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

from . import config

__all__ = ['TREATMENT_COLORS', 'TREATMENT_LABELS', 'DT_RAMP', 'INK', 'MUTED', 'GRID', 'SURFACE', 'apply_style', 'save']

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
}
TREATMENT_LABELS = {
    'reference': 'Reference',
    'standard': 'Standard',
    'substeps': 'Substeps',
    'substeps-xadj': 'Substeps, x adjusted',
    'subcycles': 'Subcycles',
    'stabilized': 'Stabilized',
    'inflate-k': 'Inflate k',
    'merge': 'Merge',
}
DT_RAMP = ('#86b6ef', '#3987e5', '#1c5cab', '#0d366b')  # light to dark, for up to four time steps


def apply_style() -> None:
    """Set the matplotlib defaults every figure script uses."""
    if set(TREATMENT_COLORS) != set(config.TREATMENTS):
        raise ValueError('every treatment needs a color')
    if len(set(TREATMENT_COLORS.values())) != len(TREATMENT_COLORS):
        raise ValueError('treatment colors must be distinct')
    mpl.rcParams.update({
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
        'font.family': 'sans-serif', 'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
        'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'axes.linewidth': 0.6,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'grid.linestyle': '-',
        'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.labelcolor': INK, 'ytick.labelcolor': INK,
        'lines.linewidth': 1.5, 'lines.solid_capstyle': 'round', 'lines.solid_joinstyle': 'round',
        'legend.frameon': False, 'legend.fontsize': 8, 'figure.dpi': 150, 'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })
    return


def save(figure: plt.Figure, name: str) -> None:
    """Write a figure as PNG and PDF to the figures directory and close it."""
    if not name or '/' in name:
        raise ValueError('a figure is saved by a bare name')
    config.FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(config.FIGURES / f'{name}.png')
    figure.savefig(config.FIGURES / f'{name}.pdf')
    plt.close(figure)
    if not (config.FIGURES / f'{name}.png').exists():
        raise OSError(f'{name}.png was not written')
    return
