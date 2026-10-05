"""
Paths, constants, and the simulation matrix shared by every script of the experiment.

The matrix is the product of a forcing scenario, a treatment of the rivers whose Muskingum coefficients are negative,
and a routing time step. Every cell of the matrix writes its outputs to
``data/results/<scenario>/<treatment>__dt<NNNN>s``.
"""

from pathlib import Path

__all__ = [
    'ROOT',
    'DATA',
    'INPUTS',
    'RESULTS',
    'FIGURES',
    'TABLES',
    'MANUSCRIPT',
    'RFS_ROOT',
    'ERA5_ROOT',
    'REGION',
    'COLUMBIA_OUTLET',
    'NETWORK_FILE',
    'METADATA_FILE',
    'WEIGHTS_FILE',
    'ERA5_CATCHMENT_DIR',
    'DT_ROUTING',
    'DT_RUNOFF',
    'TREATMENTS',
    'CHANGE_TREATMENTS',
    'REFERENCE',
    'ERA5_YEARS',
    'ERA5_SPINUP',
    'ERA5_SCENARIO',
    'ERA5_FULL_SCENARIO',
    'ERA5_FULL_YEARS',
    'CHANGES',
    'CODES',
    'run_dir',
]

ROOT = Path(__file__).resolve().parents[2]  # the project, from scripts/stability/config.py
DATA = ROOT / 'data'
INPUTS = DATA / 'inputs'  # the prepared network, forcing, observation reaches, and event windows
RESULTS = DATA / 'results'  # the cells of the simulation matrix
FIGURES = ROOT / 'figures'
TABLES = ROOT / 'tables'
MANUSCRIPT = ROOT / 'paper-claude' / 'manuscript.md'

# the River Forecast System v3 hydrofabric and routing parameters, and the ERA5 runoff
RFS_ROOT = Path.home() / 'data' / 'rfsv3'
ERA5_ROOT = Path.home() / 'data' / 'era5'
REGION = 7020014250  # the TDX-Hydro region holding the Columbia basin
COLUMBIA_OUTLET = 720207783  # riverId of the Columbia at the Pacific

# inputs prepared by scripts/01_prepare_columbia.py and scripts/02_aggregate_era5.py
NETWORK_FILE = INPUTS / 'columbia_network.parquet'
METADATA_FILE = INPUTS / 'columbia_metadata.parquet'
WEIGHTS_FILE = INPUTS / 'columbia_gridweights_era5.nc'
ERA5_CATCHMENT_DIR = INPUTS / 'era5_catchment_runoff'

# the base routing time steps of the matrix, seconds. Every one divides the hourly runoff step. The figures and
# tables show only these steps, and the reference, of what the analysis tables hold.
DT_ROUTING = (60, 300, 600, 900, 1800, 3600)
DT_RUNOFF = 3600
# the ERA5 period of the matrix: ten years holding all four event windows, after three months of spin-up that are routed
# from a dry network and then discarded. The first cells were routed over 2000-2019, from a spin-up beginning 1999-10;
# scripts/15_restrict_era5_period.py cuts those to this period.
ERA5_YEARS = tuple(range(2002, 2012))
ERA5_SPINUP = ('2001-10', '2001-12')
ERA5_SCENARIO = 'era5-2002-2011'
ERA5_FULL_SCENARIO = 'era5-2000-2019'
ERA5_FULL_YEARS = tuple(range(2000, 2020))

# the treatments of rivers whose coefficients are negative at a time step
TREATMENTS = {
    'standard': 'every river routed whole at the base time step, as delineated',
    'substeps': 'rivers too long for dt split into ceil(2kx/dt) equal sub-reaches in series, x unchanged, as '
    'river-route did before this study',
    'substeps-xadj': 'as substeps, with x of each sub-reach lowered so the cascade keeps the reach diffusion',
    'subcycles': 'rivers too short for dt routed in shorter steps of their own (river-route)',
    'stabilized': 'substeps and subcycles together: river-route network_type stabilized before this study',
    'inflate-k': 'k of rivers too short for dt raised to the smallest k that keeps c3 non-negative',
    'merge': 'rivers too short for dt removed, their upstreams and runoff joined to the river downstream',
    'reference': 'dt 60 s with every river subcycled to a step of at most k/5: the time-converged solution',
}
REFERENCE = ('reference', 60)
# the treatments of the comparison of river-route before and after the changes this study made to it
# (scripts/17_river_route_changes.py), which are not cells of the matrix
CHANGE_TREATMENTS = {
    'substeps-fewest': 'rivers too long for dt split into the fewest equal sub-reaches whose coefficients are '
    'non-negative with x lowered as in substeps-xadj: two at x = 0.2',
    'river-route': "river-route's own network_type 'stabilized', as the river-route imported routes it",
}
# the results of that comparison, one folder per version of river-route: before is the code that routed the
# matrix, copied into vendor/river-route-matrix, and after is ../river-route with the changes
CHANGES = DATA / 'river-route-changes'
CODES = {'before': ROOT / 'vendor' / 'river-route-matrix', 'after': ROOT.parent / 'river-route'}


def run_dir(scenario: str, treatment: str, dt: int, root: Path = RESULTS) -> Path:
    """The output directory of one cell of the simulation matrix, or of the comparison under ``root``."""
    if treatment not in TREATMENTS and treatment not in CHANGE_TREATMENTS:
        raise ValueError(f'unknown treatment {treatment!r}, expected one of {sorted(TREATMENTS | CHANGE_TREATMENTS)}')
    if dt <= 0 or DT_RUNOFF % dt != 0:
        raise ValueError(f'dt must be a positive divisor of {DT_RUNOFF}, got {dt}')
    return root / scenario / f'{treatment}__dt{dt:04d}s'
