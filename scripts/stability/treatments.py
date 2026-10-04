"""
The treatments of rivers whose Muskingum coefficients are negative at a routing time step.

Every treatment is a river-route Network routed with one of its two network types. The treatments that change only
how a river is routed (substeps, subcycles, or both) are a ``ConditionedNetwork``, whose ``conditioning`` gives the
substeps and subcycles of the treatment rather than of river-route's own analysis, routed as a stabilized network.
The river-route treatment is river-route's own stabilized network, whichever version of river-route is imported.
The treatments that edit the network table (inflate-k, merge) are an ordinary Network of the edited table, routed as a
standard network, and merge also maps the runoff of each removed river onto the river that receives it.
"""

from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
import river_route as rr
import scipy.sparse

from . import config

__all__ = [
    'ConditionedNetwork', 'Treated', 'treat', 'merge_short_rivers', 'inflate_short_rivers', 'matrix_substeps',
    'fewest_substeps',
]

REFERENCE_COURANT = 0.2  # the reference routes every river at a step of at most this fraction of its k


def matrix_substeps(k: np.ndarray, x: np.ndarray, dt: int) -> np.ndarray:
    """
    The substeps of the matrix: ceil(2kx/dt) equal sub-reaches for every river too long for dt, the fewest whose
    coefficients are all non-negative when they keep the river's x, as river-route counted them before this study. A
    river whose sub-reaches would be too short for dt is left whole, as river-route left it.
    """
    if dt <= 0 or k.shape != x.shape:
        raise ValueError('substeps need a positive dt and one x per k')
    count = np.maximum(1, np.ceil(2 * k * x / dt))
    keeps_c3 = k / count >= dt / (2 * (1 - x))
    return np.where((2 * k * x > dt) & keeps_c3, count, 1).astype(np.int64)


def fewest_substeps(k: np.ndarray, x: np.ndarray, dt: int) -> np.ndarray:
    """
    The fewest equal sub-reaches with x_N = 1/2 - N (1/2 - x) whose coefficients are all non-negative at dt, as
    river-route counts them after this study: ceil(k / (dt + (1 - 2x) k)), two for every river too long at x = 0.2.
    """
    if dt <= 0 or k.shape != x.shape:
        raise ValueError('substeps need a positive dt and one x per k')
    count = np.maximum(1, np.ceil(k / (dt + (1 - 2 * x) * k)))
    keeps_c3 = k / count >= dt - (1 - 2 * x) * k
    return np.where((2 * k * x > dt) & keeps_c3, count, 1).astype(np.int64)


SUBSTEPS = {  # the substeps of each treatment that divides rivers too long for dt, and whether it adjusts x
    'substeps': (matrix_substeps, False),
    'substeps-xadj': (matrix_substeps, True),
    'stabilized': (matrix_substeps, False),
    'substeps-fewest': (fewest_substeps, True),
}


class ConditionedNetwork(rr.Network):
    """
    A Network whose substeps, subcycles, and x at one routing time step are chosen by a treatment rather than by the
    river-route it is routed with. ``prepare_routing`` reads them through ``conditioning``, and through ``x`` or
    ``substep_x``, when the network is routed with network_type 'stabilized'; river-route before this study reads the
    x of every sub-reach from ``x``, and after it from ``substep_x``, so the treatment routes the same in both.

    A treatment that adjusts x gives the sub-reaches of a river divided into N the x_N = 1/2 - N (1/2 - x), which keeps
    the travel time variance of the river, k^2 (1 - 2x), that equal sub-reaches of the same x would divide by N.

    Rivers in the optional ``untreated`` mask are routed whole at the base time step, which isolates the effect of
    leaving those rivers untreated: routing is linear, so the difference from the fully treated network downstream is
    exactly what they contribute.
    """

    def __init__(
        self, network_file: Path | pd.DataFrame, treatment: str, dt: int, untreated: np.ndarray | None = None
    ) -> None:
        self._x_routed = None  # read by the x property, so set before the Network reads its table
        super().__init__(network_file)
        if treatment not in (*SUBSTEPS, 'subcycles', 'reference'):
            raise ValueError(f'{treatment} is not a treatment of substeps or subcycles')
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        self.treatment = treatment
        self.dt = dt
        k, x = self.k.astype(np.float64), self.x.astype(np.float64)
        ones = np.ones(self.size, dtype=np.int64)
        _, too_short = rr.Network.unstable_mask(self, float(dt))
        subcycles = np.where(too_short, rr.Network.subcycles(self, float(dt))[0], 1)
        if treatment not in ('subcycles', 'stabilized'):
            subcycles = ones
        substeps = SUBSTEPS[treatment][0](k, x, dt) if treatment in SUBSTEPS else ones
        if treatment == 'reference':
            own_step_limit = REFERENCE_COURANT * k
            subcycles = np.maximum(1, np.ceil(dt / own_step_limit)).astype(np.int64)
        if untreated is not None:
            if untreated.shape != (self.size,) or untreated.dtype != bool:
                raise ValueError('untreated must be a boolean mask with one value per river')
            substeps = np.where(untreated, 1, substeps)
            subcycles = np.where(untreated, 1, subcycles)
        self._substeps, self._subcycles = substeps, subcycles
        if treatment in SUBSTEPS and SUBSTEPS[treatment][1]:
            self._x_routed = (0.5 - substeps * (0.5 - x)).astype(np.float32)
        return

    @property
    def x(self) -> np.ndarray:
        """(n,) Muskingum x each river is routed with: the x of every sub-reach for a treatment that adjusts x."""
        return self._df['muskingumX'].to_numpy() if self._x_routed is None else self._x_routed

    def substep_x(self, substeps: np.ndarray) -> np.ndarray:
        """The x of every sub-reach is the treatment's, which ``x`` already gives: the river's own, or adjusted."""
        if substeps.shape != (self.size,) or not np.array_equal(substeps, self._substeps):
            raise ValueError('the sub-reach x of a treatment is for its own substeps')
        return self.x.astype(np.float64)

    def conditioning(self, dt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The substeps and subcycles of the treatment, at the one time step it was built for."""
        if int(dt) != self.dt:
            raise ValueError(f'this network was conditioned for dt={self.dt}, not {dt}')
        if self._substeps.shape[0] != self.size:
            raise ValueError('the substeps must have one count per river')
        return self._substeps, self._subcycles, np.ones(self.size, dtype=bool)


class Treated(NamedTuple):
    """A treated network, the network type it is routed as, and the runoff map of a merge (None otherwise)."""

    network: rr.Network
    network_type: str
    runoff_map: scipy.sparse.csr_matrix | None  # (n_routed, n_original) sums each original river's runoff onto a row
    kept: np.ndarray  # (n_original,) True for each original river with a row of routed discharge


def inflate_short_rivers(table: pd.DataFrame, dt: int) -> pd.DataFrame:
    """The network table with k of each river too short for dt raised to the smallest k at which c3 is not negative."""
    if dt <= 0:
        raise ValueError(f'dt must be positive, got {dt}')
    k = table['muskingumK'].to_numpy(dtype=np.float64)
    x = table['muskingumX'].to_numpy(dtype=np.float64)
    k_min = dt / (2 * (1 - x)) * (1 + 1e-6)  # a hair above the bound so float32 rounding cannot leave c3 negative
    inflated = table.copy()
    inflated['muskingumK'] = np.maximum(k, k_min).astype(np.float32)
    if np.any(inflated['muskingumK'].to_numpy() < table['muskingumK'].to_numpy()):
        raise ValueError('inflating k must never shorten a river')
    return inflated


def merge_short_rivers(table: pd.DataFrame, dt: int) -> tuple[pd.DataFrame, scipy.sparse.csr_matrix, np.ndarray]:
    """
    Remove every river too short for dt that is not a basin outlet. The rivers upstream of a removed river drain into
    the first kept river downstream of it, which also receives its runoff; at a confluence this joins the confluences
    at either end of the short river into one. Deleting rows of a DFS-ordered table leaves it in DFS order, so only
    upstreamCount and riverIndex are renumbered.

    Returns:
        tuple: (the merged table, the (n_kept, n_rivers) runoff map, the (n_rivers,) mask of kept rivers)
    """
    if dt <= 0:
        raise ValueError(f'dt must be positive, got {dt}')
    network = rr.Network(table)
    _, too_short = network.unstable_mask(float(dt))
    kept = ~too_short | (network.downstream_indices < 0)
    receiver = np.arange(network.size)
    for i in range(network.size - 1, -1, -1):  # downstream before upstream: a removed river's receiver is resolved
        d = network.downstream_indices[i]
        receiver[i] = i if kept[i] else receiver[d]
    if not np.all(kept[receiver]):
        raise ValueError('every river must drain into a kept river')

    new_row = np.cumsum(kept) - 1
    removed_before = np.concatenate(([0], np.cumsum(~kept)))
    upstream_counts = network.upstream_counts.astype(np.int64)
    rows = np.arange(network.size)
    removed_upstream = removed_before[rows] - removed_before[rows - upstream_counts]
    merged = table.copy()
    merged['upstreamCount'] = (upstream_counts - removed_upstream).astype(np.int32)
    down = network.downstream_indices
    next_ids = np.where(down >= 0, network.river_ids[receiver[np.clip(down, 0, None)]], -1)
    merged['nextRiverId'] = next_ids
    merged = merged[kept].reset_index(drop=True)
    merged['riverIndex'] = (int(table['riverIndex'].iloc[0]) + np.arange(len(merged))).astype(np.int32)

    runoff_map = scipy.sparse.csr_matrix(
        (np.ones(network.size, dtype=np.float32), (new_row[receiver], rows)), shape=(len(merged), network.size)
    )
    return merged, runoff_map, kept


def treat(treatment: str, dt: int, table: pd.DataFrame) -> Treated:
    """Build the network of a treatment at a routing time step from the original network table."""
    if treatment not in config.TREATMENTS and treatment not in config.CHANGE_TREATMENTS:
        raise ValueError(f'unknown treatment {treatment!r}')
    if dt <= 0:
        raise ValueError(f'dt must be positive, got {dt}')
    every_river = np.ones(len(table), dtype=bool)
    if treatment == 'standard':
        return Treated(rr.Network(table), 'standard', None, every_river)
    if treatment == 'river-route':
        return Treated(rr.Network(table), 'stabilized', None, every_river)
    if treatment == 'inflate-k':
        return Treated(rr.Network(inflate_short_rivers(table, dt)), 'standard', None, every_river)
    if treatment == 'merge':
        merged, runoff_map, kept = merge_short_rivers(table, dt)
        return Treated(rr.Network(merged), 'standard', runoff_map, kept)
    return Treated(ConditionedNetwork(table, treatment, dt), 'stabilized', None, every_river)
