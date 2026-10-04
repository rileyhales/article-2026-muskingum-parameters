"""
The census of reaches outside the window of non-negative Muskingum coefficients: the travel time law of RFS v3, which
gives every reach of every hydrofabric a k from its length and Strahler order, the classes of reaches too long and too
short for a routing time step, and the place of each reach in the network.
"""

import numpy as np

__all__ = ['X', 'POSITIONS', 'velocity', 'travel_time', 'coefficient_signs', 'positions']

X = 0.2  # the Muskingum x of every RFS v3 reach, and of every reach the census assigns a travel time
# the place of a reach in the network, by the confluences at its two ends
POSITIONS = ('between_confluences', 'below_confluence_only', 'above_confluence_only', 'headwater', 'chain')


def velocity(order: np.ndarray) -> np.ndarray:
    """The celerity of the RFS v3 travel time law (Eq. 17), m/s: 0.5 at order 2, changing 0.0375 per order."""
    if order.size == 0:
        raise ValueError('no reaches')
    if order.min() < 1:
        raise ValueError('Strahler orders start at 1')
    return 0.5 + 0.0375 * (order.astype(np.float64) - 2)


def travel_time(length_m: np.ndarray, order: np.ndarray) -> np.ndarray:
    """The travel time k = round(L / v) of Eq. 17, seconds."""
    if length_m.shape != order.shape:
        raise ValueError('length and order must be one per reach')
    if not np.all(np.isfinite(length_m)) or np.any(length_m < 0):
        raise ValueError('lengths must be finite and non-negative')
    return np.round(length_m.astype(np.float64) / velocity(order))


def coefficient_signs(k: np.ndarray, x: np.ndarray | float, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Masks of the reaches too long (c1 < 0, dt < 2kx) and too short (c3 < 0, dt > 2k(1-x)) for dt."""
    if dt <= 0:
        raise ValueError('dt must be positive')
    if np.ndim(x) != 0 and np.shape(x) != k.shape:
        raise ValueError('x must be one value or one per reach')
    return 2 * k * x > dt, 2 * k * (1 - x) < dt


def positions(down: np.ndarray) -> np.ndarray:
    """The place of every reach, as an index into POSITIONS, from the row of its downstream reach (-1 at an outlet).

    A reach whose upstream end is a confluence has two or more inflowing reaches; one whose downstream end is a
    confluence flows into a reach that has. Headwaters have no inflowing reach, chains one and no confluence below.
    """
    n = down.shape[0]
    has_down = down >= 0
    if np.any(down >= n) or np.any(down < -1):
        raise ValueError('downstream rows must be rows of the network or -1')
    if np.any(down[has_down] == np.flatnonzero(has_down)):
        raise ValueError('a reach cannot flow into itself')
    n_up = np.bincount(down[has_down], minlength=n)
    above = np.zeros(n, dtype=bool)
    above[has_down] = n_up[down[has_down]] >= 2
    return np.select(
        [(n_up >= 2) & above, n_up >= 2, n_up == 0, above],
        [POSITIONS.index('between_confluences'), POSITIONS.index('below_confluence_only'),
         POSITIONS.index('headwater'), POSITIONS.index('above_confluence_only')],
        default=POSITIONS.index('chain'),
    )
