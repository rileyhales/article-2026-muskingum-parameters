"""
Analytic properties of the Muskingum scheme in float64, written in the Courant number C = dt / k.

The coefficients of Q[n+1] = c1 I[n+1] + c2 I[n] + c3 Q[n] + c4 L[n] are rational in C and x, so the impulse responses,
frequency responses, and travel time moments of a reach follow from them directly. These functions are the theory the
manuscript derives; ``route_reach`` repeats river-route's recurrence for one reach in float64, an independent check of
the discharge river-route writes.
"""

import numpy as np

__all__ = [
    'coefficients',
    'positivity_window',
    'upstream_impulse_response',
    'lateral_impulse_response',
    'frequency_response',
    'nyquist_gains',
    'response_moments',
    'cascade_variance',
    'diffusion_preserving_x',
    'route_reach',
]


def coefficients(courant: np.ndarray | float, x: np.ndarray | float) -> tuple[np.ndarray, ...]:
    """The Muskingum coefficients c1, c2, c3, and the lateral inflow coefficient c4 = c1 + c2."""
    courant = np.asarray(courant, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    if np.any(courant <= 0):
        raise ValueError('the Courant number dt/k must be positive')
    if np.any(x > 0.5):
        raise ValueError('x must not exceed 1/2')
    denominator = courant + 2 * (1 - x)
    c1 = (courant - 2 * x) / denominator
    c2 = (courant + 2 * x) / denominator
    c3 = (2 * (1 - x) - courant) / denominator
    return c1, c2, c3, c1 + c2


def positivity_window(x: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """The range of Courant numbers 2x <= C <= 2(1 - x) over which every coefficient is non-negative."""
    x = np.asarray(x, dtype=np.float64)
    if np.any(x > 0.5):
        raise ValueError('x must not exceed 1/2')
    lower, upper = 2 * x, 2 * (1 - x)
    if np.any(upper < lower):
        raise ValueError('the window must not be empty')
    return lower, upper


def upstream_impulse_response(courant: float, x: float, n: int) -> np.ndarray:
    """Discharge at levels 0..n-1 after a unit pulse of upstream inflow at level 0: c1, then (c2 + c1 c3) c3^(j-1)."""
    if n < 2:
        raise ValueError('the response needs at least two levels')
    c1, c2, c3, _ = coefficients(courant, x)
    response = np.empty(n)
    response[0] = c1
    response[1:] = (c2 + c1 * c3) * c3 ** np.arange(n - 1)
    if not np.isclose(c1 + (c2 + c1 * c3) / (1 - c3), 1):
        raise ValueError('the whole response must hold the unit volume')
    return response


def lateral_impulse_response(courant: float, x: float, n: int) -> np.ndarray:
    """Discharge at levels 0..n-1 after lateral inflow of one unit during step 0 only: 0, then c4 c3^(j-1)."""
    if n < 2:
        raise ValueError('the response needs at least two levels')
    _, _, c3, c4 = coefficients(courant, x)
    response = np.zeros(n)
    response[1:] = c4 * c3 ** np.arange(n - 1)
    if np.any(~np.isfinite(response)):
        raise ValueError('the response must be finite')
    return response


def frequency_response(courant: float, x: float, omega: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    The complex gain of a reach at angular frequencies omega (radians per step) for upstream inflow,
    (c1 + c2 z^-1) / (1 - c3 z^-1), and for lateral inflow held over each step, c4 z^-1 / (1 - c3 z^-1).
    """
    omega = np.asarray(omega, dtype=np.float64)
    if np.any(omega < 0) or np.any(omega > np.pi):
        raise ValueError('frequencies must lie in [0, pi]')
    c1, c2, c3, c4 = coefficients(courant, x)
    z_inverse = np.exp(-1j * omega)
    upstream = (c1 + c2 * z_inverse) / (1 - c3 * z_inverse)
    lateral = c4 * z_inverse / (1 - c3 * z_inverse)
    return upstream, lateral


def nyquist_gains(courant: np.ndarray | float, x: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """
    The gains at the Nyquist frequency, a series alternating every step: -x/(1-x) for upstream inflow at every C,
    and -C/(2(1-x)) for lateral inflow, whose magnitude exceeds 1 exactly when c3 < 0.
    """
    courant = np.asarray(courant, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    if np.any(courant <= 0) or np.any(x >= 1):
        raise ValueError('C must be positive and x below 1')
    upstream = np.broadcast_to(-x / (1 - x), np.broadcast(courant, x).shape).astype(np.float64)
    lateral = -courant / (2 * (1 - x))
    return upstream, lateral


def response_moments(response: np.ndarray, dt: float) -> tuple[float, float]:
    """The mean and variance, in seconds and seconds squared, of a discrete response read at levels j * dt."""
    if response.ndim != 1 or dt <= 0:
        raise ValueError('a response is a 1D series at a positive time step')
    volume = response.sum()
    if not np.isfinite(volume) or abs(volume) < 1e-12:
        raise ValueError('the response must hold volume')
    times = np.arange(response.shape[0]) * dt
    mean = float((times * response).sum() / volume)
    variance = float(((times - mean) ** 2 * response).sum() / volume)
    return mean, variance


def cascade_variance(k: float, x: float, pieces: int) -> float:
    """Travel time variance of a reach of travel time k split into equal pieces of the same x: k^2 (1 - 2x) / N."""
    if k <= 0 or pieces < 1:
        raise ValueError('k must be positive and pieces at least one')
    if x > 0.5:
        raise ValueError('x must not exceed 1/2')
    return k**2 * (1 - 2 * x) / pieces


def diffusion_preserving_x(x: np.ndarray | float, pieces: np.ndarray | int) -> np.ndarray:
    """The x of each of N equal pieces whose cascade keeps the variance k^2 (1 - 2x) of the whole: 1/2 - N (1/2 - x)."""
    x = np.asarray(x, dtype=np.float64)
    pieces = np.asarray(pieces)
    if np.any(pieces < 1) or np.any(x > 0.5):
        raise ValueError('pieces must be at least one and x at most 1/2')
    return 0.5 - pieces * (0.5 - x)


def route_reach(
    inflow: np.ndarray, lateral: np.ndarray, k: float, x: float, dt: float, q0: float, steps_per_value: int = 1
) -> np.ndarray:
    """
    Route one reach in float64 exactly as river-route's recurrence does, negative discharge included.
    ``inflow`` holds the upstream discharge at the levels 0..n of the routing steps, ``lateral`` the lateral inflow
    rate held over each of the n // steps_per_value runoff steps. Returns the discharge at levels 0..n.
    """
    n = inflow.shape[0] - 1
    if lateral.shape[0] * steps_per_value != n:
        raise ValueError('each lateral value must span steps_per_value routing steps')
    if k <= 0 or dt <= 0:
        raise ValueError('k and dt must be positive')
    c1, c2, c3, c4 = (float(c) for c in coefficients(dt / k, x))
    forcing = c1 * inflow[1:] + c2 * inflow[:-1] + c4 * np.repeat(lateral, steps_per_value)
    discharge = np.empty(n + 1)
    discharge[0] = q0
    for step in range(n):
        discharge[step + 1] = c3 * discharge[step] + forcing[step]
    return discharge
