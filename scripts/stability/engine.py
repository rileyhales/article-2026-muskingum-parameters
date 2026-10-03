"""
Route one cell of the simulation matrix with river-route's own coefficients and kernels.

``route`` is the routing loop of ``rr.Router`` reduced to what the experiment needs: it builds the layout and
coefficients with ``static_muskingum.prepare_routing`` and routes each chunk of forcing with ``route_region``, carrying
the channel state from chunk to chunk, exactly as the Router does for sequential runoff files. It takes the forcing as
arrays rather than files, so the synthetic forcings never touch the disk, and hands each chunk of routed discharge to a
recorder rather than a writer. ``scripts/verify_engine.py`` checks that it reproduces ``rr.Router`` bit for bit.
"""

import time
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from typing import NamedTuple, Protocol

import numpy as np
from river_route.router import static_muskingum
from river_route.router._numba_kernels import route_region
from river_route.runoff import CatchmentRunoffVolumes

from .treatments import Treated

__all__ = ['Chunk', 'Recorder', 'Timer', 'route', 'steady_state', 'routed_reaches']

MAX_CHUNKS = 100_000  # far more than any run of the matrix holds
HOURS_PER_YEAR = 8760


class Chunk(NamedTuple):
    """One block of forcing: the start of each runoff step and each original river's runoff volume in it (m3)."""

    dates: np.ndarray  # (n_steps,) datetime64[s]
    volumes: np.ndarray  # (n_original_rivers, n_steps) float32


class Recorder(Protocol):
    """Receives the routed discharge of every chunk, a (river, runoff step) array of step means, negatives kept."""

    def record(self, dates: np.ndarray, discharge: np.ndarray) -> None: ...


def routed_reaches(treated: Treated, dt: int, dt_runoff: int) -> tuple[np.ndarray, np.ndarray]:
    """The (n + 1,) reach offsets and (n,) subcycles a treated network is routed with at dt, empty when standard."""
    layout, _ = static_muskingum.prepare_routing(treated.network, treated.network_type, dt, dt_runoff, 'ignore')
    if layout.reach_indptr.shape[0] and layout.reach_indptr.shape[0] != treated.network.size + 1:
        raise ValueError('the layout must have one reach range per river')
    if layout.subcycles.shape[0] and layout.subcycles.shape[0] != treated.network.size:
        raise ValueError('the subcycles must have one count per river')
    return layout.reach_indptr, layout.subcycles


def steady_state(treated: Treated, dt: int, dt_runoff: int, volumes: np.ndarray) -> np.ndarray:
    """
    The channel state of every routed reach at steady state under a constant runoff volume per step: each river
    discharges its own lateral inflow plus everything upstream, and the i-th of its N sub-reaches carries the upstream
    inflow plus i/N of its lateral inflow, since each sub-reach takes 1/N of it.
    """
    if volumes.ndim != 1 or np.any(volumes < 0):
        raise ValueError('steady state needs one non-negative runoff volume per original river')
    lateral = volumes.astype(np.float64) / dt_runoff
    if treated.runoff_map is not None:
        lateral = treated.runoff_map @ lateral
    if lateral.shape[0] != treated.network.size:
        raise ValueError('the runoff must map onto the rivers of the network')
    down = treated.network.downstream_indices
    outflow = lateral.copy()
    for i in range(treated.network.size):  # DFS order puts every river before the river it drains into
        if down[i] >= 0:
            outflow[down[i]] += outflow[i]
    inflow = outflow - lateral
    reach_indptr, _ = routed_reaches(treated, dt, dt_runoff)
    if not reach_indptr.shape[0]:
        return outflow.astype(np.float32)
    pieces = np.diff(reach_indptr)
    river = np.repeat(np.arange(treated.network.size), pieces)
    position = np.arange(reach_indptr[-1]) - np.repeat(reach_indptr[:-1], pieces) + 1
    return (inflow[river] + lateral[river] * position / pieces[river]).astype(np.float32)


class Timer:
    """
    The seconds a run spends in each stage, measured as it goes: preparing each chunk (reading its forcing, mapping it
    onto the network, and allocating its discharge), routing it with river-route's kernels, and recording it. The routing
    seconds of every chunk are kept with its hours, so the cost per simulated year is known from any part of a run.
    """

    def __init__(self) -> None:
        self.seconds = {'preparing': 0.0, 'routing': 0.0, 'recording': 0.0}
        self.chunk_routing_seconds: list[float] = []
        self.chunk_hours: list[int] = []

    def add(self, stage: str, since: float) -> float:
        """Add the time since ``since`` to a stage and return the time now, the start of the next stage."""
        if stage not in self.seconds:
            raise ValueError(f'unknown stage {stage}')
        now = time.perf_counter()
        if now < since:
            raise ValueError('time cannot run backward')
        self.seconds[stage] += now - since
        return now

    def add_chunk(self, since: float, hours: int) -> float:
        """Add the routing time of one chunk of ``hours`` simulated hours, and return the time now."""
        now = self.add('routing', since)
        self.chunk_routing_seconds.append(round(now - since, 4))
        self.chunk_hours.append(int(hours))
        return now

    def summary(self) -> dict:
        """Seconds per stage, simulated hours, routing seconds per simulated year, and the per-chunk record."""
        hours = int(sum(self.chunk_hours))
        if hours <= 0:
            raise ValueError('nothing was routed')
        # the first chunk also loads or compiles the kernels, so the steady rate leaves it out when there is another
        later_hours = int(sum(self.chunk_hours[1:]))
        steady = sum(self.chunk_routing_seconds[1:]) / later_hours if later_hours else self.seconds['routing'] / hours
        return {
            **{f'{stage}_seconds': round(value, 3) for stage, value in self.seconds.items()},
            'simulated_hours': hours,
            'routing_seconds_per_simulated_year': round(self.seconds['routing'] * HOURS_PER_YEAR / hours, 3),
            'steady_routing_seconds_per_simulated_year': round(steady * HOURS_PER_YEAR, 4),
            'chunk_routing_seconds': self.chunk_routing_seconds,
            'chunk_hours': self.chunk_hours,
        }


def route(
    treated: Treated,
    dt: int,
    dt_runoff: int,
    chunks: Iterable[Chunk],
    recorder: Recorder,
    threads: int = 1,
    initial_state: np.ndarray | None = None,
) -> dict:
    """
    Route every chunk of forcing in order over a treated network at routing step dt, handing each chunk's discharge to
    the recorder. Returns a summary of the run: the reaches and subcycles routed, and the time each stage took.
    """
    if dt <= 0 or dt_runoff % dt != 0:
        raise ValueError(f'dt={dt} must be a positive divisor of dt_runoff={dt_runoff}')
    if threads < 1:
        raise ValueError(f'threads must be at least 1, got {threads}')
    network = treated.network
    layout, method = static_muskingum.prepare_routing(network, treated.network_type, dt, dt_runoff, 'ignore')
    n_states = int(layout.reach_indptr[-1]) if layout.reach_indptr.shape[0] else network.size
    q_t = np.zeros(n_states, dtype=np.float32) if initial_state is None else initial_state.astype(np.float32)
    if q_t.shape != (n_states,):
        raise ValueError(f'the initial state has {q_t.shape[0]} values for {n_states} routed reaches')
    n_out = network.original_river_ids.shape[0]
    timer = Timer()
    started = time.perf_counter()
    remaining = iter(chunks)
    with ThreadPoolExecutor(threads) as pool:
        for _ in range(MAX_CHUNKS):
            mark = time.perf_counter()
            chunk = next(remaining, None)
            if chunk is None:
                break
            volumes = chunk.volumes if treated.runoff_map is None else treated.runoff_map @ chunk.volumes
            runoff = CatchmentRunoffVolumes(np.ascontiguousarray(volumes, dtype=np.float32), network.river_ids)
            discharge = np.zeros((n_out, chunk.dates.shape[0]), dtype=np.float32)
            mark = timer.add('preparing', mark)
            route_region(network, method, layout, dt_runoff // dt, q_t, discharge, runoff, pool, threads)
            mark = timer.add_chunk(mark, chunk.dates.shape[0])
            recorder.record(chunk.dates, discharge)
            timer.add('recording', mark)
        else:
            raise ValueError(f'a run may hold at most {MAX_CHUNKS} chunks')
    subcycles = layout.subcycles if layout.subcycles.shape[0] else np.ones(network.size, dtype=np.int64)
    pieces = np.diff(layout.reach_indptr) if layout.reach_indptr.shape[0] else np.ones(network.size, dtype=np.int64)
    return {
        'rivers': int(network.size),
        'reaches': int(n_states),
        'subcycled_rivers': int(np.count_nonzero(subcycles > 1)),
        'reach_steps_per_routing_step': int(np.sum(pieces * subcycles)),
        'wall_seconds': round(time.perf_counter() - started, 3),
        **timer.summary(),
        'final_state_sum': float(q_t.astype(np.float64).sum()),
    }
