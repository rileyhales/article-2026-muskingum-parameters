"""
Figures of the theory section: the window of Δt / k with non-negative coefficients, impulse responses on
either side of it, the gain of each reach at the highest frequency a routing step can carry, and the travel time
variance substeps keep or lose.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/06_theory_figures.py
"""


import matplotlib.pyplot as plt
import numpy as np

from stability import config, plotting, theory

X = 0.2  # the x of every river of the hydrofabric
COURANT_EXAMPLES = ((0.1, 'Δt = 0.1k: c1 < 0'), (1.0, 'Δt = k: all ≥ 0'), (8.0, 'Δt = 8k: c3 < 0'))


def plot_window() -> None:
    """The (Δt / k, x) plane split by the sign of c1 and c3, with the x of every reach of the hydrofabric."""
    figure, window = plt.subplots(figsize=(plotting.WIDTH, 4.0))
    courant = np.logspace(-2, 2, 400)
    x = np.linspace(0, 0.5, 200)
    grid_c, grid_x = np.meshgrid(courant, x)
    sign = np.where(grid_c < 2 * grid_x, 0, np.where(grid_c > 2 * (1 - grid_x), 2, 1))
    window.contourf(grid_c, grid_x, sign, levels=[-0.5, 0.5, 1.5, 2.5], colors=['#fbe1d6', '#dbe9fa', '#e2def3'])
    window.plot(2 * x, x, color='#eb6834', linewidth=1.5)
    window.plot(2 * (1 - x), x, color='#4a3aa7', linewidth=1.5)
    window.axhline(X, color=plotting.INK, linewidth=1, linestyle='--')
    window.text(0.012, 0.42, 'c₁ < 0', color=plotting.INK)
    window.text(0.25, 0.015, 'all ≥ 0', color=plotting.INK, ha='center', va='bottom')
    window.text(10, 0.42, 'c₃ < 0', color=plotting.INK)
    window.text(0.012, X + 0.012, 'x = 0.2', color=plotting.INK)
    window.set_xscale('log')
    window.set_xlim(courant[0], courant[-1])
    window.set_xlabel('Δt / k')
    window.set_ylabel('Muskingum x')
    window.grid(False)
    plotting.save(figure, 'theory_coefficient_signs')
    return


def plot_impulse_responses() -> None:
    """Upstream and lateral impulse responses of one reach at a Courant number below, inside, and above the window."""
    figure, axes = plt.subplots(2, 3, figsize=(plotting.WIDTH, 7.0), sharex=True)
    steps = 12
    for column, (courant, title) in enumerate(COURANT_EXAMPLES):
        for row, (response, label) in enumerate(
            ((theory.upstream_impulse_response(courant, X, steps), 'Upstream inflow'),
             (theory.lateral_impulse_response(courant, X, steps), 'Lateral inflow'))
        ):
            axis = axes[row, column]
            colors = np.where(response < 0, '#e34948', '#2a78d6')
            axis.bar(np.arange(steps), response, width=0.6, color=colors)
            axis.axhline(0, color=plotting.MUTED, linewidth=0.6)
            if row == 0:
                axis.set_title(title)
            if column == 0:
                axis.set_ylabel(f'{label}\nresponse per unit')
            if row == 1:
                axis.set_xlabel('Routing step')
    plotting.save(figure, 'theory_impulse_responses')
    return


def plot_gains() -> None:
    """The gain at the Nyquist frequency of upstream and lateral inflow, and the lateral gain across frequencies."""
    figure, (nyquist, spectrum) = plt.subplots(1, 2, figsize=(plotting.WIDTH, 4.0))
    courant = np.logspace(-2, 2, 300)
    upstream, lateral = theory.nyquist_gains(courant, X)
    nyquist.plot(courant, np.abs(upstream), color='#2a78d6', label='Upstream inflow, x/(1−x)')
    nyquist.plot(courant, np.abs(lateral), color='#eb6834', label='Lateral inflow, Δt/(2k(1−x))')
    nyquist.axhline(1, color=plotting.MUTED, linewidth=0.8)
    nyquist.axvline(2 * (1 - X), color='#4a3aa7', linewidth=0.8, linestyle='--')
    nyquist.text(2 * (1 - X) * 0.85, 20, 'c₃ = 0', color=plotting.INK, ha='right')
    nyquist.set_xscale('log')
    nyquist.set_yscale('log')
    nyquist.set_xlabel('Δt / k')
    nyquist.set_ylabel('|gain| at period 2Δt')
    nyquist.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3))  # no quadrant of the panel is clear of the lines
    omega = np.linspace(0, np.pi, 300)
    for courant_value, color in zip((0.5, 1.6, 4.0, 16.0), plotting.DT_RAMP, strict=True):
        _, response = theory.frequency_response(courant_value, X, omega)
        spectrum.plot(omega / (2 * np.pi), np.abs(response), color=color, label=f'Δt = {courant_value:g}k')
    spectrum.set_xlabel('Frequency (cycles per routing step)')
    spectrum.set_ylabel('|gain| of lateral inflow')
    spectrum.set_yscale('log')
    spectrum.legend(loc='upper left')
    plotting.save(figure, 'theory_nyquist_gain')
    return


def route_cascade(inflow: np.ndarray, k: float, x_piece: float, pieces: int, dt: float) -> np.ndarray:
    """Route an inflow series through ``pieces`` equal sub-reaches of x_piece in series, starting at rest."""
    if pieces < 1 or k <= 0:
        raise ValueError('a cascade needs a positive k and at least one piece')
    flow = inflow.copy()
    for _ in range(pieces):
        flow = theory.route_reach(flow, np.zeros(flow.shape[0] - 1), k / pieces, x_piece, dt, flow[0])
    if not np.isclose(flow.sum(), inflow.sum(), rtol=1e-3):
        raise ValueError('the cascade must conserve the volume of a drained pulse')
    return flow


def plot_substeps() -> None:
    """Travel time variance of N substeps with x held and with x adjusted, and a sharp and a broad pulse routed."""
    figure, axes = plt.subplots(1, 3, figsize=(plotting.WIDTH, 4.0))
    pieces = np.arange(1, 21)
    k = 20_000.0
    axes[0].plot(pieces, [theory.cascade_variance(k, X, n) / theory.cascade_variance(k, X, 1) for n in pieces],
                 color=plotting.TREATMENT_COLORS['substeps'], marker='o', markersize=3, label='x held')
    axes[0].plot(pieces, np.ones_like(pieces, dtype=float), color=plotting.TREATMENT_COLORS['substeps-xadj'],
                 label='x adjusted')
    axes[0].set_xlabel('Substeps N')
    axes[0].set_ylabel('Variance / k²(1 − 2x)')
    axes[0].text(6, 0.92, 'x adjusted', color=plotting.INK, va='top')  # direct labels: the panel is too narrow
    axes[0].text(8, 0.2, 'x held', color=plotting.INK, va='bottom')  # for a legend clear of the curves
    dt = 60.0
    hours = np.arange(0, 48 * 3600, dt)
    curves = (('reference', 1, X), ('substeps', 8, X), ('substeps-xadj', 8, float(theory.diffusion_preserving_x(X, 8))))
    for axis, width_h, title in ((axes[1], 0.5, 'Sharp pulse'), (axes[2], 3.0, 'Broad pulse')):
        inflow = 100 * np.exp(-0.5 * ((hours - 14 * 3600) / (width_h * 3600)) ** 2)
        axis.plot(hours / 3600, inflow, color=plotting.MUTED, linewidth=1, label='Inflow')
        for treatment, n_pieces, x_piece in curves:
            label = 'One reach' if n_pieces == 1 else f'{plotting.TREATMENT_LABELS[treatment]}, N = 8'
            axis.plot(hours / 3600, route_cascade(inflow, k, x_piece, n_pieces, dt),
                      color=plotting.TREATMENT_COLORS[treatment], label=label)
        axis.set_title(title)
        axis.set_xlabel('Hour')
        axis.set_xlim(0, 44)
    axes[1].set_ylabel('Discharge (m³/s)')
    handles, labels = axes[2].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2)
    plotting.save(figure, 'theory_substep_diffusion')
    return


if __name__ == '__main__':
    plotting.apply_style()
    plot_window()
    plot_impulse_responses()
    plot_gains()
    plot_substeps()
    print(f'theory figures written to {config.FIGURES}')
