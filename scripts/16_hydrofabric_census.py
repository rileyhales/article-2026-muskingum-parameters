"""
Census of reaches outside the window of non-negative Muskingum coefficients in the hydrofabrics of large-scale routing:
TDX-Hydro as delineated, HydroRIVERS, MERIT-Basins, NHDPlus V2, and NHDPlus HR, and, in the tables only, TDX-Hydro as
processed for RFS v3.
Every reach of every hydrofabric is given a travel time by the RFS v3 velocity law (Eq. 17) from its length and Strahler
order, and x = 0.2, so the census compares how each segmentation places reach lengths relative to the window,
independently of each model's own parameters. The figures show the hydrofabrics as published, not routing networks
derived from them.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/16_hydrofabric_census.py
"""

from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyogrio

from stability import census, config, plotting

TDX_DIR = Path.home() / 'data' / 'TDXHydroGeoParquet'  # the TDX-Hydro stream networks, one file per region
TDX_REGIONS = 62
RFS_METADATA = config.RFS_ROOT / 'hydrography' / 'global' / 'metadata.parquet'
HYDROGRAPHY = config.ROOT / 'hydrography'
HYDRORIVERS_FILE = HYDROGRAPHY / 'HydroRIVERS_v10.gdb' / 'HydroRIVERS_v10.gdb'
MERIT_DIR = HYDROGRAPHY / 'merit-hydro-rivers' / 'MERIT-Hydro_v07_Basins_v01_bugfix1' / 'pfaf_level_01'
MERIT_REGIONS = 9  # the Pfafstetter level 1 basins that together cover the globe
NHDPLUS_V2_FILE = HYDROGRAPHY / 'NHDPlusV21_National_Seamless_Flattened_Lower48.gdb'
NHDPLUS_V2_COASTLINE = 'Coastline'  # the NHDPlus V2 feature type of coastlines
NHDPLUS_HR_FILE = HYDROGRAPHY / 'NHDPlus_H_National_Release_2_GDB.gdb'
NHD_COASTLINE = 566  # the NHD feature type of coastlines, network flowlines that are not channels and have no order
NHD_NO_DATA = -9998  # the value NHDPlus HR gives a missing attribute
DT_CENSUS = (30, 60, 300, 600, 900, 1800, 3600, 7200, 10800)  # routing time steps, seconds
QUANTILES = (0.01, 0.10, 0.50, 0.90, 0.99)  # of reach length and travel time
K_WINDOWS = (300, 3600)  # the steps whose window of k the travel time figure shades


def downstream_rows(ids: np.ndarray, next_ids: np.ndarray, outlet: np.ndarray) -> np.ndarray:
    """The row of the reach downstream of each reach, -1 at outlets."""
    index = pd.Index(ids)
    if not index.is_unique:
        raise ValueError('reach ids must be unique')
    if not ids.shape == next_ids.shape == outlet.shape:
        raise ValueError('ids, downstream ids, and outlets must be one per reach')
    down = np.full(ids.shape[0], -1, dtype=np.int64)
    down[~outlet] = index.get_indexer(next_ids[~outlet])
    if np.any(down[~outlet] < 0):
        raise ValueError('a downstream id is not a reach of the network')
    return down


def network(length_m: np.ndarray, order: np.ndarray, down: np.ndarray) -> pd.DataFrame:
    """The table the census reads: length, Strahler order, and downstream row of every reach."""
    if not length_m.shape == order.shape == down.shape:
        raise ValueError('length, order, and downstream row must be one per reach')
    if not np.all(np.isfinite(length_m)) or np.any(length_m < 0):
        raise ValueError('lengths must be finite and non-negative')
    return pd.DataFrame({
        'length_m': length_m.astype(np.float64), 'order': order.astype(np.int64), 'down': down.astype(np.int64),
    })


def load_tdxhydro() -> pd.DataFrame:
    """Every stream link of TDX-Hydro, as delineated, with the length and order of its stream network file."""
    files = sorted(TDX_DIR.glob('TDX_streamnet_*_01.parquet'))
    if len(files) != TDX_REGIONS:
        raise FileNotFoundError(f'expected {TDX_REGIONS} TDX-Hydro regions in {TDX_DIR}, found {len(files)}')
    columns = ['LINKNO', 'DSLINKNO', 'strmOrder', 'Length', 'TDXHydroLinkNo']
    links = pd.concat([pd.read_parquet(file, columns=columns) for file in files], ignore_index=True)
    # LINKNO and DSLINKNO number the links of one region; TDXHydroLinkNo adds the region's prefix to LINKNO
    prefix = links['TDXHydroLinkNo'].to_numpy() - links['LINKNO'].to_numpy()
    outlet = links['DSLINKNO'].to_numpy() < 0
    next_ids = np.where(outlet, -1, links['DSLINKNO'].to_numpy() + prefix)
    down = downstream_rows(links['TDXHydroLinkNo'].to_numpy(), next_ids, outlet)
    return network(links['Length'].to_numpy(), links['strmOrder'].to_numpy(), down)


def load_rfs_v3() -> pd.DataFrame:
    """Every reach of the RFS v3 network, which carries its own k; Eq. 17 must reproduce it to within a second."""
    reaches = pd.read_parquet(RFS_METADATA, columns=['riverId', 'nextRiverId', 'strahlerOrder', 'Length', 'musk_k'])
    outlet = reaches['nextRiverId'].to_numpy() == -1
    down = downstream_rows(reaches['riverId'].to_numpy(), reaches['nextRiverId'].to_numpy(), outlet)
    table = network(reaches['Length'].to_numpy(), reaches['strahlerOrder'].to_numpy(), down)
    k = census.travel_time(table['length_m'].to_numpy(), table['order'].to_numpy())
    if np.max(np.abs(k - reaches['musk_k'].to_numpy())) > 1:
        raise ValueError('Eq. 17 does not reproduce the k of RFS v3')
    return table


def load_hydrorivers() -> pd.DataFrame:
    """Every reach of HydroRIVERS v1.0."""
    reaches = pyogrio.read_dataframe(
        HYDRORIVERS_FILE, columns=['HYRIV_ID', 'NEXT_DOWN', 'LENGTH_KM', 'ORD_STRA'], read_geometry=False,
    )
    # LENGTH_KM is published to 0.01 km and stored as float32; rounding to 0.01 km removes the float32 error
    stored_km = reaches['LENGTH_KM'].to_numpy(np.float64)
    length_km = np.round(stored_km, 2)
    if np.max(np.abs(length_km - stored_km)) > 1e-4:
        raise ValueError('HydroRIVERS lengths are not given to 0.01 km')
    outlet = reaches['NEXT_DOWN'].to_numpy() == 0
    down = downstream_rows(reaches['HYRIV_ID'].to_numpy(), reaches['NEXT_DOWN'].to_numpy(), outlet)
    return network(length_km * 1000, reaches['ORD_STRA'].to_numpy(), down)


def load_merit_basins() -> pd.DataFrame:
    """Every river reach of MERIT-Basins (MERIT Hydro v0.7, Basins v1, bugfix 1), from its nine level 1 basins."""
    files = sorted(MERIT_DIR.glob('riv_pfaf_*_MERIT_Hydro_v07_Basins_v01_bugfix1.shp'))
    if len(files) != MERIT_REGIONS:
        raise FileNotFoundError(f'expected {MERIT_REGIONS} MERIT-Basins river files in {MERIT_DIR}, found {len(files)}')
    columns = ['COMID', 'NextDownID', 'lengthkm', 'order']
    reaches = pd.concat(
        [pyogrio.read_dataframe(file, columns=columns, read_geometry=False) for file in files], ignore_index=True,
    )
    # some first-order reaches of zero or near-zero length flow into themselves; nothing flows into them, so they are
    # kept as isolated reaches, each its own outlet
    self_loop = reaches['NextDownID'].to_numpy() == reaches['COMID'].to_numpy()
    if np.any(reaches['lengthkm'].to_numpy()[self_loop] > 0.001):
        raise ValueError('a MERIT-Basins reach longer than 1 m flows into itself')
    if np.any(np.isin(reaches['NextDownID'].to_numpy()[~self_loop], reaches['COMID'].to_numpy()[self_loop])):
        raise ValueError('a MERIT-Basins reach flows into a reach that flows into itself')
    print(f'MERIT-Basins: {self_loop.sum()} reaches flow into themselves and are kept as isolated reaches')
    outlet = (reaches['NextDownID'].to_numpy() == 0) | self_loop
    down = downstream_rows(reaches['COMID'].to_numpy(), reaches['NextDownID'].to_numpy(), outlet)
    return network(reaches['lengthkm'].to_numpy(np.float64) * 1000, reaches['order'].to_numpy(), down)


def load_nhdplus_v2() -> pd.DataFrame:
    """Every network flowline of NHDPlus Version 2.1 (the conterminous United States) but the coastlines."""
    columns = ['Hydroseq', 'DnHydroseq', 'TerminalFl', 'FTYPE', 'StreamOrde', 'LENGTHKM']
    flowlines = pyogrio.read_dataframe(
        NHDPLUS_V2_FILE, layer='NHDFlowline_Network', columns=columns, read_geometry=False,
    )
    coast = flowlines['FTYPE'].to_numpy() == NHDPLUS_V2_COASTLINE
    reaches = flowlines[~coast]
    if np.any(reaches['StreamOrde'].to_numpy() < 1):
        raise ValueError('an NHDPlus V2 flowline that is not a coastline has no stream order')
    # a flowline drains along the main path, to the flowline whose Hydroseq is its DnHydroseq, so a minor divergence
    # begins with no inflowing flowline; DnHydroseq is 0 at the outlet of a network, and flowlines that drain to a
    # coastline end there, as the terminal flag of each says
    next_ids = reaches['DnHydroseq'].to_numpy()
    to_coast = np.isin(next_ids, flowlines['Hydroseq'].to_numpy()[coast])
    outlet = (next_ids == 0) | to_coast
    if not np.array_equal(outlet, reaches['TerminalFl'].to_numpy() == 1):
        raise ValueError('the outlets of NHDPlus V2 must be the flowlines its terminal flag marks')
    print(f'NHDPlus V2: {coast.sum()} coastline flowlines left out, {to_coast.sum()} flowlines that drain to them kept '
          f'as outlets')
    down = downstream_rows(reaches['Hydroseq'].to_numpy(), next_ids, outlet)
    return network(reaches['LENGTHKM'].to_numpy() * 1000, reaches['StreamOrde'].to_numpy(), down)


def load_nhdplus_hr() -> pd.DataFrame:
    """Every network flowline of NHDPlus High Resolution (national release 2) but the coastlines."""
    columns = ['hydroseq', 'dnhydroseq', 'terminalpa', 'ftype', 'streamorde', 'lengthkm']
    flowlines = pyogrio.read_dataframe(
        NHDPLUS_HR_FILE, layer='NetworkNHDFlowline', columns=columns, read_geometry=False,
    )
    coast = flowlines['ftype'].to_numpy() == NHD_COASTLINE
    reaches = flowlines[~coast]
    if np.any(reaches['streamorde'].to_numpy() < 1):
        raise ValueError('an NHDPlus HR flowline that is not a coastline has no stream order')
    # a flowline drains along the main path, to the flowline whose hydroseq is its dnhydroseq, so a minor divergence
    # begins with no inflowing flowline; dnhydroseq is 0 at the outlet of a network, missing or -9998 at a few outlets,
    # and flowlines that drain to a coastline end there
    missing = reaches['dnhydroseq'].isna().to_numpy() | (reaches['dnhydroseq'].to_numpy() == NHD_NO_DATA)
    if np.any(reaches['terminalpa'].to_numpy()[missing] != reaches['hydroseq'].to_numpy()[missing]):
        raise ValueError('an NHDPlus HR flowline with no dnhydroseq is not the outlet of its network')
    next_ids = np.where(missing, 0, reaches['dnhydroseq'].to_numpy())
    to_coast = np.isin(next_ids, flowlines['hydroseq'].to_numpy()[coast])
    print(f'NHDPlus HR: {coast.sum()} coastline flowlines left out, {to_coast.sum()} flowlines that drain to them and '
          f'{missing.sum()} with no dnhydroseq kept as outlets')
    outlet = (next_ids == 0) | to_coast
    down = downstream_rows(reaches['hydroseq'].to_numpy(), next_ids, outlet)
    return network(reaches['lengthkm'].to_numpy() * 1000, reaches['streamorde'].to_numpy(), down)


HYDROFABRICS: dict[str, Callable[[], pd.DataFrame]] = {
    'TDX-Hydro': load_tdxhydro,
    'RFS v3': load_rfs_v3,
    'HydroRIVERS': load_hydrorivers,
    'MERIT-Basins': load_merit_basins,
    'NHDPlus V2': load_nhdplus_v2,
    'NHDPlus HR': load_nhdplus_hr,
}
# the hydrofabrics as published
FIGURE_HYDROFABRICS = ('TDX-Hydro', 'HydroRIVERS', 'MERIT-Basins', 'NHDPlus V2', 'NHDPlus HR')


def sign_table(name: str, k: np.ndarray) -> pd.DataFrame:
    """Counts and shares of reaches too long, too short, and inside the window at each routing time step."""
    rows = []
    for dt in DT_CENSUS:
        too_long, too_short = census.coefficient_signs(k, census.X, dt)
        rows.append({
            'hydrofabric': name, 'dt_s': dt, 'reaches': k.shape[0], 'too_long': int(too_long.sum()),
            'too_short': int(too_short.sum()), 'valid': int((~too_long & ~too_short).sum()),
        })
    table = pd.DataFrame(rows)
    if not np.all(table[['too_long', 'too_short', 'valid']].sum(axis=1) == table['reaches']):
        raise ValueError('every reach must fall in exactly one class')
    for column in ('too_long', 'too_short', 'valid'):
        table[f'{column}_pct'] = 100 * table[column] / table['reaches']
    return table


def length_summary(name: str, table: pd.DataFrame, k: np.ndarray) -> pd.DataFrame:
    """Reach count, zero-length reaches, and the distribution of reach length and travel time."""
    length = table['length_m'].to_numpy()
    if k.shape != length.shape:
        raise ValueError('k must be one per reach')
    row = {
        'hydrofabric': name, 'reaches': length.shape[0], 'zero_length': int((length == 0).sum()),
        'order_min': int(table['order'].min()), 'order_max': int(table['order'].max()),
        'length_mean_m': float(length.mean()), 'length_total_km': float(length.sum() / 1000),
    }
    quantiles = zip(QUANTILES, np.quantile(length, QUANTILES), np.quantile(k, QUANTILES))
    for q, length_q, k_q in quantiles:
        row[f'length_p{100 * q:02.0f}_m'] = float(length_q)
        row[f'k_p{100 * q:02.0f}_s'] = float(k_q)
    return pd.DataFrame([row])


def order_table(name: str, table: pd.DataFrame, k: np.ndarray) -> pd.DataFrame:
    """Reaches, celerity, length, and travel time of each Strahler order."""
    frame = table.assign(k_s=k, velocity_m_s=census.velocity(table['order'].to_numpy()))
    summary = frame.groupby('order').agg(
        reaches=('length_m', 'size'), velocity_m_s=('velocity_m_s', 'first'),
        length_median_m=('length_m', 'median'), length_mean_m=('length_m', 'mean'), k_median_s=('k_s', 'median'),
    ).reset_index()
    if summary['reaches'].sum() != len(table):
        raise ValueError('every reach must have an order')
    return summary.assign(hydrofabric=name)[['hydrofabric', *summary.columns]]


def topology_table(name: str, table: pd.DataFrame, k: np.ndarray) -> pd.DataFrame:
    """Where the reaches too short for each step sit in the network, by the confluences at their two ends."""
    place = census.positions(table['down'].to_numpy())
    length = table['length_m'].to_numpy()
    rows = []
    for dt in DT_CENSUS:
        _, short = census.coefficient_signs(k, census.X, dt)
        counts = np.bincount(place[short], minlength=len(census.POSITIONS))
        rows.append({
            'hydrofabric': name, 'dt_s': dt, 'too_short': int(short.sum()),
            **{position: int(count) for position, count in zip(census.POSITIONS, counts, strict=True)},
            'zero_length': int((short & (length == 0)).sum()),
            'median_length_m': float(np.median(length[short])) if short.any() else np.nan,
        })
    result = pd.DataFrame(rows)
    if not np.all(result[list(census.POSITIONS)].sum(axis=1) == result['too_short']):
        raise ValueError('every short reach must have one position')
    return result


def plot_signs(signs: pd.DataFrame) -> None:
    """Share of reaches too long, inside the window, and too short at each step, one line per hydrofabric."""
    panels = (('too_long_pct', 'c₁ < 0 (too long)'), ('valid_pct', 'All coefficients ≥ 0'),
              ('too_short_pct', 'c₃ < 0 (too short)'))
    figure, grid = plt.subplots(2, 2, figsize=(plotting.WIDTH, 7.0))
    axes, key = grid.flat[:3], grid.flat[3]  # three panels, and the fourth cell holds the legend
    steps = np.arange(len(DT_CENSUS))  # the steps in order, evenly spaced
    for axis, (column, title) in zip(axes, panels, strict=True):
        for name, rows in signs.groupby('hydrofabric', sort=False):
            if not np.array_equal(rows['dt_s'].to_numpy(), DT_CENSUS):
                raise ValueError(f'{name} must have one row per census step, in order')
            # unclipped, so the markers at 0 and 100% show whole
            axis.plot(steps, rows[column], marker=plotting.HYDROFABRIC_MARKERS[name], markersize=5,
                      color=plotting.HYDROFABRIC_COLORS[name], label=name, clip_on=False)
        plotting.step_ticks(axis, DT_CENSUS)
        axis.set_title(title)
        axis.set_ylim(0, 100)
        axis.set_ylabel('Share of reaches (%)')
    key.axis('off')
    key.legend(*axes[0].get_legend_handles_labels(), loc='center')
    plotting.save(figure, 'hydrofabric_census_signs')
    return


def plot_distributions(curves: dict[str, np.ndarray]) -> None:
    """Distribution of travel time by Eq. 17, read also as the length of a second-order reach, with two windows of k."""
    if not curves:
        raise ValueError('no hydrofabrics')
    figure, axis = plt.subplots(figsize=(plotting.WIDTH, 4.0))
    share = np.linspace(0, 100, next(iter(curves.values())).shape[0])  # percent
    marked = np.searchsorted(share, (10, 30, 50, 70, 90))  # where each curve carries its marker
    for name, k in curves.items():
        if k.shape != share.shape:
            raise ValueError('every curve must be sampled at the same shares')
        axis.plot(k, share, color=plotting.HYDROFABRIC_COLORS[name], marker=plotting.HYDROFABRIC_MARKERS[name],
                  markersize=5, markevery=marked.tolist(), label=name)
    for dt in K_WINDOWS:
        low, high = dt / (2 * (1 - census.X)), dt / (2 * census.X)
        axis.axvspan(low, high, color=plotting.GRID, zorder=0, linewidth=0)
        axis.text(np.sqrt(low * high), 75, plotting.step_label(dt), ha='center', va='center', color=plotting.MUTED)
    axis.set_xscale('log')
    axis.set_xlim(10, 1e6)
    axis.set_ylim(0, 100)
    axis.set_xlabel('k by Eq. 17 (s)')
    axis.set_ylabel('Cumulative share of reaches (%)')
    axis.legend(loc='lower right')
    v = float(census.velocity(np.array([2]))[0])  # the celerity of a second-order reach, at which L = v k
    length = axis.secondary_xaxis('top', functions=(lambda k: v * k / 1000, lambda km: 1000 * km / v))
    length.set_xlabel(f'Length of a second-order reach, v = {v:g} m s⁻¹ (km)')
    plotting.save(figure, 'hydrofabric_census_distributions')
    return


if __name__ == '__main__':
    plotting.apply_style()
    config.TABLES.mkdir(parents=True, exist_ok=True)
    drawn = set(FIGURE_HYDROFABRICS)
    if not drawn <= set(HYDROFABRICS) or not drawn <= set(plotting.HYDROFABRIC_COLORS):
        raise ValueError('every hydrofabric of the figures must be censused and have a color')
    results = {'signs': [], 'lengths': [], 'orders': [], 'topology': []}
    curves = {}
    grid = np.linspace(0, 1, 2001)
    for name, load in HYDROFABRICS.items():
        table = load()
        k = census.travel_time(table['length_m'].to_numpy(), table['order'].to_numpy())
        results['signs'].append(sign_table(name, k))
        results['lengths'].append(length_summary(name, table, k))
        results['orders'].append(order_table(name, table, k))
        results['topology'].append(topology_table(name, table, k))
        if name in FIGURE_HYDROFABRICS:
            curves[name] = np.quantile(k, grid)
        print(f'{name}: {len(table):,} reaches')
    frames = {key: pd.concat(parts, ignore_index=True) for key, parts in results.items()}
    for key, frame in frames.items():
        frame.to_csv(config.TABLES / f'hydrofabric_census_{key}.csv', index=False)
        print(f'\n{key}\n{frame.to_string(index=False)}')
    plot_signs(frames['signs'][frames['signs']['hydrofabric'].isin(FIGURE_HYDROFABRICS)])
    plot_distributions(curves)
