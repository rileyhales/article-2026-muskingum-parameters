"""
Map of the study area: (a) the Columbia basin's streams with the main stems of the major tributaries and the
observation reaches, and (b) every reach by the sign of its Muskingum coefficients at a routing step of one hour. The
line work is rasterized at the figure resolution, since tens of thousands of vector lines make a PDF of many
megabytes, while the text stays vector.

Run with the river-route environment:  ../river-route/.venv/bin/python scripts/13_study_area_map.py
"""


import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import river_route as rr
from matplotlib import patheffects
from matplotlib.lines import Line2D

from stability import config, hydrofabric, plotting

STREAMS_FILE = config.RFS_ROOT / 'hydrography' / f'region={config.REGION}' / f'streams_{config.REGION}.geo.parquet'
CRS = 'EPSG:5070'  # an equal-area projection that covers the whole basin, the Canadian part included
DT = 3600
CLASS_COLORS = {'too long': '#eb6834', 'valid': '#2a78d6', 'too short': '#4a3aa7'}
STREAM_GRAY = '#b9b8b3'
# label offsets in points from the middle of each main stem, chosen so no two names collide
LABEL_OFFSETS = {'Okanogan': (-72, 4), 'Kettle': (5, -12), 'Kootenay': (6, -2), 'Upper Columbia': (6, 0),
                 'Willamette': (-70, -14), 'Cowlitz': (-60, 6), 'John Day': (4, 8), 'Deschutes': (-40, -18)}


def load_streams(table: pd.DataFrame) -> gpd.GeoDataFrame:
    """The stream lines of the Columbia's reaches in network order, projected, with order and coefficient class."""
    streams = gpd.read_parquet(STREAMS_FILE, columns=['riverId', 'strahlerOrder', 'geometry'])
    streams = streams[streams['riverId'].isin(table['riverId'])].set_index('riverId').loc[table['riverId']]
    if len(streams) != len(table):
        raise ValueError('every Columbia reach must have a stream line')
    too_long, too_short = rr.Network(config.NETWORK_FILE).unstable_mask(float(DT))
    streams['class'] = np.where(too_long, 'too long', np.where(too_short, 'too short', 'valid'))
    streams = streams.reset_index().to_crs(CRS)
    if streams.geometry.is_empty.any():
        raise ValueError('stream lines must not be empty')
    return streams


def draw_streams(axis: plt.Axes, streams: gpd.GeoDataFrame, color) -> None:
    """Streams as lines whose width grows with Strahler order."""
    widths = 0.15 + 0.22 * (streams['strahlerOrder'].to_numpy() - 2)
    streams.plot(ax=axis, color=color, linewidth=widths, rasterized=True)
    axis.set_axis_off()
    axis.set_aspect('equal')
    return


def label_stems(axis: plt.Axes, table: pd.DataFrame, streams: gpd.GeoDataFrame) -> None:
    """Each study main stem drawn over the streams and named at the middle of its length."""
    rows = pd.Series(np.arange(len(table)), index=table['riverId'].to_numpy())
    for name, mouth in hydrofabric.TRIBUTARIES.items():
        stem = hydrofabric.main_stem(table, int(rows[mouth]))
        lines = streams.iloc[stem]
        lines.plot(ax=axis, color=plotting.TREATMENT_COLORS['standard'], linewidth=1.2, rasterized=True)
        middle = lines.geometry.iloc[len(stem) // 2].interpolate(0.5, normalized=True)
        if name != 'Columbia':
            axis.annotate(name, (middle.x, middle.y), color=plotting.INK,
                          xytext=LABEL_OFFSETS.get(name, (4, 3)), textcoords='offset points',
                          path_effects=[patheffects.withStroke(linewidth=2, foreground=plotting.SURFACE)])
    return


if __name__ == '__main__':
    plotting.apply_style()
    columbia = hydrofabric.load()
    lines = load_streams(columbia)
    observed = pd.read_parquet(config.INPUTS / 'observation_rivers.parquet')
    figure, (basin, signs) = plt.subplots(2, 1, figsize=(plotting.WIDTH, 8.5))
    draw_streams(basin, lines, STREAM_GRAY)
    label_stems(basin, columbia, lines)
    points = gpd.GeoSeries(gpd.points_from_xy(observed['lon'], observed['lat']), crs='EPSG:4326').to_crs(CRS)
    basin.scatter(points.x, points.y, s=0.8, color=plotting.INK, zorder=4, rasterized=True)
    basin.set_title('(a) Study main stems and observation reaches')
    basin.legend(handles=[Line2D([], [], color=plotting.TREATMENT_COLORS['standard'], label='Main stem'),
                          Line2D([], [], color=plotting.INK, marker='o', markersize=2, linestyle='',
                                 label='Observation reach')], loc='lower left', bbox_to_anchor=(1.0, 0.0))
    for name in ('too long', 'valid'):
        draw_streams(signs, lines[lines['class'] == name], CLASS_COLORS[name])
    short = lines[lines['class'] == 'too short'].geometry.interpolate(0.5, normalized=True)
    signs.scatter(short.x, short.y, s=0.8, color=CLASS_COLORS['too short'], zorder=4, rasterized=True)
    counts = lines['class'].value_counts()
    signs.legend(handles=[
        Line2D([], [], color=CLASS_COLORS['too long'], label=f'c₁ < 0 ({counts["too long"]:,})'),
        Line2D([], [], color=CLASS_COLORS['valid'], label=f'All ≥ 0 ({counts["valid"]:,})'),
        Line2D([], [], color=CLASS_COLORS['too short'], marker='o', markersize=2, linestyle='',
               label=f'c₃ < 0 ({counts["too short"]:,})'),
    ], loc='lower left', bbox_to_anchor=(1.0, 0.0))
    signs.set_title('(b) Sign of the coefficients at Δt = 1 h')
    plotting.save(figure, 'study_area')
    print(counts.to_string())
