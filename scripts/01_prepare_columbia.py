"""
Subset the TDX-Hydro region holding the Columbia to the Columbia basin: its river-route network table, the stream
metadata (length, Strahler order, contributing area, coordinates), and the ERA5 grid weight table.

Run with the project environment:  uv run python scripts/01_prepare_columbia.py
"""


import numpy as np
import pandas as pd
import river_route as rr
import xarray as xr
from river_route.network.streams import subset_network_to_river

from stability import config

REGION_DIR = config.RFS_ROOT / 'routing' / f'region={config.REGION}'
HYDRO_DIR = config.RFS_ROOT / 'hydrography' / f'region={config.REGION}'
METADATA_COLUMNS = [
    'riverId', 'nextRiverId', 'strahlerOrder', 'shreveOrder', 'USContArea', 'DSContArea', 'areaM2', 'Length',
    'musk_k', 'musk_x', 'velocity_factor', 'lat', 'lon',
]


def prepare_network() -> pd.DataFrame:
    """Write the Columbia network table and return it."""
    subset_network_to_river(config.COLUMBIA_OUTLET, REGION_DIR / 'network.parquet', config.NETWORK_FILE)
    network = pd.read_parquet(config.NETWORK_FILE)
    if network['riverId'].iloc[-1] != config.COLUMBIA_OUTLET:
        raise ValueError('the Columbia outlet must be the last row of its DFS-ordered subset')
    if (network['nextRiverId'] == -1).sum() != 1:
        raise ValueError('the Columbia subset must have exactly one outlet')
    rr.Network(config.NETWORK_FILE)  # checks the topology and DFS order
    return network


def prepare_metadata(network: pd.DataFrame) -> pd.DataFrame:
    """Write the stream metadata of the Columbia rivers in network order and return it."""
    metadata = pd.read_parquet(HYDRO_DIR / f'metadata_{config.REGION}.parquet', columns=METADATA_COLUMNS)
    metadata = metadata.set_index('riverId').reindex(network['riverId'].to_numpy()).reset_index()
    if metadata['Length'].isna().any():
        raise ValueError('every Columbia river must have stream metadata')
    if not np.allclose(metadata['musk_k'].to_numpy(), network['muskingumK'].to_numpy()):
        raise ValueError('the network k must be the metadata musk_k')
    metadata['nextRiverId'] = network['nextRiverId'].to_numpy()  # the outlet drains nowhere in the subset
    metadata.to_parquet(config.METADATA_FILE, index=False)
    return metadata


def prepare_weights(network: pd.DataFrame) -> None:
    """Write the ERA5 weight table of the Columbia rivers in network order, with the column names river-route reads."""
    with xr.open_dataset(REGION_DIR / f'gridweights_ERA5_{config.REGION}.nc') as ds:
        weights = ds.to_dataframe().rename(columns={'river_id': 'riverId'})
    order = pd.Series(np.arange(len(network)), index=network['riverId'].to_numpy())
    weights = weights[weights['riverId'].isin(order.index)].copy()
    weights['row'] = order.reindex(weights['riverId'].to_numpy()).to_numpy()
    weights = weights.sort_values(['row'], kind='stable').drop(columns='row').reset_index(drop=True)
    covered = np.isin(network['riverId'].to_numpy(), weights['riverId'].to_numpy())
    if not covered.all():
        raise ValueError(f'{np.count_nonzero(~covered)} Columbia rivers have no grid weights')
    proportion_sums = weights.groupby('riverId')['proportion'].sum()
    if not np.allclose(proportion_sums.to_numpy(), 1, atol=1e-4):
        raise ValueError('the weight proportions of every river must sum to 1')
    weights.index.name = 'index'
    weights.to_xarray().to_netcdf(config.WEIGHTS_FILE)
    return


if __name__ == '__main__':
    config.INPUTS.mkdir(parents=True, exist_ok=True)
    columbia = prepare_network()
    prepare_metadata(columbia)
    prepare_weights(columbia)
    rr.Configs(network_file=config.NETWORK_FILE, grid_weights_file=config.WEIGHTS_FILE, forcing='grid').deep_validate()
    print(f'Columbia basin: {len(columbia):,} rivers written to {config.INPUTS}')
