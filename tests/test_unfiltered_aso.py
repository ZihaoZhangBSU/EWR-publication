"""Guard the author's raw-ASO choice and spatial identity checks."""
from pathlib import Path
import sys
import numpy as np
import pytest
import rasterio
import yaml
from affine import Affine
from pyproj import CRS

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/src'))
sys.path.insert(0, str(ROOT / 'analysis/aso'))
from prepare_unfiltered_aso import read_model, read_weights
from east_river_workflow.data_access import read_aso_raster
from east_river_workflow.grids import grid_from_centers, regrid_area_average
from east_river_workflow.config import load_config


def test_raw_extremes_and_negative_values_are_retained(tmp_path):
    path = tmp_path / 'original.tif'
    values = np.array([[6.0, 16.0], [-0.002, -9999.0]], dtype='float32')
    with rasterio.open(path, 'w', driver='GTiff', width=2, height=2, count=1,
                       dtype='float32', crs='EPSG:32613',
                       transform=Affine(50, 0, 1000, 0, -50, 2000), nodata=-9999) as ds:
        ds.write(values, 1)
        ds.update_tags(units='[m]')
    raw, grid, _ = read_aso_raster(path)
    assert raw[0, 0] == 6000
    assert raw[0, 1] == 16000
    assert raw[1, 0] < 0
    assert np.isnan(raw[1, 1])
    result, _ = regrid_area_average(raw, grid, grid, minimum_valid_coverage=.5)
    np.testing.assert_allclose(result, raw, equal_nan=True)


def test_shifted_model_grid_rejected(tmp_path):
    grid = grid_from_centers('test', np.array([25., 75.]), np.array([75., 25.]), CRS.from_epsg(32613))
    path = tmp_path / 'shifted.tif'
    with rasterio.open(path, 'w', driver='GTiff', width=2, height=2, count=1,
                       dtype='float32', crs=grid.crs,
                       transform=grid.transform * Affine.translation(1, 0)) as ds:
        ds.write(np.ones((2, 2), dtype='float32'), 1)
    with pytest.raises(ValueError, match='configured grid'):
        read_model(path, grid)


def test_weight_coordinate_identity_required(tmp_path):
    grid = grid_from_centers('test', np.array([25., 75.]), np.array([75., 25.]), 'EPSG:32613')
    path = tmp_path / 'weights.npz'
    np.savez(path, weights=np.ones((2, 2)), x=grid.x + 50, y=grid.y)
    with pytest.raises(ValueError, match='coordinates'):
        read_weights(path, grid)


def test_publication_configuration_disables_aso_qc(tmp_path):
    # The loader creates output directories. Keep this test off production paths.
    source = ROOT / 'workflow/config/east_river_config.yaml'
    settings = yaml.safe_load(source.read_text(encoding='utf-8'))
    settings['paths']['output_dir'] = str(tmp_path / 'output')
    config = tmp_path / 'publication.yaml'
    config.write_text(yaml.safe_dump(settings), encoding='utf-8')
    cfg = load_config(config)
    assert cfg.section('aso')['quality_control']['enabled'] is False
