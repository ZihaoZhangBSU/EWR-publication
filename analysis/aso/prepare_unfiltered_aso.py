"""Rebuild ASO comparisons from original rasters without magnitude filtering.

Reuse the existing model states, never the archived ASO or difference rasters.
Write a separate publication_unfiltered_aso directory; leave output1 intact.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/src'))
import numpy as np
import pandas as pd
import rasterio
from east_river_workflow.aso import distribution_statistics, elevation_band_mask, _aso_density_summary
from east_river_workflow.config import load_config
from east_river_workflow.constants import ELEVATION_BANDS
from east_river_workflow.data_access import read_aso_raster, read_noah_grid, read_isnobal_grid
from east_river_workflow.grids import regrid_area_average, write_geotiff


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_model(path, grid):
    with rasterio.open(path) as src:
        if (src.shape != grid.shape or src.crs != grid.crs
                or not src.transform.almost_equals(grid.transform, precision=1e-5)):
            raise ValueError(f'{path}: model raster does not match the configured grid.')
        return src.read(1, masked=True).astype(float).filled(np.nan)


def read_weights(path, grid):
    with np.load(path, allow_pickle=False) as src:
        weights = np.asarray(src['weights'], dtype=float)
        if weights.shape != grid.shape:
            raise ValueError(f'{path}: watershed weights have the wrong shape.')
        for axis in ('x', 'y'):
            if not np.allclose(src[axis], getattr(grid, axis), rtol=0, atol=1e-5):
                raise ValueError(f'{path}: watershed {axis} coordinates do not match.')
        if not np.isfinite(weights).all() or np.any((weights < 0) | (weights > 1)):
            raise ValueError(f'{path}: invalid fractional watershed weights.')
        return weights


def prepare(config, output_root, *, native=False):
    config, output_root = Path(config).resolve(), Path(output_root).resolve()
    cfg = load_config(config)
    if Path(cfg.output_dir).resolve() != output_root:
        raise ValueError('Configuration output_dir must match the supplied output root.')
    dest = output_root / 'publication_unfiltered_aso'
    (dest / 'tables').mkdir(parents=True, exist_ok=True)
    # Invalidate a previous successful-run marker before touching products.
    marker = dest / 'lineage.json'
    marker.unlink(missing_ok=True)
    grid, dem, *_ = read_noah_grid(cfg.data['paths']['noahmp_geo_file'], cfg.section('noahmp'))
    weight_path = output_root / 'cache/east_river_fraction_noahmp.npz'
    weights = read_weights(weight_path, grid)
    lineage = {'aso_treatment': 'original observations; no magnitude clipping or outlier removal',
               'minimum_valid_coverage': float(cfg.section('aso')['minimum_valid_coverage']),
               'model_states': 'existing archived states; unchanged time selection',
               'native_products': native, 'inputs': []}
    def record(path):
        lineage['inputs'].append({'path': str(path), 'sha256': sha256(path)})
    for path in (config, weight_path, Path(cfg.data['paths']['noahmp_geo_file'])):
        record(path)
    if native:
        native_path = Path(cfg.data['paths']['isnobal_topo_file'])
        ig, _, _ = read_isnobal_grid(native_path, cfg.section('isnobal')['dem_variable'])
        iw_path = output_root / 'cache/east_river_fraction_isnobal.npz'
        iw = read_weights(iw_path, ig)
        record(native_path)
        record(iw_path)
    rows = []
    for date in cfg.section('aso')['dates']:
        date = str(date)
        observations = {}
        for variable, slug, pattern in (
            ('SWE', 'swe', 'swe_pattern'), ('Snow depth', 'snow_depth', 'snow_depth_pattern')
        ):
            raw_path = Path(cfg.data['paths']['aso_directory']) / cfg.section('aso')[pattern].format(date=date)
            raw, raw_grid, _ = read_aso_raster(raw_path)
            record(raw_path)
            obs, coverage = regrid_area_average(raw, raw_grid, grid,
                minimum_valid_coverage=lineage['minimum_valid_coverage'])
            obs = np.where(weights > 0, obs, np.nan)
            observations[variable] = obs
            prefix = f'{date}_{slug}_'
            write_geotiff(dest / 'regridded' / (prefix + 'aso_on_noah.tif'), obs, grid)
            write_geotiff(dest / 'regridded' / (prefix + 'aso_coverage_noah.tif'), coverage, grid)
            for key, difference, product in (
                ('isnobal_on_noah', 'difference_isnobal_noah', 'iSnobal-ASO on Noah-MP grid'),
                ('noah_native', 'difference_noah', 'Noah-MP-ASO on Noah-MP grid')
            ):
                model_path = output_root / 'regridded' / (prefix + key + '.tif')
                model = read_model(model_path, grid)
                record(model_path)
                shutil.copy2(model_path, dest / 'regridded' / model_path.name)
                diff = np.where(weights > 0, model - obs, np.nan)
                write_geotiff(dest / 'regridded' / (prefix + difference + '.tif'), diff, grid)
                for band in ELEVATION_BANDS:
                    mask = (weights > 0) & elevation_band_mask(dem, band)
                    rows.append({'date': pd.Timestamp(date), 'variable': variable,
                        'product': product, 'elevation_band': band,
                        **distribution_statistics(diff[mask])})
            if native:
                obs_i, cov_i = regrid_area_average(raw, raw_grid, ig,
                    minimum_valid_coverage=lineage['minimum_valid_coverage'])
                obs_i = np.where(iw > 0, obs_i, np.nan)
                model_path = output_root / 'regridded' / (prefix + 'isnobal_native.tif')
                model_i = read_model(model_path, ig)
                record(model_path)
                shutil.copy2(model_path, dest / 'regridded' / model_path.name)
                for key, values in [('aso_on_isnobal', obs_i), ('aso_coverage_isnobal', cov_i),
                                    ('difference_isnobal_native', model_i - obs_i)]:
                    write_geotiff(dest / 'regridded' / (prefix + key + '.tif'), values, ig)
        for band in ELEVATION_BANDS:
            density = _aso_density_summary(observations['SWE'], observations['Snow depth'],
                (weights > 0) & elevation_band_mask(dem, band))
            for row in rows:
                if row['date'] == pd.Timestamp(date) and row['elevation_band'] == band:
                    row.update(density)
        print(f'Rebuilt unfiltered ASO: {date}', flush=True)
    pd.DataFrame(rows).to_csv(dest / 'tables/table6_aso_validation.csv', index=False)
    marker.write_text(json.dumps(lineage, indent=2), encoding='utf-8')
    return dest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    parser.add_argument('--native', action='store_true', help='Also rebuild native-grid ASO for the density-scale diagnostic.')
    args = parser.parse_args()
    print(prepare(args.config, args.output_root, native=args.native))
