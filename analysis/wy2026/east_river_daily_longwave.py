"""Existing-output, approximate daily net longwave; positive toward surface.

Daily UTC support is [D 00:00, D+1 00:00]. Never subtract a 24-point solar
mean from the 23-hour em record. Bounds describe solar quadrature only,
conditional on continuous snow, not total uncertainty or confidence limits.
"""
from __future__ import annotations

import configparser
from pathlib import Path
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer

from east_river_workflow.data_access import decode_noah_times, normalized_units, read_isnobal_grid, read_noah_grid
from east_river_workflow.grids import (load_weights, fractional_polygon_weights,
                                     select_watershed, regrid_area_average)
from east_river_workflow.utils import open_dataset_robust

FLUX_UNITS = {'w/m2', 'w m-2', 'w m^-2', 'watt/m2', 'watt m-2'}
SWE_UNITS = {'kg m-2', 'kg/m2', 'kg m^-2', 'mm'}


def solar_quadrature(endpoints):
    """Return trapezoidal estimate and conservative left-substep enclosure.

    For each hour, effective solar = midpoint - delta/2 * sum(dt_hour**2).
    Since the nonnegative fractional timestep lengths sum to one, the mean
    lies between its left endpoint and midpoint. Bound each hour separately:
    using only the two daily means can hide errors that cancel over a day.
    """
    s = np.asarray(endpoints, dtype=float)
    if s.ndim != 3 or s.shape[0] != 25:
        raise ValueError('Solar requires 25 hourly endpoints on a 2D grid.')
    if np.any(np.isfinite(s) & (s < 0)):
        raise ValueError('Absorbed solar radiation must be nonnegative.')
    midpoint = (s[:-1] + s[1:]) / 2
    return (midpoint.mean(axis=0), np.minimum(s[:-1], midpoint).mean(axis=0),
            np.maximum(s[:-1], midpoint).mean(axis=0))


def estimate_isnobal(rn23, rn_next00, solar25, snow_checkpoints):
    solar, solar_lo, solar_hi = solar_quadrature(solar25)
    snow = np.asarray(snow_checkpoints, dtype=float)
    if snow.shape != (3, *solar.shape):
        raise ValueError('Require snow checkpoints at D00, D23, and next D00.')
    if np.shape(rn23) != solar.shape or np.shape(rn_next00) != solar.shape:
        raise ValueError('Radiation grids differ.')
    good = (np.isfinite(snow).all(axis=0) & (snow > 0).all(axis=0)
            & np.isfinite(solar25).all(axis=0)
            & np.isfinite(rn23) & np.isfinite(rn_next00))
    rn = (23 * np.asarray(rn23) + np.asarray(rn_next00)) / 24
    return tuple(np.where(good, a, np.nan) for a in
                 (rn-solar, rn-solar_hi, rn-solar_lo)), good


def average_noah(fira25, mode='sampled-trapezoid'):
    a = np.asarray(fira25, dtype=float)
    if a.ndim != 3 or a.shape[0] != 25:
        raise ValueError('Require FIRA records at all 25 hourly boundaries.')
    if mode == 'sampled-trapezoid':
        return -((a[:-1]+a[1:])/2).mean(axis=0)
    if mode == 'hourly-end-mean':
        return -a[1:].mean(axis=0)
    raise ValueError(f'Unknown Noah temporal convention: {mode}')


def shared_means(noah, isnobal_fields, weights, mask):
    """Exactly the same Noah cells and fractional watershed weights for all fluxes."""
    w = np.asarray(weights, dtype=float)
    if not np.isfinite(w).all() or (w < 0).any() or (w > 1+1e-6).any() or w.sum() <= 0:
        raise ValueError('Invalid watershed fractional weights.')
    if np.shape(mask) != w.shape:
        raise ValueError('Mask/weights shape mismatch.')
    use = np.asarray(mask, bool) & (w > 0)
    fields = [noah, *isnobal_fields]
    if any(np.shape(v) != w.shape for v in fields):
        raise ValueError('Flux/weights shape mismatch.')
    if any(not np.isfinite(v[use]).all() for v in fields):
        raise ValueError('Nonfinite flux on selected common support.')
    used = w[use].sum()
    means = [float(np.sum(v[use]*w[use])/used) if used else np.nan for v in fields]
    return means, float(used/w.sum()), int(use.sum())


def target_inside_source(source, target):
    """Reject edge targets whose footprint extends outside the native domain.

    Rasterio average normalizes over intersecting pixels; a domain-edge target
    could otherwise appear fully valid despite being only partly covered.
    Check corners and edge midpoints in the native projection.
    """
    rr, cc = np.indices(target.shape)
    convert = Transformer.from_crs(target.crs, source.crs, always_xy=True)
    left, bottom, right, top = source.bounds
    good = np.ones(target.shape, bool)
    for dc, dr in [(0,0), (1,0), (0,1), (1,1), (.5,0), (.5,1), (0,.5), (1,.5)]:
        xx, yy = target.transform * (cc+dc, rr+dr)
        x, y = convert.transform(xx, yy)
        good &= (x >= left-1e-6) & (x <= right+1e-6) & (y >= bottom-1e-6) & (y <= top+1e-6)
    return good


def _day_folder(cfg, date):
    wy = date.year + int(date.month >= 10)
    root = cfg.data['paths']['isnobal_root_template'].format(wy=wy, year=wy, start_year=wy-1)
    return Path(root).expanduser() / cfg.section('isnobal')['daily_folder_pattern'].format(date=date.to_pydatetime(), wy=wy)


def _check_run(folder, date, config_name):
    """Do not infer 23h/1h weights solely from two output timestamps."""
    path = folder / config_name
    if not path.is_file():
        raise FileNotFoundError(f'Missing daily-run configuration: {path}; required to validate averaging intervals.')
    cp = configparser.ConfigParser(interpolation=None)
    cp.read(path, encoding='utf-8')
    if cp.getfloat('time', 'time_step') != 60 or cp.getfloat('awsm system', 'output_frequency') != 24:
        raise ValueError(f'{path}: expected hourly input and 24-hour output frequency.')
    if not cp.get('files', 'init_file').strip():
        raise ValueError(f'{path}: cold start is not supported by the daily-restart calculation.')
    for key, expected in [('start_date', date), ('end_date', date+pd.Timedelta(hours=23))]:
        actual = pd.Timestamp(cp.get('time', key))
        if actual.tzinfo is not None:
            actual = actual.tz_convert('UTC').tz_localize(None)
        if actual != expected:
            raise ValueError(f'{path}: {key} must be {expected}, found {actual}.')
    # Updates reset accumulators in the inspected implementation.
    if any('update' in section.lower() or 'assimilat' in section.lower() for section in cp.sections()):
        raise ValueError(f'{path}: assimilation/update settings require a separate interval audit.')
    return path


def _check_axis_coordinates(stored, reference, axis, path):
    """Allow only exact equality or the file's provable storage rounding.

    AWSM create_netCDF writes x/y as f4, whereas the configured topography
    retains f8 coordinates. Compare against the reference cast to the stored
    floating dtype, rather than using an arbitrary spatial allclose tolerance.
    Always retain the reference grid for georeferencing; no data are shifted.
    """
    actual = np.asarray(stored)
    expected = np.asarray(reference)
    if (actual.ndim != 1 or expected.ndim != 1 or actual.shape != expected.shape
            or actual.size < 2):
        raise ValueError(f'{path}: {axis} coordinate shape {actual.shape} differs from reference {expected.shape}.')
    if not np.isfinite(actual).all() or not np.isfinite(expected).all():
        raise ValueError(f'{path}: nonfinite {axis} coordinate.')
    for values in (actual, expected):
        diffs = np.diff(values.astype(np.float64))
        if not (np.all(diffs > 0) or np.all(diffs < 0)):
            raise ValueError(f'{path}: {axis} coordinates must be strictly monotonic.')
    if np.array_equal(actual, expected):
        return
    if np.issubdtype(actual.dtype, np.floating) and np.array_equal(actual, expected.astype(actual.dtype)):
        return
    max_delta = float(np.max(np.abs(actual.astype(np.float64)-expected.astype(np.float64))))
    raise ValueError(
        f'{path}: {axis} coordinates differ beyond storage rounding from configured topography; '
        f'max absolute offset={max_delta:.9g} coordinate units; stored dtype={actual.dtype}; '
        f'file endpoints=({actual[0]}, {actual[-1]}), '
        f'topography endpoints=({expected[0]}, {expected[-1]}). '
        'Check the configured isnobal_topo_file; do not bypass this mismatch.'
    )


def _native(path, variable, wanted_times, grid, accepted_units, expected_times=None):
    with open_dataset_robust(path, decode_times=True) as ds:
        v = ds[variable]
        if set(v.dims) != {'time', 'y', 'x'}:
            raise ValueError(f'{path}:{variable}: expected time,y,x dimensions, got {v.dims}.')
        if normalized_units(v.attrs.get('units', '')) not in accepted_units:
            raise ValueError(f'{path}:{variable}: unexpected units {v.attrs.get("units")}.')
        for axis in ('x', 'y'):
            _check_axis_coordinates(np.asarray(ds[axis]), getattr(grid, axis), axis, path)
        mapping = v.attrs.get('grid_mapping', 'projection')
        if mapping not in ds:
            raise ValueError(f'{path}: missing projection metadata.')
        attrs = ds[mapping].attrs
        wkt = attrs.get('spatial_ref') or attrs.get('crs_wkt')
        if not wkt or CRS.from_wkt(wkt) != grid.crs:
            raise ValueError(f'{path}: projection differs from configured topography.')
        zone = str(ds.time.attrs.get('time_zone', '')).lower()
        if zone not in {'utc', 'etc/utc'}:
            raise ValueError(f'{path}: explicit UTC time_zone is required, found {zone!r}.')
        times = pd.DatetimeIndex(ds.time.values)
        if times.has_duplicates or not times.is_monotonic_increasing:
            raise ValueError(f'{path}: duplicate or unordered timestamps.')
        if expected_times is not None and not times.equals(pd.DatetimeIndex(expected_times)):
            raise ValueError(f'{path}: unexpected output schedule: {times.tolist()}')
        ix = times.get_indexer(pd.DatetimeIndex(wanted_times))
        if (ix < 0).any():
            raise ValueError(f'{path}: missing requested interval boundary or hourly record.')
        return np.asarray(v.isel(time=ix).transpose('time','y','x').values, dtype=float)


def load_isnobal_day(cfg, date, grid, config_name):
    sec = cfg.section('isnobal')
    tomorrow = date + pd.Timedelta(days=1)
    folders = [_day_folder(cfg, d) for d in (date, tomorrow)]
    inputs = [_check_run(f, d, config_name) for f,d in zip(folders,(date,tomorrow))]
    rn, sw, swe = [], [], []
    for f,d in zip(folders,(date,tomorrow)):
        endpoints = pd.DatetimeIndex([d, d+pd.Timedelta(hours=23)])
        take = endpoints if d == date else endpoints[:1]
        for filename, name, units, sink, wanted, expected in [
            (sec['energy_filename'], 'net_rad', FLUX_UNITS, rn, take, endpoints),
            (sec['snow_filename'], sec['swe_variable'], SWE_UNITS, swe, take, endpoints),
            (sec['absorbed_sw_filename'], sec['absorbed_sw_variable'], FLUX_UNITS, sw,
             pd.date_range(d, periods=24 if d == date else 1, freq='h'), pd.date_range(d,periods=24,freq='h')),
        ]:
            path = f/filename
            sink.append(_native(path, name, wanted, grid, units, expected))
            inputs.append(path)
    fields, good = estimate_isnobal(rn[0][1], rn[1][0], np.concatenate(sw), np.concatenate(swe))
    return fields, good, inputs


def _noah_array(ds, name, indices, units):
    v = ds[name]
    if normalized_units(v.attrs.get('units', '')) not in units:
        raise ValueError(f'Noah {name}: unexpected units {v.attrs.get("units")}.')
    if set(v.dims) != {'Time','south_north','west_east'}:
        raise ValueError(f'Noah {name}: unexpected dimensions {v.dims}.')
    return np.asarray(v.isel(Time=indices).transpose('Time','south_north','west_east').values, dtype=float)


def build_daily_longwave(cfg, dates, noah_mode='sampled-trapezoid',
                        config_name='awsm_config_backup.ini', minimum_noah_fsno=0.0):
    if not 0 <= minimum_noah_fsno <= 1:
        raise ValueError('minimum_noah_fsno must be between 0 and 1.')
    paths = cfg.data['paths']
    noah_grid, *_ = read_noah_grid(paths['noahmp_geo_file'], cfg.section('noahmp'))
    ig, _, topo_valid = read_isnobal_grid(paths['isnobal_topo_file'], cfg.section('isnobal')['dem_variable'])
    weight_path = cfg.output_dir/'cache/east_river_fraction_noahmp.npz'
    inputs = {Path(paths[k]).expanduser() for k in ['noahmp_geo_file','isnobal_topo_file','noahmp_output_file']}
    if weight_path.exists():
        weights = load_weights(weight_path, noah_grid)
        inputs.add(weight_path)
    else:
        sec = cfg.section('watershed')
        watershed = select_watershed(paths['watershed_shapefile'], sec['name_field'], sec['name_contains'])
        weights = fractional_polygon_weights(noah_grid, watershed, chunk_size=int(sec['fractional_chunk_size']))
        shp = Path(paths['watershed_shapefile'])
        inputs.update(shp.parent.glob(shp.stem+'.*'))
    weights = np.asarray(weights, float)
    inside = target_inside_source(ig, noah_grid)
    rows = []
    with open_dataset_robust(paths['noahmp_output_file'], decode_times=False) as ds:
        times = decode_noah_times(ds, cfg.section('noahmp'))
        if times.has_duplicates or not times.is_monotonic_increasing:
            raise ValueError('Noah output times must be unique and increasing.')
        for date in pd.DatetimeIndex(dates):
            print(f'Computing estimated daily net LW: {date.date()}', flush=True)
            expected = pd.date_range(date, periods=25, freq='h')
            ix = times.get_indexer(expected)
            if (ix < 0).any():
                raise ValueError(f'Noah missing hourly/boundary records for {date.date()}.')
            fira = _noah_array(ds, 'FIRA', ix, FLUX_UNITS)
            snow = _noah_array(ds, cfg.section('noahmp')['swe_variable'], ix, SWE_UNITS)
            fsno = _noah_array(ds, 'FSNO', ix, {'-', '1', 'fraction'})
            if fira.shape[1:] != noah_grid.shape:
                raise ValueError('Noah output and geogrid shapes differ.')
            if np.any(np.isfinite(fsno) & ((fsno < 0) | (fsno > 1))):
                raise ValueError('Noah FSNO outside [0,1].')
            n = average_noah(fira, noah_mode)
            ngood = (np.isfinite(fira).all(axis=0) & np.isfinite(snow).all(axis=0)
                     & np.isfinite(fsno).all(axis=0) & (snow > 0).all(axis=0)
                     & (fsno > 0).all(axis=0) & (fsno >= minimum_noah_fsno).all(axis=0))
            native, good, used_paths = load_isnobal_day(cfg, date, ig, config_name)
            inputs.update(used_paths)
            good &= topo_valid
            mapped = []
            coverages = []
            for a in native:
                field, coverage = regrid_area_average(np.where(good,a,np.nan), ig, noah_grid,
                                                     minimum_valid_coverage=1-1e-6)
                mapped.append(field); coverages.append(coverage)
            common = inside & ngood & np.logical_and.reduce([np.isfinite(a) for a in mapped])
            common &= np.logical_and.reduce([a >= 1-1e-6 for a in coverages])
            (mn, mi, lo, hi), fraction, cells = shared_means(n, mapped, weights, common)
            use = common & (weights > 0)
            area = float(weights[use].sum()*noah_grid.cell_area_m2/1e6)
            fsno_mean = float(np.average(fsno.mean(axis=0)[use],weights=weights[use])) if cells else np.nan
            rows.append(dict(date=date, water_year=date.year+int(date.month>=10),
                             net_lw_noah_w_m2=mn, estimated_net_lw_isnobal_w_m2=mi,
                             isnobal_lw_solar_lower_w_m2=lo, isnobal_lw_solar_upper_w_m2=hi,
                             delta_net_lw_n_minus_i_w_m2=mn-mi,
                             delta_lw_solar_lower_w_m2=mn-hi, delta_lw_solar_upper_w_m2=mn-lo,
                             common_area_fraction=fraction, common_area_km2=area,
                             common_noah_cell_count=cells, common_noah_mean_fsno=fsno_mean,
                             interval_start_utc=date, interval_end_utc=date+pd.Timedelta(days=1),
                             duration_hours=24, status='estimated' if cells else 'no_common_snow_support'))
    # Size/mtime for large raw inputs: avoid hashing the 147 GB production file.
    lineage = [{'path':str(p.resolve()),'size_bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}
               for p in sorted(inputs)]
    return pd.DataFrame(rows), lineage
