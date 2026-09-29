"""East River two-panel study map, revised to plotting protocol v6.

The detailed map uses the supplied cluster script's raw inputs. The regional
locator is rebuilt from bundled Natural Earth data and the supplied watershed.
See README.md for inputs, scientific assumptions, and provenance.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.transform import from_bounds
from rasterio.features import geometry_mask
from shapely.geometry import mapping
from pyproj import Transformer, CRS

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LightSource, ListedColormap
from matplotlib.patches import Patch

from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from PIL import Image
import argparse
import json


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', type=Path, default=Path(
        '/bsuscratch/zihaozhang/iSnobal/post_process/boundary/upload2Chat'))
    p.add_argument('--output-dir', type=Path, default=Path.cwd())
    p.add_argument('--font', type=Path, default=Path(__file__).resolve().parent /
                   'fonts/NimbusSans-Regular.ttf')
    p.add_argument('--basemap-dir', type=Path, default=Path(__file__).resolve().parent / 'basemap')
    p.add_argument('--check-only', action='store_true',
                   help='Validate inputs/font and report DEM metadata; do not render.')
    return p.parse_args()


def configure_font(path):
    if path.is_file():
        fm.fontManager.addfont(str(path))
        name = fm.FontProperties(fname=str(path)).get_name()
        if name != 'Nimbus Sans':
            raise ValueError(f'Expected Nimbus Sans, got {name!r} in {path}')
        resolved = path.resolve()
    else:
        resolved = Path(fm.findfont(fm.FontProperties(family='Nimbus Sans',
                            style='normal', weight='normal'), fallback_to_default=False))
    if fm.FontProperties(fname=str(resolved)).get_name() != 'Nimbus Sans':
        raise ValueError(f'Font verification failed: {resolved}')
    plt.rcParams.update({'font.family': 'Nimbus Sans', 'font.size': 18,
                         'font.weight': 'normal', 'axes.titleweight': 'normal',
                         'axes.labelweight': 'normal', 'pdf.fonttype': 42,
                         'ps.fonttype': 42})
    print(f'Verified Nimbus Sans: {resolved}')
    return str(resolved)


def check_axis(values, name):
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError(f'iSnobal {name} must be a finite 1-D coordinate axis.')
    diff = np.diff(values)
    if not ((diff > 0).all() or (diff < 0).all()):
        raise ValueError(f'iSnobal {name} must be strictly monotonic.')
    step = float(np.median(np.abs(diff)))
    if not np.allclose(np.abs(diff), step, rtol=1e-5, atol=1e-5):
        raise ValueError(f'iSnobal {name} coordinates are not uniformly spaced.')
    return step


def check_isnobal_crs(handle):
    # EPSG:32613 is explicit in the supplied source script, not newly inferred.
    # Check any available WKT/proj4 metadata for conflict with that source.
    found = []
    for obj in [handle] + list(handle.variables.values()):
        for key in ('spatial_ref', 'crs_wkt', 'proj4', 'proj4string'):
            if key in obj.ncattrs():
                val = obj.getncattr(key)
                if isinstance(val, bytes):
                    val = val.decode()
                crs = CRS.from_user_input(str(val))
                if not crs.equals(CRS.from_epsg(32613), ignore_axis_order=True):
                    raise ValueError(f'iSnobal CRS {key} conflicts with source EPSG:32613.')
                found.append(key)
    return ('metadata verified: ' + ', '.join(found)) if found else (
        'EPSG:32613 retained from supplied script; no WKT/proj4 metadata found')


def draw_locator(ax, east, detail_bounds, basemap_dir):
    # Equirectangular overview with latitude-dependent aspect correction.
    # Background data only; no Natural Earth geometry defines the watershed.
    extent = (-123., -103., 33., 48.)
    ax.set_facecolor('#c8dce8')
    for name, color, edge, lw in [
        ('ne_110m_land', '#f0eee3', '0.35', 0.65),
        ('ne_110m_admin_1_states_provinces', 'none', '0.30', 0.7),
        ('ne_110m_lakes', '#c8dce8', '0.35', 0.5),
    ]:
        data = gpd.read_file(basemap_dir / (name + '.shp'))
        if data.crs is None:
            raise ValueError(f'Basemap CRS missing: {name}')
        data.to_crs('EPSG:4326').plot(ax=ax, color=color, edgecolor=edge,
                                    linewidth=lw, zorder=1)
    east_geo = east.to_crs('EPSG:4326')
    east_geo.plot(ax=ax, facecolor='black', edgecolor='black', linewidth=0.7, zorder=5)
    # Transform a densely sampled UTM map frame rather than assuming its edges
    # remain straight in longitude/latitude. This outlines the extent of (b).
    x0, y0, x1, y1 = detail_bounds
    xx = np.r_[np.linspace(x0,x1,40), np.full(40,x1), np.linspace(x1,x0,40), np.full(40,x0)]
    yy = np.r_[np.full(40,y0), np.linspace(y0,y1,40), np.full(40,y1), np.linspace(y1,y0,40)]
    lon, lat = Transformer.from_crs(32613,4326,always_xy=True).transform(xx,yy)
    ax.plot(np.r_[lon,lon[0]], np.r_[lat,lat[0]], color='0.35', lw=1.2, zorder=4)
    center = east.geometry.union_all().centroid
    cl, ct = Transformer.from_crs(32613,4326,always_xy=True).transform(center.x,center.y)
    ax.annotate('East River', xy=(cl,ct), xytext=(-106.2,36.6),
                ha='right', va='center', fontsize=18,
                arrowprops={'arrowstyle':'-', 'color':'black', 'lw':1.0}, zorder=6)
    ax.text(-105.5, 40.5, 'Colorado', ha='center', fontsize=18, color='0.25')
    ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:])
    ax.set_aspect(1 / np.cos(np.deg2rad(40.5)))
    ax.set_anchor('N')
    ax.set_xticks([-120,-115,-110,-105], ['120°W','115°W','110°W','105°W'])
    ax.set_yticks([35,40,45], ['35°N','40°N','45°N'])
    ax.tick_params(length=4, width=1, pad=5, labelsize=18)
    ax.grid(color='0.5', alpha=0.3, linewidth=0.5, linestyle=':')


def run(args):
    # paths
    base = args.base
    noah_file = base/'Noah-MP boundary.nc'
    isnobal_file = base/'iSnobal_boundary.nc'
    watershed_file = base/'East River shape file'/'WBDHU10.shp'
    snotel_file = base/'SNOTEL Location.txt'
    tif_files = sorted((base/'USGS topo').glob('*.tif'))
    out_dir = args.output_dir

    out_png = out_dir / 'east_river_study_domains_two_panel_v0.3.0.png'
    out_pdf = out_dir / 'east_river_study_domains_two_panel_v0.3.0.pdf'

    required = [noah_file, isnobal_file, watershed_file, snotel_file]
    for name in ['ne_110m_land', 'ne_110m_admin_1_states_provinces', 'ne_110m_lakes']:
        required.extend(args.basemap_dir / (name + ext) for ext in ['.shp','.shx','.dbf','.prj'])
    missing = [str(p) for p in required if not p.is_file()]
    if not tif_files:
        missing.append(str(base / 'USGS topo' / '*.tif'))
    if missing:
        raise FileNotFoundError('Required map inputs are missing:\n  ' + '\n  '.join(missing))
    font_path = configure_font(args.font)

    # Noah-MP domain
    with Dataset(noah_file, 'r') as f:
        latc = np.ma.asarray(f['XLAT_C'][0], dtype=float).filled(np.nan)
        lonc = np.ma.asarray(f['XLONG_C'][0], dtype=float).filled(np.nan)
    if latc.ndim != 2 or latc.shape != lonc.shape or min(latc.shape) < 2:
        raise ValueError('XLAT_C/XLONG_C must be matching 2-D corner arrays after selecting time 0.')
    if not np.isfinite(latc).all() or not np.isfinite(lonc).all():
        raise ValueError('Noah-MP domain corners contain missing coordinates.')

    noah_lon = np.concatenate([lonc[0,:], lonc[1:,-1], lonc[-1,-2::-1], lonc[-2:0:-1,0], lonc[:1,:1].ravel()])
    noah_lat = np.concatenate([latc[0,:], latc[1:,-1], latc[-1,-2::-1], latc[-2:0:-1,0], latc[:1,:1].ravel()])

    to_utm = Transformer.from_crs('EPSG:4326','EPSG:32613',always_xy=True)
    noah_x, noah_y = to_utm.transform(noah_lon, noah_lat)

    # iSnobal rectangle
    with Dataset(isnobal_file,'r') as f:
        ix = np.ma.asarray(f['x'][:], dtype=float).filled(np.nan)
        iy = np.ma.asarray(f['y'][:], dtype=float).filled(np.nan)
        isnobal_crs_status = check_isnobal_crs(f)
    dx = check_axis(ix, 'x')
    dy = check_axis(iy, 'y')
    isnobal_xmin = ix.min() - dx/2
    isnobal_xmax = ix.max() + dx/2
    isnobal_ymin = iy.min() - dy/2
    isnobal_ymax = iy.max() + dy/2
    isnobal_boundary_x = np.array([isnobal_xmin,isnobal_xmax,isnobal_xmax,isnobal_xmin,isnobal_xmin])
    isnobal_boundary_y = np.array([isnobal_ymin,isnobal_ymin,isnobal_ymax,isnobal_ymax,isnobal_ymin])

    # watershed
    w = gpd.read_file(watershed_file)
    if w.crs is None or 'name' not in w.columns:
        raise ValueError('Watershed file needs a CRS and the source column name.')
    east = w[w['name'].str.contains('East River', case=False, na=False)].to_crs('EPSG:32613')
    if east.empty or east.geometry.is_empty.any() or not east.geometry.is_valid.all():
        raise ValueError('East River watershed selection is empty or invalid.')
    geom = east.geometry.union_all()

    # stations
    snotel = pd.read_csv(snotel_file, sep=r'\s{2,}|\t+', engine='python')
    required_columns = {'Longitude', 'Latitude', 'Station Name'}
    if not required_columns.issubset(snotel.columns):
        raise ValueError(f'SNOTEL table requires {sorted(required_columns)}; found {list(snotel.columns)}')
    coords = snotel[['Longitude', 'Latitude']].to_numpy(dtype=float)
    if not len(coords) or not np.isfinite(coords).all():
        raise ValueError('SNOTEL coordinates are empty or nonfinite.')
    sx, sy = to_utm.transform(snotel['Longitude'].to_numpy(), snotel['Latitude'].to_numpy())
    # gauge
    gauge_name = 'Almont'
    gx, gy = to_utm.transform(-106.847500, 38.664444)

    # DEM mosaic/reproject
    pad = 7000.0
    xmin = float(np.min(noah_x) - pad)
    xmax = float(np.max(noah_x) + pad)
    ymin = float(np.min(noah_y) - pad)
    ymax = float(np.max(noah_y) + pad)
    bounds_utm = (xmin, ymin, xmax, ymax)
    # Read terrain CRS from the files; do not silently assume NAD83.
    with rasterio.open(tif_files[0]) as first:
        dem_crs = first.crs
        dem_count = first.count
    if dem_crs is None or dem_count != 1:
        raise ValueError('Terrain inputs must be single-band elevation rasters with a CRS.')
    for path in tif_files:
        with rasterio.open(path) as ds:
            if ds.crs != dem_crs or ds.count != 1:
                raise ValueError(f'Inconsistent terrain CRS/band count: {path}')
            if ds.units[0] not in (None, '', 'm', 'meter', 'metre', 'meters', 'metres'):
                raise ValueError(f'Terrain vertical units are not metres: {path}: {ds.units}')
    print(f'DEM: {len(tif_files)} file(s), CRS={dem_crs}; elevation metres retained from source script.')
    print(f'iSnobal CRS: {isnobal_crs_status}')
    if args.check_only:
        print('Input checks passed; no figure written.')
        return
    bounds_geo = transform_bounds('EPSG:32613',dem_crs,*bounds_utm,densify_pts=21)

    srcs = [rasterio.open(p) for p in tif_files]
    try:
        mosaic, mosaic_transform = merge(srcs, bounds=bounds_geo, nodata=-999999.0, dtype='float32')
    finally:
        for ds in srcs:
            ds.close()

    source_dem = mosaic[0]
    res = 30.0
    width = int(np.ceil((xmax-xmin)/res))
    height = int(np.ceil((ymax-ymin)/res))
    transform = from_bounds(xmin, ymin, xmax, ymax, width, height)
    dem = np.full((height,width), np.nan, dtype=np.float32)
    reproject(source=source_dem, destination=dem, src_transform=mosaic_transform, src_crs=dem_crs, src_nodata=-999999.0,
              dst_transform=transform, dst_crs='EPSG:32613', dst_nodata=np.nan, resampling=Resampling.bilinear)

    inside = geometry_mask([mapping(geom)], out_shape=dem.shape, transform=transform, invert=True)
    valid = np.isfinite(dem)
    if not valid.any() or not (inside & valid).any():
        raise ValueError('No valid terrain overlaps the map/watershed.')
    filled = np.where(valid, dem, np.nanmedian(dem[valid]))

    # hillshade + HRUs
    ls = LightSource(azdeg=315, altdeg=45)
    hillshade = ls.hillshade(filled, vert_exag=1.3, dx=abs(transform.a), dy=abs(transform.e))
    hillshade = np.ma.masked_where(~valid, hillshade)

    lower_middle = 2896.0
    middle_upper = 3353.0
    hru = np.full(dem.shape, np.nan)
    hru[inside & valid & (dem < lower_middle)] = 0
    hru[inside & valid & (dem >= lower_middle) & (dem < middle_upper)] = 1
    hru[inside & valid & (dem >= middle_upper)] = 2
    hru_colors = ['#f2e28d', '#b8ddd9', '#d9c7e8']
    cmap = ListedColormap(hru_colors)
    rgba = cmap(np.nan_to_num(hru, nan=0).astype(int))
    rgba[..., 3] = np.where(np.isfinite(hru), 0.58, 0.0)

    xcoords = transform.c + (np.arange(width)+0.5)*transform.a
    ycoords = transform.f + (np.arange(height)+0.5)*transform.e
    masked_dem = np.ma.masked_where((~inside) | (~valid), dem)

    # plot
    fig = plt.figure(figsize=(15, 10))
    locator = fig.add_axes([0.065, 0.405, 0.365, 0.535])
    ax = fig.add_axes([0.49, 0.055, 0.495, 0.885])
    draw_locator(locator, east, bounds_utm, args.basemap_dir)
    extent = (xmin,xmax,ymin,ymax)
    ax.imshow(hillshade, extent=extent, origin='upper', cmap='gray', vmin=0.15, vmax=0.90, interpolation='bilinear', zorder=0)
    ax.imshow(rgba, extent=extent, origin='upper', interpolation='nearest', zorder=1)
    ax.contour(xcoords, ycoords, masked_dem, levels=np.arange(2500,4251,250), colors='black', linewidths=0.35, alpha=0.18, zorder=2)
    east.boundary.plot(ax=ax, color='black', linewidth=1.7, zorder=4)

    noah_color = '#D62728'
    isnobal_color = '#1F77B4'
    ax.plot(noah_x, noah_y, color=noah_color, linewidth=2.5, zorder=6)
    ax.plot(isnobal_boundary_x, isnobal_boundary_y, color=isnobal_color, linewidth=2.5, linestyle=(0, (7, 3)), zorder=5)

    snotel_color = 'black'
    gauge_color = 'black'
    ax.scatter(sx, sy, marker='^', s=120, color=snotel_color, edgecolor='white', linewidth=0.8, zorder=8)
    ax.scatter([gx],[gy], marker='D', s=105, color=gauge_color, edgecolor='white', linewidth=0.8, zorder=8)

    label_offsets = {'Butte':(900,450,'left'),'Schofield Pass':(850,650,'left'),'Upper Taylor':(-850,650,'right')}
    for name,x,y in zip(snotel['Station Name'], sx, sy):
        ox, oy, ha = label_offsets.get(name,(700,500,'left'))
        t = ax.text(x+ox,y+oy,name,ha=ha,va='bottom',fontsize=18,fontstyle='normal',weight='normal',zorder=9)
        t.set_path_effects([pe.withStroke(linewidth=3.2, foreground='white')])
    t = ax.text(gx+850, gy-120, gauge_name, ha='left', va='bottom', fontsize=18, fontstyle='normal', weight='normal', zorder=9)
    t.set_path_effects([pe.withStroke(linewidth=3.2, foreground='white')])

    ax.text((np.min(noah_x)+np.max(noah_x))/2, np.max(noah_y)+1800, 'Noah-MP domain', color=noah_color, ha='center', va='bottom', fontsize=18, weight='normal', bbox=dict(boxstyle='round,pad=0.18', facecolor='white', edgecolor='none', alpha=0.80), zorder=10)
    ax.text((isnobal_xmin+isnobal_xmax)/2, isnobal_ymin+1100, 'iSnobal domain', color=isnobal_color, ha='center', va='bottom', fontsize=18, weight='normal', bbox=dict(boxstyle='round,pad=0.18', facecolor='white', edgecolor='none', alpha=0.80), zorder=10)

    # Boxed legend under the first panel avoids covering either map.
    handles = [
        Patch(facecolor=hru_colors[0], edgecolor='0.25', label='Lower (<2896 m)'),
        Patch(facecolor=hru_colors[1], edgecolor='0.25', label='Middle (2896 to <3353 m)'),
        Patch(facecolor=hru_colors[2], edgecolor='0.25', label='Upper (≥3353 m)'),
        Line2D([], [], color='black', marker='^', linestyle='none', markersize=9, label='SNOTEL'),
        Line2D([], [], color='black', marker='D', linestyle='none', markersize=8, label='USGS gauge'),
        Line2D([], [], color='black', linewidth=1.7, label='East River watershed'),
        Line2D([], [], color=noah_color, linewidth=2.5, label='Noah-MP domain'),
        Line2D([], [], color=isnobal_color, linewidth=2.5, linestyle=(0,(7,3)), label='iSnobal domain'),
    ]
    legend = fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.065, 0.055),
               ncol=1, frameon=True, fancybox=False, facecolor='white',
               edgecolor='0.4', framealpha=0.9, fontsize=18, handlelength=1.8,
               borderpad=0.55, labelspacing=0.35)

    # scale bar
    scale_x = xmin + 4500
    scale_y = ymin + 4500
    scale_length = 10000
    ax.plot([scale_x, scale_x+scale_length], [scale_y, scale_y], color='black', linewidth=4, solid_capstyle='butt', zorder=10)
    for tx in [scale_x, scale_x+scale_length]:
        ax.plot([tx,tx], [scale_y-350, scale_y+350], color='black', linewidth=1.5, zorder=10)
    t = ax.text(scale_x+scale_length/2, scale_y+800, '10 km', ha='center', va='bottom', fontsize=18, weight='normal', zorder=10)
    t.set_path_effects([pe.withStroke(linewidth=3, foreground='white')])

    # Grid-north arrow (EPSG:32613); retain the source's vertical arrow geometry.
    north_x = xmin + 4500
    north_start_y = ymax - 15500
    north_end_y = ymax - 9500
    ax.annotate('', xy=(north_x, north_end_y), xytext=(north_x, north_start_y), arrowprops=dict(facecolor='black', edgecolor='white', linewidth=0.6, width=5, headwidth=16, headlength=18), zorder=10)
    t = ax.text(north_x, north_end_y+900, 'Grid\nN', ha='center', va='bottom', fontsize=18, weight='normal', zorder=10)
    t.set_path_effects([pe.withStroke(linewidth=3, foreground='white')])

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color('black')
    # no title
    ax.set_anchor('N')
    for label, panel in zip(['(a)', '(b)'], [locator, ax]):
        panel.annotate(label, xy=(0,1), xycoords='axes fraction', xytext=(0,8),
                       textcoords='offset points', ha='left', va='bottom',
                       fontsize=18, fontfamily='Nimbus Sans', fontweight='normal',
                       color='black', annotation_clip=False, clip_on=False)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=600, bbox_inches='tight')
    fig.savefig(out_pdf, dpi=600, bbox_inches='tight')
    # Matplotlib writes RGBA by default; final submission PNG is explicitly RGB.
    with Image.open(out_png) as im:
        rgb = im.convert('RGB')
    rgb.save(out_png, dpi=(600, 600))
    rgb.close()
    lineage = {'source_script': 'draw_east_river_domains_cluster.py',
               'protocol': 'East_River_Plotting_Requirements_JHM_AMS_v6.md',
               'font': font_path, 'map_crs': 'EPSG:32613', 'dem_source_crs': str(dem_crs),
               'isnobal_crs': isnobal_crs_status,
               'inputs': [str(p.resolve()) for p in required + tif_files],
               'watershed_selection': 'name contains East River, case insensitive (source rule)',
               'elevation_thresholds_m': [2896, 3353],
               'elevation_units': 'metres, retained from supplied source; not independently established',
               'map_dem_pixel_size_m': [abs(transform.a), abs(transform.e)],
               'elevation_shading': 'display DEM, bilinear resampling; not model-grid HRU membership',
               'gauge_lon_lat': [-106.847500, 38.664444],
               'regional_basemap': 'Natural Earth 1:110m; bundled sources.json contains URLs and checksums',
               'figure_status': 'two-panel map; locator rebuilt, detail redrawn from supplied inputs'}
    (out_dir / 'east_river_study_domains_two_panel_v0.3.0_lineage.json').write_text(
        json.dumps(lineage, indent=2), encoding='utf-8')

    plt.close(fig)
    print(f'Created: {out_png}')
    print(f'Created: {out_pdf}')


if __name__ == '__main__':
    run(parse_args())
