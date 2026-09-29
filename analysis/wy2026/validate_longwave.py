"""Numerical and file-reader regression tests using explicitly synthetic fixtures.

Run inside the v0.8.2 environment. Optionally provide --sample-dir for the
uploaded 2017 example. No synthetic output is represented as research data.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import sys
import tempfile
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyproj import CRS
from PIL import Image

import east_river_daily_longwave as lw
from east_river_workflow.config import WorkflowConfig, DEFAULTS
from east_river_workflow.grids import grid_from_centers, save_weights
from east_river_publication import configure_matplotlib, save_figure

HERE = Path(__file__).resolve().parent
CHECKS = []


def check(label, condition):
    assert condition, label
    CHECKS.append(label)


def fails(label, func, error=ValueError):
    try:
        func()
    except error:
        CHECKS.append(label)
    else:
        raise AssertionError(label)


def write_native(path, name, a, times, grid, units):
    ds = xr.Dataset({name: (('time','y','x'), np.asarray(a), {'units':units,'grid_mapping':'projection'}),
                     'projection': ((), 0, {'spatial_ref':grid.crs.to_wkt()})},
                    coords={'time':times, 'x':grid.x, 'y':grid.y})
    ds.time.attrs['time_zone'] = 'utc'
    ds.to_netcdf(path)


def setup_fixture(root):
    grid = grid_from_centers('test', np.array([500000.,501000.]), np.array([4301000.,4300000.]), CRS.from_epsg(32613))
    date = pd.Timestamp('2021-02-01')
    data = copy.deepcopy(DEFAULTS)
    data['paths'] = {'isnobal_root_template':str(root/'wy{wy}'),
                     'noahmp_output_file':str(root/'noah.nc'), 'isnobal_topo_file':str(root/'topo.nc'),
                     'noahmp_geo_file':str(root/'geo.nc'), 'output_dir':str(root/'output')}
    data['noahmp'].update(output_start_datetime=str(date), expected_record_count=25,
                         expected_first_time=None, expected_last_time=None)
    cfg = WorkflowConfig(data, root/'fixture.yaml')
    (cfg.output_dir/'cache').mkdir(parents=True)
    weights = np.array([[1., .5],[.25,1.]])
    save_weights(cfg.output_dir/'cache/east_river_fraction_noahmp.npz', weights, grid)
    # Spatial metadata is injected into build_daily_longwave only; all time/data
    # reads, native coordinate checks, raster remapping and aggregation are real.
    (root/'topo.nc').write_text('synthetic spatial-metadata placeholder')
    (root/'geo.nc').write_text('synthetic spatial-metadata placeholder')
    for day in [date, date+pd.Timedelta(days=1)]:
        folder = lw._day_folder(cfg, day);folder.mkdir(parents=True)
        (folder/'awsm_config_backup.ini').write_text(
            f'[time]\ntime_step=60\nstart_date={day}\nend_date={day+pd.Timedelta(hours=23)}\n'
            '[awsm system]\noutput_frequency=24\n[files]\ninit_file=previous/snow.nc\n')
        hours = pd.date_range(day, periods=24, freq='h')
        solar = np.broadcast_to((np.arange(24)+(0 if day==date else 24))[:,None,None],(24,2,2))
        write_native(folder/'net_solar.nc','net_solar',solar,hours,grid,'watt/m2')
        rn = np.full((2,2,2),200.) if day==date else np.full((2,2,2),224.)
        rn[0] = 999. if day==date else 224.  # Today's midnight must be ignored.
        write_native(folder/'em.nc','net_rad',rn,hours[[0,23]],grid,'W m-2')
        snow = np.full((2,2,2),100.)
        if day==date:snow[0,0,0]=0  # Remove a cell with observed new snow.
        write_native(folder/'snow.nc','specific_mass',snow,hours[[0,23]],grid,'kg m-2')
    fira = np.broadcast_to(np.arange(25)[:,None,None], (25,2,2))
    snow = np.full((25,2,2),100.)
    fsno = np.ones((25,2,2));fsno[12,0,1]=0  # Noah mid-day snow absence.
    ds = xr.Dataset({n:(('Time','south_north','west_east'),a,{'units':u}) for n,a,u in
                     [('FIRA',fira,'W/m2'),('SNEQV',snow,'kg m-2'),('FSNO',fsno,'-')]},
                    coords={'Time':np.arange(25)})
    ds.to_netcdf(root/'noah.nc')
    return cfg, grid, date, weights


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sample-dir',type=Path);ap.add_argument('--sample-topo',type=Path);args=ap.parse_args()
    ref=np.array([311004.61127129965,311054.61127129965],dtype=np.float64)
    stored=ref.astype(np.float32)
    check('Regression fixture reproduces old exact-equality failure',not np.array_equal(stored,ref))
    lw._check_axis_coordinates(stored,ref,'x','fixture')
    check('Float32 encoding of float64 coordinates accepted',True)
    lw._check_axis_coordinates(ref,ref,'x','fixture')
    check('Exact float64 coordinates accepted',True)
    fails('Full-cell translation still rejected',lambda:lw._check_axis_coordinates(stored+50,ref,'x','fixture'))
    fails('One float32 ULP beyond encoding rounding rejected',lambda:lw._check_axis_coordinates(np.nextafter(stored,np.float32(np.inf)),ref,'x','fixture'))
    fails('Reversed coordinate order still rejected',lambda:lw._check_axis_coordinates(stored[::-1],ref,'x','fixture'))
    fails('Mismatched coordinate lengths rejected',lambda:lw._check_axis_coordinates(stored[:1],ref,'x','fixture'))
    fails('Nonfinite coordinates rejected',lambda:lw._check_axis_coordinates(np.array([np.nan,2]),ref,'x','fixture'))
    solar = np.broadcast_to(np.arange(25)[:,None,None], (25,2,2))
    q,low,high=lw.solar_quadrature(solar)
    check('Trapezoidal 25-endpoint daily integral', np.all(q==12))
    check('Conservative rising-forcing bounds', np.all(low==11.5) and np.all(high==12))
    fields,good=lw.estimate_isnobal(np.full((2,2),200),np.full((2,2),224),solar,np.ones((3,2,2)))
    check('23:1 midnight weights and matched solar subtraction',np.all(fields[0]==189))
    check('Subtracting solar reverses bounds',np.all(fields[1]<=fields[0]) and np.all(fields[0]<=fields[2]))
    check('Noah sign and sampled integration',np.all(lw.average_noah(solar)==-12))
    check('Noah explicit interval-end averaging option',np.all(lw.average_noah(solar,'hourly-end-mean')==-12.5))
    fails('Reject wrong number of endpoints',lambda:lw.solar_quadrature(solar[:-1]))
    fails('Reject negative absorbed solar',lambda:lw.solar_quadrature(solar-1))
    fails('Reject unknown Noah convention',lambda:lw.average_noah(solar,'guess'))
    snow=np.ones((3,2,2));snow[1,0,0]=0
    f,g=lw.estimate_isnobal(np.ones((2,2)),np.ones((2,2)),solar,snow)
    check('Observed snow loss excluded rather than zero filled',not g[0,0] and np.isnan(f[0][0,0]))
    missing=solar.astype(float).copy();missing[5,1,1]=np.nan
    f,g=lw.estimate_isnobal(np.ones((2,2)),np.ones((2,2)),missing,np.ones((3,2,2)))
    check('Missing hourly solar excludes that native cell',not g[1,1])
    rng=np.random.default_rng(4)
    for _ in range(200):
        v=rng.uniform(0,500,size=(25,2,2));_,l,h=lw.solar_quadrature(v);actual=[]
        for a,b in zip(v[:-1],v[1:]):
            widths=rng.dirichlet(np.ones(int(rng.integers(1,61))))
            t=np.r_[0,np.cumsum(widths)[:-1]]
            actual.append(np.sum((a[None]+(b-a)[None]*t[:,None,None])*widths[:,None,None],axis=0))
        check_value=np.mean(actual,axis=0)
        assert np.all(check_value>=l-1e-9) and np.all(check_value<=h+1e-9)
    check('Solar bounds contain 200 randomized adaptive-step integrations',True)
    means,frac,count=lw.shared_means(np.ones((2,2)),[np.ones((2,2))*2]*3,np.ones((2,2)),np.zeros((2,2),bool))
    check('No common snow produces NaN and zero coverage',np.isnan(means).all() and frac==count==0)
    with tempfile.TemporaryDirectory(prefix='lw_validation_') as temp:
        root=Path(temp);cfg,grid,date,weights=setup_fixture(root)
        native,good,inputs=lw.load_isnobal_day(cfg,date,grid,'awsm_config_backup.ini')
        check('Real NetCDF reader excludes observed appearance',not good[0,0] and good.sum()==3)
        check('Reader ignores current-day midnight and uses next-day midnight',np.all(native[0][good]==189))
        check('Native raster fully inside itself',lw.target_inside_source(grid,grid).all())
        wider=grid_from_centers('wide',grid.x+1000,grid.y,grid.crs)
        check('Partial native-domain targets rejected',not lw.target_inside_source(grid,wider)[:,-1].any())
        with patch.object(lw,'read_noah_grid',return_value=(grid,None,None,None,None)), patch.object(lw,'read_isnobal_grid',return_value=(grid,None,np.ones(grid.shape,bool))):
            table,lineage=lw.build_daily_longwave(cfg,[date])
        r=table.iloc[0]
        check('File-to-table common area and independent masks',r.common_noah_cell_count==2 and np.isclose(r.common_area_fraction,1.25/2.75))
        check('End-to-end longwave sign and delta',r.net_lw_noah_w_m2==-12 and r.estimated_net_lw_isnobal_w_m2==189 and r.delta_net_lw_n_minus_i_w_m2==-201)
        check('Input lineage records actual files',len(lineage)>=10)
        path=lw._day_folder(cfg,date+pd.Timedelta(days=1))/'net_solar.nc'
        with xr.open_dataset(path) as ds:bad=ds.load().isel(time=slice(1,None))
        bad.to_netcdf(path,mode='w')
        fails('Missing next-midnight solar must fail',lambda:lw.load_isnobal_day(cfg,date,grid,'awsm_config_backup.ini'))
        conf=lw._day_folder(cfg,date)/'awsm_config_backup.ini'
        conf.write_text(conf.read_text().replace('output_frequency=24','output_frequency=1'))
        fails('Incompatible output-frequency configuration rejected',lambda:lw._check_run(conf.parent,date,conf.name))
        # Check actual HRLDAS variable dimension/units validation.
        with xr.open_dataset(root/'noah.nc') as ds:
            bad=ds.load();bad.FIRA.attrs['units']='MJ/m2'
            fails('Incorrect radiation units rejected',lambda:lw._noah_array(bad,'FIRA',np.arange(25),lw.FLUX_UNITS))

    if args.sample_dir:
        # Read the real single-day example without inventing the absent next day.
        from netCDF4 import Dataset
        with Dataset(args.sample_dir/'net_solar.nc') as d:
            realgrid=grid_from_centers('real',d['x'][:],d['y'][:],CRS.from_wkt(d['projection'].spatial_ref))
        if args.sample_topo:
            realgrid,_,_=lw.read_isnobal_grid(args.sample_topo,'dem')
            with Dataset(args.sample_dir/'em.nc') as ds:
                check('Real sample differs from topography only by float32 encoding',
                      all(not np.array_equal(ds[a][:],getattr(realgrid,a)) and
                          np.array_equal(ds[a][:],getattr(realgrid,a).astype(ds[a].dtype)) for a in ('x','y')))
        d=pd.Timestamp('2017-10-02')
        lw._check_run(args.sample_dir,d,'awsm_config_backup.ini')
        a=lw._native(args.sample_dir/'em.nc','net_rad',[d,d+pd.Timedelta(hours=23)],realgrid,lw.FLUX_UNITS)
        check('Real archived sample configuration, units, CRS and timestamps',a.shape[0]==2)
        if args.sample_topo:
            check('Actual float32 em.nc reads against float64 supplied topo.nc',np.isfinite(a).any())

    spec=importlib.util.spec_from_file_location('comparison_v03',HERE/'East_River_WY2021_WY2026_Comparison_Plot_v0.3.1.py')
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    check('Revised comparison imports',module.SCRIPT_VERSION=='0.3.1')
    with tempfile.TemporaryDirectory(prefix='lw_pipeline_') as temp:
        root=Path(temp)
        data=copy.deepcopy(DEFAULTS);data['paths']={'output_dir':str(root)}
        cfg=WorkflowConfig(data,root/'fixture.yaml');cfg.source_path.write_text('synthetic config') if hasattr(cfg,'source_path') else None
        config_path=root/'config.yaml';config_path.write_text('synthetic config for mocked load_config')
        dates=pd.DatetimeIndex(['2021-02-01','2021-02-02','2026-02-01','2026-02-02'])
        sep=pd.DataFrame(dict(date=dates,water_year=dates.year,fsno_n_watershed_mean=.5,
                             delta_swe_n_minus_i_mm=1.,delta_effective_albedo_n_minus_i=.01,
                             delta_absorbed_sw_n_minus_i_mj_m2_day=.1))
        forcing=pd.DataFrame(dict(date=dates,water_year=dates.year,precip_mm_day=0.,temperature_mean_c=-5.,hour_count=24))
        for dirname,name,df in [('separation_diagnostic_v0.1','separation_watershed_daily_full_period.csv',sep),
                               ('ldasin_wy2026_forcing_v0.1.1','ldasin_watershed_daily_oct1_jun30.csv',forcing)]:
            target=root/dirname/'tables';target.mkdir(parents=True);df.to_csv(target/name,index=False)
        fake=pd.concat([table]*4,ignore_index=True);fake['date']=dates;fake['water_year']=dates.year
        fake['interval_start_utc']=dates;fake['interval_end_utc']=dates+pd.Timedelta(days=1)
        parsed=module.build_parser().parse_args(['--config',str(config_path),'--focus-end','02-02'])
        def preview_save(fig,cfg,settings,stem):
            path=settings.figures_dir/(stem+'.png');fig.savefig(path,dpi=30);fig.savefig(path.with_suffix('.pdf'));plt.close(fig);return path
        with patch.object(module,'load_config',return_value=cfg), patch.object(module,'build_daily_longwave',return_value=(fake,[])), patch.object(module,'_save_figure',side_effect=preview_save):
            result=module.run(parsed)
        check('Pipeline joins both years and writes tables, plot and manifest',Path(result['figure']).is_file() and len(pd.read_csv(result['merged_table']))==4 and (Path(result['output_root'])/'comparison_manifest.json').is_file())
        check('Manifest labels approximation and includes code fingerprints',len(result['code_sha256'])==3 and 'Estimated' in result['longwave'])
    out=HERE/'validation';out.mkdir(exist_ok=True)
    cfg=WorkflowConfig(copy.deepcopy(DEFAULTS),HERE/'east_river_config_output1.yaml')
    configure_matplotlib(cfg)
    rows=[]
    for year in [2021,2026]:
        dates=pd.date_range(f'{year}-02-01',f'{year}-04-30');t=np.arange(len(dates))
        rows.append(pd.DataFrame(dict(date=dates,water_year=year,fsno_n_watershed_mean=.8-.006*t,
              delta_swe_n_minus_i_mm=10*np.sin(t/8),delta_effective_albedo_n_minus_i=.02*np.cos(t/5),
              delta_net_sw_n_minus_i_w_m2=20*np.sin(t/10),temperature_mean_c=-5+.1*t,
              precip_event=t%12==0,delta_net_lw_n_minus_i_w_m2=5*np.sin(t/9),
              delta_lw_solar_lower_w_m2=5*np.sin(t/9)-2,delta_lw_solar_upper_w_m2=5*np.sin(t/9)+3,
              common_area_fraction=.4+.05*np.sin(t/11))))
    frame=pd.concat(rows,ignore_index=True)
    fig=module.make_figure(frame,cfg,SimpleNamespace(precip_event_threshold_mm_day=1.0))
    fig.text(.5,.985,'SYNTHETIC DATA — LAYOUT TEST ONLY',ha='center',fontsize=16,color='darkred')
    fig.savefig(out/'synthetic_layout_preview.png',dpi=100)
    check('Main comparison renders estimated LW and conditional shading',len(fig.axes)==12)
    plt.close(fig)
    # A compact diagnostic checks the same publication exporter at full 600 dpi.
    fig,ax=plt.subplots(figsize=(5,3));ax.plot([0,1],[0,1],label='Synthetic test',linewidth=2.5);ax.legend(frameon=True)
    png,pdf=save_figure(fig,out/'synthetic_export_check')
    with Image.open(png) as im:check('Export is 600 dpi RGB PNG',im.mode=='RGB' and abs(im.info['dpi'][0]-600)<1)
    from pypdf import PdfReader
    fonts=[]
    for page in PdfReader(pdf).pages:
        for font in page['/Resources']['/Font'].values():
            f=font.get_object();fonts.append(str(f.get('/BaseFont')))
            if '/DescendantFonts' in f:
                desc=f['/DescendantFonts'][0].get_object()['/FontDescriptor'].get_object()
                assert '/FontFile2' in desc
    check('Vector PDF embeds Nimbus Sans',bool(fonts) and all('NimbusSans' in f for f in fonts))
    report={'checks_passed':len(CHECKS),'checks':CHECKS,'test_data':'Synthetic unless explicitly marked real archived sample',
            'not_run':'Full WY2021/WY2026 raw comparison; production LDASOUT and adjacent daily files are on HPC.',
            'fingerprints':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.py')}}
    (out/'validation_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
