from pathlib import Path
import argparse
import json, hashlib, csv
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description='Regenerate the raw ASO supplementary diagnostics without range-based filtering.')
parser.add_argument('--rasters',type=Path,required=True,help='Folder containing original paired ASO SWE and SD rasters.')
parser.add_argument('--watershed',type=Path,required=True,help='WBDHU10 shapefile containing the East River polygon.')
parser.add_argument('--font',type=Path,default=ROOT/'analysis/wy2026/fonts/NimbusSans-Regular.ttf')
parser.add_argument('--output',type=Path,default=ROOT/'outputs/aso')
args=parser.parse_args()
OUT=args.output
FIG=OUT/'Figures/Supplement'
FIG.mkdir(parents=True,exist_ok=True)
font=args.font
fm.fontManager.addfont(str(font))
assert fm.FontProperties(fname=str(font)).get_name()=='Nimbus Sans'
plt.rcParams.update({'font.family':'Nimbus Sans','font.size':18,'axes.titlesize':20,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
w=gpd.read_file(args.watershed)
col=next(c for c in w if c.lower()=='name')
w=w[w[col].str.lower()=='east river'].to_crs(32613)
assert len(w)==1
raw=args.rasters
rows=[]; points=[]; densities=[]; lineage=[]
for sf in sorted(raw.glob('ASO_50M_SWE_*.tif')):
    date=sf.stem[-8:]; df=sf.with_name(sf.name.replace('_SWE_','_SD_'))
    with rasterio.open(sf) as s, rasterio.open(df) as d:
        assert s.crs==d.crs and s.transform==d.transform and s.shape==d.shape
        assert s.crs.to_epsg()==32613
        swe=s.read(1,masked=True).filled(np.nan).astype(float)
        depth=d.read(1,masked=True).filled(np.nan).astype(float)
        inside=geometry_mask(w.geometry, out_shape=s.shape,transform=s.transform,invert=True,all_touched=False)
        paired=inside & np.isfinite(swe) & np.isfinite(depth)
        extreme=paired & ((swe>5)|(depth>15))
        rr,cc=np.where(extreme); x,y=rasterio.transform.xy(s.transform,rr,cc)
        for r,c,xx,yy in zip(rr,cc,x,y):
            points.append({'date':date,'x_m':xx,'y_m':yy,'swe_m':swe[r,c],'depth_m':depth[r,c]})
        # Density requires positive depth; negative SWE is not silently clipped.
        sel=paired & (depth>0)
        rho=1000*swe[sel]/depth[sel]
        q=np.percentile(rho,[5,25,50,75,95])
        densities.append({'date':pd.to_datetime(date),'q05':q[0],'q25':q[1],'median':q[2],'q75':q[3],'q95':q[4],'n':len(rho)})
        rows.append({'date':date,'n_paired':int(paired.sum()),'n_extreme':int(extreme.sum()),'n_swe_gt5':int((paired&(swe>5)).sum()),'n_depth_gt15':int((paired&(depth>15)).sum()),'max_swe_m':float(swe[paired].max()),'max_depth_m':float(depth[paired].max()),'extreme_area_km2':float(extreme.sum()*abs(s.transform.a*s.transform.e)/1e6)})
    for p in (sf,df): lineage.append({'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pts=pd.DataFrame(points); stats=pd.DataFrame(rows); den=pd.DataFrame(densities)
if len(pts):
    tree=cKDTree(pts[['x_m','y_m']].values)
    pts['recurs_within125m']= [any(pts.iloc[j]['date']!=pts.iloc[i]['date'] for j in tree.query_ball_point(xy,125)) for i,xy in enumerate(pts[['x_m','y_m']].values)]
stats.to_csv(OUT/'ASO_extreme_counts_v1.csv',index=False)
pts.to_csv(OUT/'ASO_extreme_locations_v1.csv',index=False)
den.to_csv(OUT/'ASO_density_distributions_v1.csv',index=False)
(OUT/'ASO_redraw_lineage_v1.json').write_text(json.dumps({'inputs':lineage,'threshold':'SWE > 5 m OR depth > 15 m','support':'paired finite SWE and depth; source NoData excluded; watershed pixel centers; native 50 m grids','recurrence':'at least one cell from a DIFFERENT flight within 125 m (Euclidean); not proof of persistence or accuracy','density':'1000*SWE/depth; paired finite, depth > 0; no clipping; percentiles are spatial variability, not confidence intervals','watershed_area_km2':float(w.area.sum()/1e6)},indent=2))

def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(FIG/(name+'.png'),dpi=600,bbox_inches='tight',facecolor='white')
    (OUT/'previews').mkdir(exist_ok=True)
    fig.savefig(OUT/'previews'/f'{name}_preview.png',dpi=120,bbox_inches='tight')
    plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(16,9),gridspec_kw={'width_ratios':[1.1,1]},layout='constrained')
ax=axs[0]; w.plot(ax=ax,facecolor='#f6f6f6',edgecolor='0.25',linewidth=1.2)
colors=plt.get_cmap('tab10').colors
for i,(date,g) in enumerate(pts.groupby('date')):
    ax.scatter(g.x_m,g.y_m,s=35,color=colors[i],label=pd.to_datetime(date).strftime('%Y-%m-%d'),zorder=3)
    recur=g[g.recurs_within125m]
    ax.scatter(recur.x_m,recur.y_m,s=95,facecolors='none',edgecolors='black',linewidths=.65,zorder=4)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}'))
ax.yaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}'))
ax.set_xlabel('UTM easting (km)');ax.set_ylabel('UTM northing (km)')
ax.set_title('(a) Threshold-exceeding cells',loc='left')
handles,labels=ax.get_legend_handles_labels()
handles.append(Line2D([],[],marker='o',markersize=9,markerfacecolor='none',markeredgecolor='black',linestyle='none'))
labels.append('Recurrence within 125 m')
ax.legend(handles,labels,loc='lower left',fontsize=14,frameon=True,facecolor='white',framealpha=1,edgecolor='0.3')
ax.annotate('N',xy=(.91,.95),xytext=(.91,.83),xycoords='axes fraction',ha='center',arrowprops={'arrowstyle':'->','lw':1.5})
ax=axs[1]; ypos=np.arange(len(stats))
ax.scatter(stats.n_extreme,ypos,c='black',s=55,zorder=3)
ax.hlines(ypos,0,stats.n_extreme,colors='.7',lw=1.2)
for y,n in zip(ypos,stats.n_extreme):ax.annotate(str(n),(n,y),xytext=(7,0),textcoords='offset points',va='center',fontsize=15)
ax.set_yticks(ypos,[pd.to_datetime(x).strftime('%Y-%m-%d') for x in stats.date]);ax.invert_yaxis()
ax.set_xlim(-.8,max(stats.n_extreme)*1.25+1)
ax.set_xlabel('Number of paired 50 m cells');ax.set_title('(b) Counts by ASO flight',loc='left')
ax.grid(axis='x',alpha=.2);save(fig,'ASO_extremes_revised_v1')

fig,axs=plt.subplots(1,2,figsize=(17,8),layout='constrained')
ax=axs[0]
for i,row in den.iterrows():
    ax.plot([i,i],[row.q05,row.q95],color='.55',lw=1)
    ax.plot([i,i],[row.q25,row.q75],color='.15',lw=7,solid_capstyle='butt')
    ax.scatter(i,row['median'],c='white',edgecolor='black',s=35,zorder=3)
ax.set_xticks(range(len(den)),den.date.dt.strftime('%Y-%m-%d'),rotation=65,ha='right')
ax.set_ylabel('ASO-implied bulk density (kg m$^{-3}$)');ax.set_title('(a) Spatial distributions',loc='left')
ax=axs[1]
for i,(year,g) in enumerate(den.groupby(den.date.dt.year)):
    t=(g.date-g.date.min()).dt.days
    ax.plot(t,g['median'],marker='o',lw=2.5,color=colors[i],label=f'WY{year}')
ax.set_xlabel('Days since first flight of the water year');ax.set_ylabel('Median density (kg m$^{-3}$)')
ax.set_title('(b) Within-year medians',loc='left');ax.legend(ncol=2,fontsize=15,frameon=True,edgecolor='.3')
for ax in axs:ax.grid(axis='y',alpha=.2)
save(fig,'ASO_density_revised_v1')
print(stats.to_string(index=False));print('Density quantiles:');print(den.to_string(index=False))
