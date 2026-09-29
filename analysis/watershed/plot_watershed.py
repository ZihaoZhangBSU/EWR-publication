"""Reproduce F2/S1 from daily model tables and the final native-ASO snapshot CSV.

The earlier notebook in this folder lacks the paper's ASO markers and relative
bias annotations. This driver retains its date-paired statistics and adds the
author's saved final ASO values. It does not regenerate native ASO aggregation.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[2]

def paired_statistics(reference,model,column):
    paired=reference[['date',column]].merge(model[['date',column]],on='date',suffixes=('_I','_N'),validate='one_to_one')
    a=paired[column+'_I'].to_numpy(float);b=paired[column+'_N'].to_numpy(float)
    valid=np.isfinite(a)&np.isfinite(b);a=a[valid];b=b[valid]
    if not len(a):return dict(n=0,pearson_r=np.nan,mean_bias_mm=np.nan,bias_percent=np.nan)
    mb=float(np.mean(b-a));avg=float(np.mean(a))
    r=float(np.corrcoef(a,b)[0,1]) if len(a)>1 and np.std(a)>0 and np.std(b)>0 else np.nan
    return dict(n=len(a),pearson_r=r,mean_bias_mm=mb,bias_percent=100*mb/avg if avg!=0 else np.nan)

def load(path,columns):
    df=pd.read_csv(path,parse_dates=['date'])
    missing=set(columns)-set(df)
    if missing:raise ValueError(f'{path}: missing columns {sorted(missing)}')
    if df.date.isna().any() or df.date.duplicated().any():raise ValueError(f'{path}: missing or duplicate dates')
    wy=df.date.dt.year+(df.date.dt.month>=10).astype(int)
    if not wy.eq(df.water_year).all():raise ValueError(f'{path}: water-year/date conflict')
    return df

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tables',type=Path,default=ROOT/'data/reference/streamflow')
    p.add_argument('--aso-table',type=Path,default=ROOT/'data/reference/watershed/aso_watershed_observations_raw_native.csv')
    p.add_argument('--output',type=Path,default=ROOT/'outputs/watershed')
    p.add_argument('--font',type=Path,default=ROOT/'analysis/wy2026/fonts/NimbusSans-Regular.ttf')
    args=p.parse_args()
    font_manager.fontManager.addfont(str(args.font))
    if font_manager.FontProperties(fname=str(args.font)).get_name()!='Nimbus Sans':raise ValueError('Nimbus Sans required')
    plt.rcParams.update({'font.family':'Nimbus Sans','font.size':18,'axes.titlesize':20,'axes.titleweight':'normal','pdf.fonttype':42,'ps.fonttype':42})
    cols=['date','water_year','swe_mm','snow_depth_mm']
    i=load(args.tables/'isnobal_watershed_daily.csv',cols)
    n=load(args.tables/'noahmp_watershed_daily.csv',cols)
    aso=load(args.aso_table,cols)
    args.output.mkdir(parents=True,exist_ok=True)
    records=[]
    for variable,column,step,stem in [('SWE','swe_mm',100,'Watershed_SWE'),('Snow depth','snow_depth_mm',500,'Watershed_SD')]:
        fig,axs=plt.subplots(3,3,figsize=(18,14),sharey=True)
        maximum=max(float(frame[column].max()) for frame in [i,n,aso])
        ymax=np.ceil(maximum*1.05/step)*step
        for panel,(ax,wy) in enumerate(zip(axs.flat,range(2018,2027))):
            ii=i[i.water_year==wy];nn=n[n.water_year==wy];aa=aso[aso.water_year==wy]
            ax.plot(ii.date,ii[column],color='#1F77B4',lw=2.5,label='iSnobal')
            ax.plot(nn.date,nn[column],color='#D62728',lw=2.5,label='Noah-MP')
            ax.scatter(aa.date,aa[column],color='black',marker='o',s=55,label='ASO',zorder=3)
            stats=paired_statistics(ii,nn,column)
            records.append(dict(water_year=wy,variable=variable,reference='iSnobal',model='Noah-MP',**stats))
            ax.set_title(f'({chr(97+panel)}) WY {wy}',loc='left')
            ax.text(.98,.97,f"r = {stats['pearson_r']:.2f}\nMB = {stats['mean_bias_mm']:.1f} mm ({stats['bias_percent']:.0f}%)",transform=ax.transAxes,ha='right',va='top',fontsize=18)
            ax.set_xlim(pd.Timestamp(wy-1,10,1),pd.Timestamp(wy,9,30));ax.set_ylim(0,ymax)
            ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10,12,2,4,6,8]));ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            ax.spines[['top','right']].set_visible(False)
        axs.flat[0].legend(loc='upper left',frameon=True,facecolor='white',edgecolor='.4',fontsize=18)
        fig.supxlabel('Month');fig.supylabel(variable+' (mm)')
        fig.subplots_adjust(left=.08,right=.99,bottom=.08,top=.96,wspace=.16,hspace=.35)
        fig.savefig(args.output/(stem+'.png'),dpi=600,bbox_inches='tight')
        fig.savefig(args.output/(stem+'.pdf'),bbox_inches='tight');plt.close(fig)
    pd.DataFrame(records).to_csv(args.output/'watershed_yearly_stats_replot_with_aso.csv',index=False)

if __name__=='__main__':main()
