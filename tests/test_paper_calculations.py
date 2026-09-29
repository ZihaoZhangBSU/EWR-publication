"""Test the actual notebook functions without executing production I/O cells."""
from pathlib import Path
import ast, json
import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
def functions(path,names):
    obj=json.loads((ROOT/path).read_text(encoding='utf-8'))
    source='\n'.join(''.join(c['source']) for c in obj['cells'] if c['cell_type']=='code')
    nodes=[n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert {n.name for n in nodes}==set(names)
    ns={'np':np,'pd':pd}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),path,'exec'),ns)
    return ns

def test_rmb_paired_sums_zero_days_and_missing_values():
    f=functions('analysis/stations/Point_Comp_combined_MB_RMB_f_density.ipynb',['mean_bias_and_rmb'])['mean_bias_and_rmb']
    out=f([0,10,np.nan,20],[5,20,500,np.inf])
    assert out['n']==2
    assert out['mean_bias']==7.5
    assert out['relative_mean_bias_pct']==150
    assert np.isnan(f([0,0],[1,2])['relative_mean_bias_pct'])
    assert f([np.nan],[1])['n']==0
    with pytest.raises(ValueError):f([1],[1,2])

def test_density_units_and_undefined_depth():
    f=functions('analysis/aso/ASO_metrics_with_bulk_density.ipynb',['bulk_density_kg_m3'])['bulk_density_kg_m3']
    np.testing.assert_allclose(f([100,200,0,np.nan],[500,0,100,100]),[200,np.nan,0,np.nan],equal_nan=True)

def test_accumulation_reference_closure():
    d=pd.read_csv(ROOT/'data/reference/accumulation/h1_xu_first_level_metrics.csv')
    gi=d.SWE_Mar15_I_mm-d.SWE_Oct01_I_mm
    gn=d.SWE_Mar15_N_mm-d.SWE_Oct01_N_mm
    si=d.Snowfall_Oct02_Mar15_I_mm
    sn=d.Snowfall_Oct02_Mar15_N_mm
    np.testing.assert_allclose(d.epsilon_total_corrected,gn/gi-1,atol=1e-12,rtol=0)
    np.testing.assert_allclose(d.epsilon_accumulation,sn/si-1,atol=1e-12,rtol=0)
    np.testing.assert_allclose(d.epsilon_ablation,((1-gi/si)-(1-gn/sn))/(gi/si),atol=1e-12,rtol=0)
    np.testing.assert_allclose(1+d.epsilon_total_corrected,(1+d.epsilon_accumulation)*(1+d.epsilon_ablation),atol=1e-12,rtol=0)

def test_watershed_final_figure_statistics():
    import importlib.util
    spec=importlib.util.spec_from_file_location('paper_watershed',ROOT/'analysis/watershed/plot_watershed.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    i=pd.read_csv(ROOT/'data/reference/streamflow/isnobal_watershed_daily.csv',parse_dates=['date'])
    n=pd.read_csv(ROOT/'data/reference/streamflow/noahmp_watershed_daily.csv',parse_dates=['date'])
    expected=pd.read_csv(ROOT/'data/reference/watershed/watershed_yearly_stats_replot_with_aso.csv')
    for row in expected.itertuples():
        col={'SWE':'swe_mm','Snow depth':'snow_depth_mm'}[row.variable]
        actual=mod.paired_statistics(i[i.water_year==row.water_year],n[n.water_year==row.water_year],col)
        for key in ['n','pearson_r','mean_bias_mm','bias_percent']:
            np.testing.assert_allclose(actual[key],getattr(row,key),rtol=0,atol=1e-9)
