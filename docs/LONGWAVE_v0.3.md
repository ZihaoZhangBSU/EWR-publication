> v0.3.1 coordinate-precision patch: use the v0.3.1 entry point in a complete package update, or replace only the helper in an existing v0.3.0 installation. See ../PATCH_v0.3.1.md. Scientific calculations below are unchanged.

# Daily longwave comparison — version 0.3.0

This replaces the incoming-longwave panel in the revised WY2021/WY2026 comparison with **estimated daily net-longwave difference, Noah-MP minus iSnobal**. No model rerun is needed. The original scripts and version 0.2.0 are preserved.

## Run

Unzip this folder into your existing East River Workflow v0.8.2 repository, for example under `examples/`. Keep the plotting script, both `east_river_*.py` helpers, and `fonts/` together. The existing workflow must be importable (installed or available in a parent `src/` directory). Use its normal scientific Python environment.

To use the existing archived forcing-table directory:

```bash
python East_River_WY2021_WY2026_Comparison_Plot_v0.3.1.py \
  --config east_river_config_output1.yaml \
  --forcing-subdir ldasin_wy2026_forcing_v0.1
```

If you have already run the revised forcing script v0.1.1, omit `--forcing-subdir`; its default is `ldasin_wy2026_forcing_v0.1.1`. This script does not rerun either simulation or require rerunning the forcing script when its existing tables are available.

The included config preserves the source archive's raw-input paths and uses your specified output root:

`/bsuhome/zihaozhang/scratch/iSnobal/East_River_Workflow_Code_v0.8.2/output/output1`

Use your current config if raw input locations have changed. No alternative raw paths are inferred. The script reads `separation_diagnostic_v0.1/tables/separation_watershed_daily_full_period.csv` and the selected forcing directory's `tables/ldasin_watershed_daily_oct1_jun30.csv`.

For each February 1–April 30 day in WY2021 and WY2026, the raw iSnobal `em.nc`, `net_solar.nc`, `snow.nc`, and `awsm_config_backup.ini` must exist for that day **and the following day**. Thus May 1 is also required. Configured native topography, Noah geography, watershed weights/shapefile, and hourly LDASOUT containing `FIRA`, `SNEQV`, and `FSNO` are required. Variable names for SWE and iSnobal solar come from the workflow config; FIRA, FSNO, and net_rad are the verified radiation/snow-fraction variables. Missing inputs raise an explicit error instead of substituting another day.

## Calculation

All longwave means cover midnight to midnight UTC. Positive net longwave means energy toward the surface.

1. For the verified hourly daily-restart configuration, iSnobal daily net all-wave is `(23 * net_rad[D,23:00] + net_rad[D+1,00:00]) / 24`. Today's 00:00 output belongs to the previous day and is not included. The code checks hourly forcing, 24-hour output frequency, restart initialization, daily start/end times, and actual NetCDF timestamps. Runs with assimilation/update sections are rejected for separate interval review.
2. Integrate 25 hourly `net_solar` endpoints from D00 through next D00 with the trapezoidal rule. Subtract that average from daily net all-wave to estimate net longwave.
3. Integrate Noah-MP `-FIRA` over the same 25 endpoints. This is an **hourly-sample approximation** to a daily mean, not an assertion that saved fluxes are exact instantaneous measurements. The inspected public HRLDAS driver writes FIRAXY directly. If your actual writer explicitly stores preceding-hour means instead, select `--noah-lw-time-support hourly-end-mean`, which averages the 24 records ending at D01 through next D00.
4. Screen native iSnobal cells for positive SWE at D00, D23, and next D00, and screen Noah cells for positive SWE and FSNO at every hourly checkpoint. Regrid the screened native radiation to the Noah grid using the workflow's area-average method. Retain only target cells entirely inside the native domain whose contributing native pixels pass the screen (valid coverage at least `1 - 1e-6`). Apply the **same exact fractional watershed weights and common Noah-cell mask** to both models. There is no comparison of a screened iSnobal mean with an unscreened Noah mean.
5. An optional `--minimum-noah-fsno` in [0,1] permits a stricter sensitivity test. Default 0 still requires strictly positive FSNO everywhere in time; it does not assert that Noah cells are fully snow covered. The table includes the selected mean Noah FSNO, area fraction, area in km², and cell count. Days with no common support retain their date and have NaN fluxes (plot gaps), never zero flux.

## Meaning of the gray range

The inspected pySnobal uses linearly interpolated solar forcing evaluated at the start of internal timesteps. Within each hour, its effective mean lies between the left endpoint and the trapezoidal mean. Summing the per-hour lower and upper limits produces a conservative **solar integration range**, conditional on uninterrupted snow. Subtracting solar reverses the bounds. These bounds are propagated through the same regridding and spatial weights as the central estimate and then into the Noah-minus-iSnobal difference.

This shading is **not a confidence interval or total uncertainty**. It excludes unobserved snow appearance/disappearance, Noah-MP temporal sampling error, spatial-remapping error, canopy differences, and model error. Snow at saved checkpoints does not prove snow was present between them. Noah FIRA describes total land/canopy net longwave; iSnobal describes snow-surface radiation. Matching the sample area does not make those surfaces physically identical. Interpret the result as an exploratory snow-conditioned model diagnostic, especially in mixed or forested cells.

Other rows retain their original watershed-wide support. The longwave row alone uses the common snow-present support, which can change each day. Inspect area coverage alongside any apparent longwave differences. The figure's footnote identifies this distinction.

## Outputs

Results go into `output1/wy2021_wy2026_comparison_v0.3/`:

- `figures/figure_wy2021_wy2026_daily_comparison.png` — 600-dpi RGB.
- Matching PDF — vector lines/text with embedded Nimbus Sans.
- `tables/daily_estimated_net_longwave.csv` — each model's estimate, N−I difference, conditional solar bounds, interval, common area and snow-fraction diagnostics.
- `tables/wy2021_wy2026_daily_merged.csv` — all plotted quantities.
- `comparison_manifest.json` — method, selected options, limitations and input provenance. Tables have SHA-256 hashes; large raw files have paths, sizes and modification times to avoid hashing the 147 GB production output.

No old longwave cache is reused. There is one comparison figure, preserving the original 3×2 arrangement, LDASIN temperature, FSNO/SWE and albedo/shortwave rows. Main lines are 2.5 pt, base font 18 pt, WY titles 20 pt, with boxed left-panel legends and shared labels. Neutral gray represents differences; Noah quantities use the established Noah color.

## Validation and limits

`validate_longwave.py` checks the daily weights and signs, 25-endpoint integration, adaptive-step bounds over 200 randomized cases, snow screening, shared area weights, actual NetCDF readers and regridding on synthetic fixtures, missing midnight rejection, incompatible run settings, and PNG/PDF font/export behavior. The uploaded October 2, 2017 sample also passed actual configuration, timestamp, unit and CRS checks.

The layout preview and compact export files in `validation/` are explicitly **synthetic test data**, not WY2021/WY2026 research figures. Full production execution was not possible locally: the complete LDASOUT and daily native archives are on your HPC system. The real October sample alone lacks the next day's files and was not represented as a complete daily model comparison.

## Source references

- [AWSM output accumulator reset and restart handling](https://github.com/iSnobal/awsm/blob/6ebc9482440ffb6f632d6cb653593ceb14d1b5c4/awsm/interface/ipysnobal.py#L358-L435)
- [pySnobal interpolation](https://github.com/iSnobal/pysnobal/blob/849a61796990eeb8dfc774eb5e9ae4d553d50587/pysnobal/c_snobal/libsnobal/do_data_tstep.c#L91-L108) and [averaging before advancing inputs](https://github.com/iSnobal/pysnobal/blob/849a61796990eeb8dfc774eb5e9ae4d553d50587/pysnobal/c_snobal/libsnobal/_do_tstep.c#L114-L165)
- [Snow-free zero flux behavior](https://github.com/iSnobal/pysnobal/blob/849a61796990eeb8dfc774eb5e9ae4d553d50587/pysnobal/c_snobal/libsnobal/_e_bal.c#L24-L62)
- [Public HRLDAS history output of FIRA](https://github.com/NCAR/hrldas/blob/cd96df470220f7d7133cdbccd5f9c5355cf173e2/hrldas/IO_code/module_NoahMP_hrldas_driver.F#L1141-L1177). This source and your sample both use the v20220701 title; that title alone does not prove the exact installed commit. The sampling choice is therefore explicit and recorded.

The installed pySnobal match is based on your confirmation; the AWSM hash is recorded in your supplied configuration.
