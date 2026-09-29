# Manuscript corrections for unfiltered ASO

The author selected original, unfiltered ASO observations on 29 September 2026.
The existing manuscript is preserved. This report identifies updates needed
because earlier metrics tables used filtered archived rasters while the
spatial-figure code already used original observations.

## Tables to replace

- Main Table 4 (SWE) and supplementary Table S3 (depth): use the `All` rows
  of `data/reference/aso_unfiltered/aso_spatial_metrics_by_elevation.csv`.
  Use `aso_reference_mean_isnobal_support_mm` for the equally weighted
  ASO table mean; it matches the Noah-MP support mean for these flights.
  Do not substitute the fractional-area-weighted watershed mean.
- Supplementary Table S4 (density): use `aso_bulk_density_metrics_all.csv`
  in the same directory. The original depth screen and common support remain.
- Refresh bold/underlined extrema after replacing the table values.
  Check narrative numbers quoted from these tables against the corrected CSVs.

The machine-readable `validation/aso_all_metric_changes.csv` includes changes
to reference means, correlations, relative errors, density and sample counts.
Some changes disappear at manuscript display precision. The table below
lists the 23 whole-domain MB/RMSE entries that change at one decimal place.

| Date | Variable | Model | Metric | Archived (mm) | Unfiltered (mm) |
|---|---|---|---|---:|---:|
| 2019-04-07 | SWE | iSnobal | MB | -97.9 | -98.5 |
| 2019-04-07 | SWE | iSnobal | RMSE | 222.0 | 223.2 |
| 2019-04-07 | SWE | Noah-MP | MB | -129.8 | -130.3 |
| 2019-04-07 | SWE | Noah-MP | RMSE | 237.1 | 238.3 |
| 2019-04-07 | Snow depth | iSnobal | MB | -141.9 | -142.5 |
| 2019-04-07 | Snow depth | iSnobal | RMSE | 378.1 | 379.6 |
| 2019-04-07 | Snow depth | Noah-MP | MB | -268.7 | -269.2 |
| 2019-04-07 | Snow depth | Noah-MP | RMSE | 463.2 | 464.5 |
| 2019-06-10 | SWE | iSnobal | MB | -56.5 | -56.8 |
| 2019-06-10 | SWE | iSnobal | RMSE | 182.0 | 182.8 |
| 2019-06-10 | SWE | Noah-MP | MB | -58.2 | -58.5 |
| 2019-06-10 | SWE | Noah-MP | RMSE | 183.6 | 184.1 |
| 2019-06-10 | Snow depth | iSnobal | MB | -188.5 | -188.7 |
| 2019-06-10 | Snow depth | iSnobal | RMSE | 430.3 | 430.4 |
| 2019-06-10 | Snow depth | Noah-MP | MB | -123.7 | -123.8 |
| 2022-05-18 | SWE | Noah-MP | MB | -4.3 | -4.4 |
| 2022-05-18 | SWE | Noah-MP | RMSE | 125.6 | 125.5 |
| 2023-04-01 | SWE | iSnobal | RMSE | 216.8 | 216.9 |
| 2023-04-01 | SWE | Noah-MP | RMSE | 222.9 | 223.0 |
| 2023-05-23 | SWE | iSnobal | MB | -28.5 | -28.6 |
| 2023-05-23 | SWE | iSnobal | RMSE | 151.6 | 151.8 |
| 2023-05-23 | SWE | Noah-MP | MB | -14.4 | -14.5 |
| 2023-05-23 | SWE | Noah-MP | RMSE | 168.2 | 168.3 |

## Scope of the correction

Seven of the 28 common-grid ASO rasters differ beyond export roundoff:
SWE/depth on 7 April and 10 June 2019, and SWE on 18 May 2022,
1 April 2023 and 23 May 2023. No common-grid support cells changed in
the raw-versus-archived raster audit. Density statistics can change
when either constituent changes; their corrected tables are included.

The WY2026 ASO rasters and the 24 May 2018 density example are unchanged.
The native extreme-cell map retains its SWE >5 m OR depth >15 m
flagging criteria; those flags do not remove observations.
The model time series, station, forcing and streamflow calculations were
not changed by this ASO correction.

## Verification

Both ASO metrics notebooks completed with original observations.
All 56 SWE date/elevation rows and seven numerical columns agree with
the author's independent unfiltered-ASO CSV within 1e-5 absolute tolerance.
The largest difference is below 5e-7, consistent with the CSV rounding.
Input fingerprints are recorded in `validation/aso_input_fingerprints.json`.
This verifies these archived model-state comparisons; it does not prove
that inaccessible current cluster files are identical.
