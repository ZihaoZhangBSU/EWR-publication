# Scientific audit and release status

## ASO lineage resolved: use unfiltered observations

The supplied archive's configuration enables ASO quality control. Its workflow
reads cleaned ASO products when that setting is enabled. The Desktop
`ASO_metrics_with_bulk_density.ipynb` reads the saved common-grid ASO rasters
and validates against that archive's Table 6. By contrast, the spatial-figure
notebook regrids original ASO rasters, and manuscript v4 states that extremes
were retained.

A new check used all 28 original SWE/depth rasters, the actual workflow's
reader and area-average implementation, the 0.50 coverage setting, and the
provided Noah grid. The Noah grid coordinates were verified against the
archived watershed-weight coordinates. Seven saved ASO rasters differ from
fresh raw-ASO aggregation by more than 0.002 mm. Differences occur for SWE
and/or depth on 7 April 2019, 10 June 2019, 18 May 2022, 1 April 2023, and
23 May 2023. This cannot be explained as ordinary float32 export roundoff.

For 7 April 2019, the manuscript's iSnobal SWE MB/RMSE of −97.9/222.0 mm
matches the saved archive's −97.942296/221.959840 mm. The author's raw-ASO
comparison CSV instead reports −98.513363/223.189877 mm. Merely renaming a
script or disabling a configuration flag cannot reconcile already cached
filtered rasters with a raw-ASO description.

The author selected **unfiltered ASO** on 29 September 2026. The publication
configurations now disable ASO QC. Both metrics notebooks explicitly regenerate
ASO rasters from original observations into `output1/publication_unfiltered_aso`,
recompute differences against unchanged model states, and use fresh validation
tables. They never read the old ASO/difference rasters. Grid identity and watershed
weight coordinates are checked. A successful run records source SHA-256 hashes.
No magnitude clipping is applied; NoData, valid-coverage and metric-specific
depth/support requirements remain part of the original calculation.

The common-grid notebook was fully rerun on all 28 original rasters. Its 56 SWE
date/elevation rows agree with the author's independent raw-ASO CSV
within 1e-5 absolute tolerance across seven numerical columns. Corrected SWE,
depth and density tables are included under `data/reference/aso_unfiltered`.
The native density/aggregation-order notebook also completed successfully,
using a local topography verified against archived model rasters and weights.
The native density/aggregation-order notebook also completed successfully,
using a local topography verified against archived model rasters and weights.
See `MANUSCRIPT_CORRECTIONS.md`; the existing manuscript has not been overwritten.

## Checks completed during repository preparation

- Syntax check of 48 candidate scripts/notebooks passed.
- Representative-year selection reproduces WY2018, WY2022, WY2019 and WY2026.
- All 54 station/model/year SWE bias comparisons match archived Table 4.
- The standalone accumulation calculation reproduces every supplied CSV
  column within floating-point precision.
- Targeted paired-sum RMB tests retain zero-reference days, reject unequal
  array lengths, and handle missing and zero-denominator cases.
- Density-unit tests verify the factor of 1000 and positive-depth requirement.
- 30 workflow and publication tests passed (22 workflow tests, four original
  publication checks, and four unfiltered-ASO/grid-identity checks).
- 37 longwave checks passed, including 23:1 interval weights, sign, integration
  bounds, coordinate encoding, snow screening and PNG/PDF export.
- Full execution of the accumulation, streamflow and postpeak-timing notebooks
  passed. Streamflow metric and plotting-data CSVs and timing CSVs agree with
  the reference outputs within 1e-9 absolute tolerance, with zero relative
  tolerance. Calculations were not changed to obtain agreement.
- The final watershed figure-data CSVs were located alongside the author's
  watershed PDFs. A new explicit plotting driver includes their ASO markers,
  preserves all 18 saved statistical rows within 1e-9, and exports Nimbus Sans
  PNG/PDF figures with panel labels. The older watershed notebook remains a
  source reference, not the final figure driver.

## Remaining limits

The live cluster `output1` and full hourly forcing/model/MODIS archives have
not been accessed. The static source and synthetic tests cannot establish
their final data lineage. Archived caches are regression inputs, not proof
that the current production files are the same. Exact original paper PDF
binary equality is not claimed: PDF metadata and plotting-library versions
can change file bytes without changing plotted values.

The preserved source-reference watershed notebook uses the old workflow's
plotting configuration; use the new final-figure driver instead. No raw
reprocessing should silently adopt archived ASO QC settings while claiming
the manuscript's unfiltered treatment.

## Numerical definitions retained

- Station MB/RMB use model-specific finite paired dates; RMB is the ratio of
  paired sums, not a mean of daily percentages.
- Watershed means use fractional-area weights. Common-grid ASO spatial
  distributions give each retained cell equal weight.
- ASO bulk-density tables compute density after aggregation and use ASO depth
  >100 mm plus positive model depths on three-source support. Native density
  panels and station annual densities have different supports.
- Postpeak scatter retains simultaneous cell-day losses >1 mm in both models.
- Snowpack-release metrics use unsmoothed March–June daily values; only the
  hydrograph is smoothed. SWI/QSNBOT are not routed discharge.
- Forcing climatology shows the previous eight-year mean and full range;
  the separate FSNO-bin diagnostic shows the median and IQR.
- MODIS aggregation is an area fraction of binary snow detections, not
  NDSI divided by 100.
