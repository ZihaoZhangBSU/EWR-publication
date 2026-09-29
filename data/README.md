# Reference data

The small CSV snapshots under `reference/` include the author's archived outputs,
final watershed figure data, and newly verified unfiltered-ASO metric tables.
They support numerical regression checks;
they are not raw meteorological or model archives and do not reproduce every
figure on their own. SHA-256 fingerprints and origins are in
`docs/source_manifest.json`.

`reference/aso_unfiltered/` was recomputed from all 14 original ASO flights and
the unchanged archived model states. These tables supersede the earlier
filtered-ASO numerical baseline. The original model/forcing time series are
unaffected by this observation-processing correction.

Full reproduction requires the author's `output1/tables`, `output1/cache`,
`output1/regridded`, native model output, forcing, original ASO rasters, MODIS
granules, and map inputs described in `docs/REPRODUCTION.md`. No model run is
implied by running these plotting scripts.

Observation sources used in the study:

- [MOD10A1 Collection 6.1](https://doi.org/10.5067/MODIS/MOD10A1.061).
- [ASO 50-m snow depth](https://doi.org/10.5067/STOT5I0U1WVI) and
  [50-m SWE](https://doi.org/10.5067/M4TUH28NHL4Z) for 2018–2019.
- [ASO AVA](https://ava.airbornesnowobservatories.com/) for later flights.
- [NRCS SNOTEL](https://wcc.sc.egov.usda.gov/nwcc/): Butte 380,
  Schofield Pass 737, Upper Taylor 1141.
- [USGS East River at Almont, 09112500](https://waterdata.usgs.gov/nwis/uv?site_no=09112500).
- [NOAA HRRR archive](https://registry.opendata.aws/noaa-hrrr-pds/).

Provider terms apply to source observations. No additional license is asserted
for third-party data or bundled fonts; font notices accompany the font files.
