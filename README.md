# East River publication analyses

Code for **Snow-Model Performance across Contrasting Snow Years: A Cross-Community
Comparison of Noah-MP and iSnobal in the East River Watershed**.

The repository is organized around the 12 main-text figures, nine supplementary
figures, and associated numerical tables in manuscript v4. It compares the
author's existing Noah-MP, iSnobal, SNOTEL, ASO, MODIS and USGS products for
WY2018–WY2026. WY2017 is the model spin-up year; WY2026 ends on 30 June.

**ASO treatment:** original, unfiltered observations, as confirmed by the author.
The ASO metrics notebooks rebuild their observation inputs to avoid reusing
filtered archive products. Some earlier manuscript table entries therefore
require correction; see [the correction report](docs/MANUSCRIPT_CORRECTIONS.md).
[The scientific audit](docs/SCIENTIFIC_AUDIT.md) distinguishes verified numerical
regressions from analyses requiring the full production archives.

## Start here

1. [Find the code for each figure and table](docs/FIGURE_INDEX.md).
2. [Install dependencies and reproduce results](docs/REPRODUCTION.md).
3. [Read what has and has not been validated](docs/SCIENTIFIC_AUDIT.md).
4. [Locate input datasets and reference snapshots](data/README.md).

```text
analysis/
  study_area/     F1: model domains and stations
  watershed/     F2, S1: watershed SWE and depth
  stations/      F3, S2, S3, S6 and station tables
  aso/           F4, F5, S7–S9 and spatial/density tables
  accumulation/  F6: October–March snowfall/storage accounting
  postpeak/      F7, F8: conditional daily losses and timing
  streamflow/    F9, S4: snowpack release and observed discharge
  wy2026/        F10–F12, S5: forcing, radiation and MODIS
workflow/        source dependency from the author's v0.8.2 output archive
data/reference/  small archived regression tables
tests/          targeted numerical checks
tools/          notebook execution and source-integrity checks
docs/           figure index, methods, provenance and validation evidence
```

Original research calculations are preserved. The two ASO metrics notebooks now
rebuild unfiltered observation inputs; their metric formulas are unchanged.
Stored notebook outputs were cleared to avoid presenting stale results. Later documented coordinate,
MODIS metadata and v4 figure-label revisions supersede the older comparison ZIP.
[The source manifest](docs/source_manifest.json) records file origins, original
hashes, packaged hashes and any changes. The full manuscript, emails, raw model
archives and personal files are not part of this repository.

No new software license has been chosen for the author's code. Public access
alone does not grant an open-source license. Existing third-party font and map
asset notices remain with those files.
