# Reproduction instructions

## Environment

From the repository root:

```bash
conda env create -f environment.yml
conda activate erw-publication
python tools/check_sources.py
python -m pytest workflow/tests tests -q
```

An existing compatible scientific environment can instead use
`python -m pip install -r requirements.txt`. MODIS HDF4 reading requires a GDAL
build with HDF4 support; a successful Python installation alone does not prove
that this driver is present. Nimbus Sans and its notices are bundled alongside
the map and comparison scripts. The notebook runner registers the included font.

## Inputs and configuration

Keep your existing v0.8.2 workflow directory and its production configuration.
The configurations included here retain the earlier run's paths and settings,
with ASO QC explicitly disabled to implement the author's unfiltered choice.
Adapt filesystem paths to your installation. Existing filtered ASO caches must
not be relabeled as raw observations; use the rebuild described below.

The original research output root is:

```text
/bsuhome/zihaozhang/scratch/iSnobal/East_River_Workflow_Code_v0.8.2/output/output1
```

It contains required daily `tables/`, spatial `cache/`, `regridded/` rasters and
diagnostics. Scripts that compare radiation additionally need native hourly
Noah output and adjacent-day iSnobal `em.nc`, `snow.nc`, `net_solar.nc` and run
configuration backups. Full forcing climatology and MODIS validation require
their respective archives. These large inputs are not committed to GitHub.

## Notebook analyses

F2/S1 can be reproduced directly from the included final figure-data snapshots:

```bash
python analysis/watershed/plot_watershed.py
```

This saves both figures and their statistics under `outputs/watershed/`.
Use `--tables` and `--aso-table` to explicitly select different inputs.
The ASO snapshot table was found with the author's final watershed figures
under `Desktop/Paper structure/Watershed-av`; it is not inferred from a
different gridded ASO summary.

Run from this repository root, using the real existing workflow and output paths:

```bash
WORKFLOW=/bsuhome/zihaozhang/scratch/iSnobal/East_River_Workflow_Code_v0.8.2
OUTPUT="$WORKFLOW/output/output1"

python tools/run_notebook.py analysis/postpeak/SWE_peak_SWE50_revised_v5.ipynb \
  --workflow-root "$WORKFLOW" --output-root "$OUTPUT"
```

Substitute another notebook from `FIGURE_INDEX.md` to reproduce its results.
Run `Pick_snow_year.ipynb` before the station comparisons if you need to verify
the selected case years. Figures and derived tables go to the destinations
defined in each notebook beneath the configured output root. They may replace
same-named derived products; preserve your published baseline before rerunning.

The runner executes the actual notebook cells in order, without rewriting
formulas. It overrides only explicitly named top-level path settings; it never
silently changes sample dates, thresholds, masks or missing-data treatment.
`--set-path NAME PATH` and `--set-json NAME JSON` allow documented notebook
settings to be supplied explicitly. The paths inside loaded YAML files must
also match the intended dataset. The watershed notebook uses the YAML output
directory directly.

For F6, first compute the source CSV:

```bash
python analysis/accumulation/h1_xu_first_level_standalone.py \
  --isnobal "$OUTPUT/tables/isnobal_watershed_daily.csv" \
  --noah "$OUTPUT/tables/noahmp_watershed_daily.csv" \
  --output "$WORKFLOW/h1_xu_first_level_metrics.csv"

python tools/run_notebook.py analysis/accumulation/Dia_Pre_melt.ipynb \
  --workflow-root "$WORKFLOW" --output-root "$OUTPUT" \
  --set-path XU_SOURCE_PATH analysis/accumulation/h1_xu_first_level_standalone.py
```

## Unfiltered ASO metrics and density

The two ASO metrics notebooks automatically call `prepare_unfiltered_aso.py`.
They read original SWE/depth TIFFs from the supplied configuration's
`paths.aso_directory`, validate the grids, reuse existing model states, and
rebuild all ASO and difference rasters. They write separate products under
`$OUTPUT/publication_unfiltered_aso/`, including fresh tables, figures, and a
`lineage.json` containing input hashes. Previously filtered ASO files in
`$OUTPUT/regridded/` are ignored. Rebuilding is intentional on every execution.

```bash
python tools/run_notebook.py analysis/aso/ASO_metrics_with_bulk_density.ipynb \
  --workflow-root "$WORKFLOW" --output-root "$OUTPUT"

python tools/run_notebook.py analysis/aso/ASO_iSnobal_density_grid_scale_diagnostic.ipynb \
  --workflow-root "$WORKFLOW" --output-root "$OUTPUT"
```

The second notebook also requires the configured iSnobal topography and native
model rasters. Both notebooks retain the original formulas and sampling masks;
unfiltered means that original magnitudes are retained, including extremes.
Raster NoData, 50% valid aggregation coverage, and the density metric's depth
screen remain. Reprocessing these observations does not require a model rerun.
Main F4/F5/S7 and S8/S9 already read original ASO directly.

The small tables in `data/reference/aso_unfiltered/` are the corrected baseline.
Use `docs/MANUSCRIPT_CORRECTIONS.md` when transferring values to the paper.

## WY2026 analyses

All commands below use the production configuration, not an invented path to
the input archives. The energy comparison also needs the existing separation
diagnostic tables under `output1/separation_diagnostic_v0.1/`.

```bash
python analysis/wy2026/East_River_LDASIN_WY2026_Forcing_Climatology_v0.1.1.py \
  --config "$WORKFLOW/config/east_river_config.yaml"

python analysis/wy2026/East_River_WY2021_WY2026_Comparison_Plot_v0.3.1.py \
  --config "$WORKFLOW/config/east_river_config.yaml"

python analysis/wy2026/modis_fsno_validation_v0.1.5.py \
  --config "$WORKFLOW/config/east_river_config.yaml" \
  --modis-dir /bsuhome/zihaozhang/scratch/MODIS-SCF
```

If using the earlier compatible forcing tables instead of freshly generated
v0.1.1 tables, add `--forcing-subdir ldasin_wy2026_forcing_v0.1` to the energy
command. Run each script with `--help` for the supported options.

The energy figure's longwave panel is an **estimate**: interval-weighted net
all-wave radiation minus a matching solar estimate, compared with Noah-MP
`-FIRA` on common snow-present support. See `LONGWAVE_v0.3.md`. This is not a
closed whole-watershed energy budget and does not require a model rerun.

## Study map and native ASO supplementary figures

The study-map `--base` is the directory extracted from the author's
`upload2Chat.zip`. It must contain the two model-boundary NetCDFs, SNOTEL
location file, watershed shapefile, and USGS terrain rasters listed in the
script. Use `--check-only` before rendering:

```bash
python analysis/study_area/draw_east_river_domains_cluster_v0.3.0.py \
  --base /bsuscratch/zihaozhang/iSnobal/post_process/boundary/upload2Chat \
  --check-only
```

For S8–S9, explicitly supply the original ASO raster directory and watershed
shapefile; no default input location is assumed:

```bash
python analysis/aso/redraw_aso.py --help
```

Required options are `--rasters` and `--watershed`; `--output` defaults to
`outputs/aso` in this repository. Thresholds flag SWE >5 m **or** depth >15 m;
they do not remove cells. Density distributions use native paired pixels with
positive depth, distinct from the common-grid density metrics.

## Verification scope

`python analysis/wy2026/validate_longwave.py` runs synthetic file-reader,
integration, sign, spatial-support and export checks. It does not rerun the
full production longwave archive. Reports in `docs/validation/` distinguish
these tests from numerical regressions against real archived tables. No
claim of identical results on unavailable production data is made.
