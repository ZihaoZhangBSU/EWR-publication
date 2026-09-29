#!/usr/bin/env python3
"""
East River MODIS–Noah-MP snow-cover validation companion workflow.

Version
-------
0.1.5

Purpose
-------
Validate native Noah-MP FSNO against MODIS Terra MOD10A1 V061 snow-cover
observations over the East River watershed.

Primary scientific definition
-----------------------------
1. Read MOD10A1 NDSI_Snow_Cover.
2. Retain values 0..100 with Basic_QA in {0, 1}.
3. Convert valid MODIS pixels to binary snow presence:
       NDSI_Snow_Cover > 0 -> 1
       NDSI_Snow_Cover = 0 -> 0
       invalid/cloud/QA-filtered -> NaN
4. Area-average the binary field from native MODIS (~463 m) to the native
   Noah-MP grid. This produces MODIS-derived fSCA [0, 1].
5. Compute MODIS valid source-area coverage separately and require the
   configured per-Noah-cell threshold (default 0.75).
6. Recover MODIS source-swath times from:
       granule_pnt
       -> GRANULEPOINTERARRAY index
       -> GRANULEBEGINNINGDATETIMEARRAY / GRANULEENDINGDATETIMEARRAY
7. Match each source swath to the nearest hourly Noah-MP UTC state.
   The representative swath time is the midpoint between beginning and ending
   times. This choice is explicit in the manifest and diagnostics.
8. If a Noah cell contains valid MODIS source area from more than one granule,
   combine the corresponding hourly FSNO states using the MODIS valid-area
   fractions from those granules.
9. Compare Noah FSNO and MODIS-derived fSCA on exactly the same paired East
   River support.

Bias sign convention
--------------------
    Noah-MP - MODIS

Known project paths used as defaults
------------------------------------
Config:
    east_river_config_output1.yaml beside this script (override with --config).

MODIS:
    /bsuhome/zihaozhang/scratch/MODIS-SCF

The Noah LDASOUT, Noah geography file, watershed shapefile, and main workflow
output directory are read from the v0.8.2 YAML config rather than duplicated
here.

Outputs
-------
<main workflow output>/modis_fsno_validation_v0.1/
    cache/
        modis_fsno_wyYYYY.npz
    tables/
        modis_fsno_watershed_daily.csv
        modis_fsno_cell_statistics.csv
        modis_fsno_annual_summary.csv
    diagnostics/
        modis_file_inventory.csv
        modis_date_alignment.csv
        modis_qa_summary.csv
        modis_validation_manifest.json
    figures/
        figure01_modis_fsno_yearly_comparison.<format>
        figure02_modis_fsno_yearly_difference.<format>
        figure03_modis_fsno_scatter.<format>
        figure04_modis_fsno_spatial_bias.<format>

Notes
-----
- Primary basin comparison requires at least 75% paired East River watershed area.
- The per-Noah-cell MODIS valid-source threshold is 0.75; both primary coverage thresholds are 0.75 in v0.1.3.
- Missing MODIS dates remain missing. They are not interpolated.
- Files outside the configured analysis window are inventoried but not used.
- Duplicate MODIS dates are treated as an error.
- Algorithm_Flags_QA is retained for diagnostics but is not an additional
  primary mask in v0.1.
- The script fails loudly if granule_pnt cannot be mapped uniquely to MODIS
  acquisition metadata. It does not guess.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from pyproj import CRS
from rasterio.windows import Window, from_bounds
from rasterio.warp import Resampling, reproject, transform_bounds


# =============================================================================
# Bootstrap East River imports
# =============================================================================
def _bootstrap_east_river_import() -> None:
    """Find the v0.8.2 source tree when the package is not installed."""
    try:
        import east_river_workflow  # noqa: F401
        return
    except ImportError:
        pass

    here = Path(__file__).resolve()
    candidates: list[Path] = []

    # Exact project location from the supplied East River v0.8.2 source context.
    candidates.append(
        Path(
            "/bsuscratch/zihaozhang/iSnobal/"
            "East_River_Workflow_Code_v0.8.2/src"
        )
    )

    for parent in [here.parent, *here.parents]:
        candidates.extend(
            [
                parent / "src",
                parent / "East_River_Workflow_Code_v0.8.2" / "src",
                parent / "erw_code" / "East_River_Workflow_Code_v0.8.2" / "src",
            ]
        )

    for candidate in candidates:
        if (candidate / "east_river_workflow").is_dir():
            sys.path.insert(0, str(candidate))
            return

    raise ImportError(
        "Could not import east_river_workflow and could not locate the v0.8.2 "
        "src tree. Activate the east-river-workflow environment or install the "
        "project with `python -m pip install -e .`."
    )


_bootstrap_east_river_import()

from east_river_workflow import __version__ as ERW_VERSION  # noqa: E402
from east_river_workflow.config import WorkflowConfig, load_config  # noqa: E402
from east_river_publication import configure_matplotlib, save_figure
from east_river_workflow.constants import MODEL_COLORS  # noqa: E402
from east_river_workflow.data_access import decode_noah_times, read_noah_grid  # noqa: E402
from east_river_workflow.grids import GridSpec, load_weights, regrid_area_average  # noqa: E402
from east_river_workflow.utils import open_dataset_robust, save_json, weighted_nanmean  # noqa: E402


# =============================================================================
# Constants
# =============================================================================
SCRIPT_VERSION = "0.1.5"
DEFAULT_CONFIG = Path(__file__).with_name("east_river_config_output1.yaml")
DEFAULT_MODIS_DIR = Path("/bsuhome/zihaozhang/scratch/MODIS-SCF")
DEFAULT_SUBDIR = "modis_fsno_validation_v0.1.5"

DEFAULT_CELL_COVERAGE_THRESHOLD = 0.75
DEFAULT_BASIN_COVERAGE_THRESHOLD = 0.75

MODIS_PRODUCT_PREFIX = "MOD10A1"
MODIS_TILE = "h09v05"
MODIS_COLLECTION = "061"
MODIS_COLOR = "black"  # observational reference; no MODIS color exists in project constants

MODIS_SNOW_SDS = "NDSI_Snow_Cover"
MODIS_BASIC_QA_SDS = "NDSI_Snow_Cover_Basic_QA"
MODIS_ALGORITHM_QA_SDS = "NDSI_Snow_Cover_Algorithm_Flags_QA"
MODIS_GRANULE_SDS = "granule_pnt"
MODIS_GRID_NAME = "MOD_Grid_Snow_500m"

PRIMARY_BASIC_QA = (0, 1)
MODIS_NODATA = 255
MAX_TIME_MATCH_MINUTES = 30.0 + 1.0e-6

_FILENAME_RE = re.compile(
    r"^MOD10A1\.A(?P<year>\d{4})(?P<doy>\d{3})\."
    r"(?P<tile>h\d{2}v\d{2})\.(?P<collection>\d{3})\..*\.hdf$"
)


# =============================================================================
# Settings / paths
# =============================================================================
@dataclass(frozen=True)
class ValidationSettings:
    output_root: Path
    modis_dir: Path
    cell_coverage_threshold: float
    basin_coverage_threshold: float
    overwrite: bool
    max_files: int | None
    skip_figures: bool

    @property
    def cache_dir(self) -> Path:
        return self.output_root / "cache"

    @property
    def tables_dir(self) -> Path:
        return self.output_root / "tables"

    @property
    def diagnostics_dir(self) -> Path:
        return self.output_root / "diagnostics"

    @property
    def figures_dir(self) -> Path:
        return self.output_root / "figures"

    def create_tree(self) -> None:
        for path in (
            self.cache_dir,
            self.tables_dir,
            self.diagnostics_dir,
            self.figures_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Generic helpers
# =============================================================================
def _setup_logger() -> logging.Logger:
    logger = logging.getLogger("modis_fsno_validation")
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
        )
        logger.addHandler(handler)
    return logger


def _water_year(date: pd.Timestamp) -> int:
    return int(date.year + 1 if date.month >= 10 else date.year)


def _weighted_fraction(mask: np.ndarray, weights: np.ndarray) -> float:
    m = np.asarray(mask, dtype=bool)
    w = np.asarray(weights, dtype=float)
    valid_w = np.isfinite(w) & (w > 0)
    denom = float(np.sum(w[valid_w]))
    if denom <= 0:
        return np.nan
    return float(np.sum(w[valid_w & m]) / denom)


def _weighted_mean_masked(
    values: np.ndarray,
    weights: np.ndarray,
    mask: np.ndarray,
) -> float:
    arr = np.where(np.asarray(mask, dtype=bool), np.asarray(values, dtype=float), np.nan)
    return weighted_nanmean(arr, weights)[0]


def _pearson_r(x: np.ndarray, y: np.ndarray) -> float:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    valid = np.isfinite(xx) & np.isfinite(yy)
    if np.count_nonzero(valid) < 2:
        return np.nan
    xx = xx[valid]
    yy = yy[valid]
    if np.nanstd(xx) == 0 or np.nanstd(yy) == 0:
        return np.nan
    return float(np.corrcoef(xx, yy)[0, 1])


def _configure_publication_matplotlib() -> None:
    configure_matplotlib()


def _save_figure(fig, cfg, settings, stem) -> tuple[Path, Path]:
    return save_figure(fig, settings.figures_dir / stem,
                       close=bool(cfg.section("plotting").get("close_after_save", True)))


def _subplot_grid(years: list[int], *, sharey: bool = True):
    ncols = 3
    nrows = int(math.ceil(len(years) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(21, 5.5 * nrows),
        sharey=sharey,
        squeeze=False,
    )
    for ax in axes.flat[len(years):]:
        ax.set_visible(False)
    return fig, axes


def _apply_template_spacing(fig: plt.Figure) -> None:
    fig.subplots_adjust(
        top=0.95,
        bottom=0.15,
        left=0.07,
        hspace=0.30,
        wspace=0.16,
    )


def _style_standard_axis(ax: plt.Axes) -> None:
    """Apply the v3 standard left/bottom spine and major-tick styling."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(
        axis="both",
        which="major",
        width=1.0,
        length=5,
    )


def _format_water_year_axis(ax: plt.Axes, wy: int) -> None:
    ax.set_xlim(pd.Timestamp(wy - 1, 10, 1), pd.Timestamp(wy, 9, 30))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    _style_standard_axis(ax)


def _water_year_plot_frame(
    daily: pd.DataFrame,
    wy: int,
    cfg: WorkflowConfig,
) -> pd.DataFrame:
    """
    Return a complete daily WY axis.

    Missing/filtered MODIS dates remain NaN.  They are not interpolated.
    Figure 1 renders the retained dates as discrete observation markers rather
    than connecting them across cloud/missing gaps.
    """
    start = max(
        pd.Timestamp(wy - 1, 10, 1),
        cfg.analysis_start.normalize(),
    )
    end = min(
        pd.Timestamp(wy, 9, 30),
        cfg.analysis_end.normalize(),
    )
    if end < start:
        return pd.DataFrame(
            columns=[
                "date",
                "modis_fsca_basin",
                "noah_fsno_basin",
                "delta_fsno_noah_minus_modis",
                "primary_usable",
            ]
        )

    frame = pd.DataFrame({"date": pd.date_range(start, end, freq="D")})
    columns = [
        "date",
        "modis_fsca_basin",
        "noah_fsno_basin",
        "delta_fsno_noah_minus_modis",
        "primary_usable",
    ]
    available = daily.loc[daily["water_year"].eq(wy), columns].copy()
    frame = frame.merge(available, on="date", how="left")

    usable = frame["primary_usable"].eq(True)
    for column in [
        "modis_fsca_basin",
        "noah_fsno_basin",
        "delta_fsno_noah_minus_modis",
    ]:
        frame.loc[~usable, column] = np.nan

    frame["primary_usable"] = usable
    return frame


# =============================================================================
# MODIS inventory / metadata
# =============================================================================
def _parse_modis_filename(path: Path) -> dict[str, Any] | None:
    match = _FILENAME_RE.match(path.name)
    if match is None:
        return None

    year = int(match.group("year"))
    doy = int(match.group("doy"))
    date = pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=doy - 1)

    return {
        "path": path,
        "filename": path.name,
        "date": date.normalize(),
        "year": year,
        "doy": doy,
        "tile": match.group("tile"),
        "collection": match.group("collection"),
    }


def _build_modis_inventory(
    cfg: WorkflowConfig,
    settings: ValidationSettings,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []

    for path in sorted(settings.modis_dir.glob("MOD10A1*.hdf")):
        parsed = _parse_modis_filename(path)
        if parsed is None:
            rows.append(
                {
                    "filename": path.name,
                    "path": str(path),
                    "date": pd.NaT,
                    "tile": None,
                    "collection": None,
                    "parse_ok": False,
                    "in_analysis_window": False,
                    "status": "unparsed_filename",
                }
            )
            continue

        date = parsed["date"]
        in_window = (
            cfg.analysis_start.normalize()
            <= date
            <= cfg.analysis_end.normalize()
        )

        status = "usable_candidate"
        if parsed["tile"] != MODIS_TILE:
            status = "wrong_tile"
        elif parsed["collection"] != MODIS_COLLECTION:
            status = "wrong_collection"
        elif not in_window:
            status = "outside_analysis_window"

        rows.append(
            {
                "filename": parsed["filename"],
                "path": str(parsed["path"]),
                "date": date,
                "tile": parsed["tile"],
                "collection": parsed["collection"],
                "parse_ok": True,
                "in_analysis_window": in_window,
                "status": status,
            }
        )

    inventory = pd.DataFrame(rows)
    if inventory.empty:
        raise FileNotFoundError(
            f"No MOD10A1 HDF files found in {settings.modis_dir}"
        )

    candidates = inventory[inventory["status"].eq("usable_candidate")].copy()
    candidates["date"] = pd.to_datetime(candidates["date"])

    duplicate = candidates["date"].duplicated(keep=False)
    if duplicate.any():
        dup = (
            candidates.loc[duplicate, ["date", "filename"]]
            .sort_values(["date", "filename"])
        )
        raise ValueError(
            "Duplicate MOD10A1 dates detected inside the analysis window. "
            "The script will not choose silently:\n"
            + dup.to_string(index=False)
        )

    expected = pd.DataFrame(
        {
            "date": pd.date_range(
                cfg.analysis_start.normalize(),
                cfg.analysis_end.normalize(),
                freq="D",
            )
        }
    )
    alignment = expected.merge(
        candidates[["date", "filename", "path"]],
        on="date",
        how="left",
    )
    alignment["modis_available"] = alignment["filename"].notna()

    if settings.max_files is not None:
        available_idx = alignment.index[alignment["modis_available"]]
        keep_idx = set(available_idx[: settings.max_files].tolist())
        alignment["selected_for_run"] = [
            bool(i in keep_idx) for i in alignment.index
        ]
    else:
        alignment["selected_for_run"] = alignment["modis_available"]

    return inventory, alignment


def _split_metadata_csv(value: str | None) -> list[str]:
    if value is None:
        return []
    return [item.strip() for item in str(value).split(",")]


def _container_tags(path: Path) -> dict[str, str]:
    """
    Read HDF container metadata.

    Rasterio/GDAL normally exposes the metadata through tags(). If the required
    granule arrays are absent, fail loudly rather than guessing.
    """
    with rasterio.open(path) as ds:
        tags = dict(ds.tags())

        # Some GDAL builds place metadata in named domains. Merge all domains
        # without overwriting a key already present in the default domain.
        try:
            domains = ds.tag_namespaces()
        except Exception:
            domains = []

        for namespace in domains:
            try:
                for key, value in ds.tags(ns=namespace).items():
                    tags.setdefault(key, value)
            except Exception:
                continue

    return tags


def _granule_time_mapping(
    path: Path,
    pointer_values: np.ndarray,
) -> dict[int, dict[str, Any]]:
    """
    Resolve granule_pnt values to source-swath acquisition intervals.

    Mapping rule:
      pointer value -> position where that value occurs in GRANULEPOINTERARRAY
      -> same position in beginning/ending datetime arrays.

    Fixed-size pointer arrays may have trailing -1 padding beyond the declared
    NUMBEROFINPUTGRANULES. Only that verified trailing padding is removed;
    unused input granules (-1 within the declared inputs) retain their positions.
    Any missing or ambiguous mapping raises an exception.
    """
    tags = _container_tags(path)

    required = [
        "GRANULEPOINTERARRAY",
        "GRANULEBEGINNINGDATETIMEARRAY",
        "GRANULEENDINGDATETIMEARRAY",
    ]
    missing = [key for key in required if key not in tags]
    if missing:
        raise KeyError(
            f"{path.name}: required MOD10A1 metadata absent: {missing}. "
            "No acquisition-time guess will be made."
        )

    raw_pointer = _split_metadata_csv(tags["GRANULEPOINTERARRAY"])
    raw_begin = _split_metadata_csv(tags["GRANULEBEGINNINGDATETIMEARRAY"])
    raw_end = _split_metadata_csv(tags["GRANULEENDINGDATETIMEARRAY"])

    if any(v == "" for values in (raw_pointer, raw_begin, raw_end) for v in values):
        raise ValueError(f"{path.name}: empty MODIS pointer/time metadata entry.")
    pointer_numbers = [float(v) for v in raw_pointer]
    if any(not np.isfinite(v) or v != int(v) for v in pointer_numbers):
        raise ValueError(f"{path.name}: non-integral MODIS pointer metadata.")
    pointers = [int(v) for v in pointer_numbers]
    n_times = len(raw_begin)
    counts = (
        f"pointers={len(pointers)}, beginning={n_times}, ending={len(raw_end)}, "
        f"NUMBEROFINPUTGRANULES={tags.get('NUMBEROFINPUTGRANULES', 'absent')}"
    )
    if n_times == 0 or len(raw_end) != n_times:
        raise ValueError(f"{path.name}: inconsistent MODIS time arrays ({counts}).")
    declared_count = tags.get("NUMBEROFINPUTGRANULES")
    if declared_count is not None:
        try:
            n_inputs = int(declared_count)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path.name}: invalid input-granule count ({counts}).") from exc
        if n_inputs != n_times:
            raise ValueError(f"{path.name}: input count conflicts with time arrays ({counts}).")
    if len(pointers) != n_times:
        # MOD10A1 sample: four inputs, 100 pointer slots, 96 trailing -1s.
        # Never filter all -1s: the first input can also be unused, and removing
        # it would shift every remaining pointer onto the wrong timestamp.
        if (
            declared_count is None
            or len(pointers) < n_times
            or any(p != -1 for p in pointers[n_times:])
        ):
            raise ValueError(
                f"{path.name}: unverified MODIS pointer/time length mismatch ({counts}); "
                "only trailing -1 padding beyond the declared inputs is supported."
            )
        pointers = pointers[:n_times]
    begin = [pd.Timestamp(v) for v in raw_begin if v != ""]
    end = [pd.Timestamp(v) for v in raw_end if v != ""]

    mapping: dict[int, dict[str, Any]] = {}

    for pointer in sorted(int(v) for v in np.unique(pointer_values)):
        if pointer == MODIS_NODATA:
            continue

        positions = [
            i
            for i, p in enumerate(pointers)
            if p == pointer and i < len(begin) and i < len(end)
        ]

        if len(positions) != 1:
            raise ValueError(
                f"{path.name}: granule_pnt={pointer} maps to {len(positions)} "
                "metadata positions; expected exactly one. "
                "No acquisition-time guess will be made."
            )

        idx = positions[0]
        b = begin[idx]
        e = end[idx]

        # MODIS metadata is UTC (Z). Normalize to naive UTC for direct matching
        # to the project's confirmed UTC Noah-MP clock.
        if b.tzinfo is not None:
            b = b.tz_convert("UTC").tz_localize(None)
        if e.tzinfo is not None:
            e = e.tz_convert("UTC").tz_localize(None)

        if pd.isna(b) or pd.isna(e) or e < b:
            raise ValueError(
                f"{path.name}: granule_pnt={pointer} has end before begin."
            )

        midpoint = b + (e - b) / 2

        mapping[pointer] = {
            "metadata_index": idx,
            "begin_utc": b,
            "end_utc": e,
            "midpoint_utc": midpoint,
        }

    return mapping


# =============================================================================
# MODIS SDS / source window
# =============================================================================
def _sds_path(hdf_path: Path, name: str) -> str:
    return (
        f"HDF4_EOS:EOS_GRID:{hdf_path}:"
        f"{MODIS_GRID_NAME}:{name}"
    )


def _clip_window_to_dataset(window: Window, width: int, height: int) -> Window:
    col0 = max(0, int(math.floor(window.col_off)))
    row0 = max(0, int(math.floor(window.row_off)))
    col1 = min(width, int(math.ceil(window.col_off + window.width)))
    row1 = min(height, int(math.ceil(window.row_off + window.height)))

    if col1 <= col0 or row1 <= row0:
        raise ValueError("Requested MODIS window does not intersect the tile.")

    return Window(col0, row0, col1 - col0, row1 - row0)


def _modis_window_for_noah_domain(
    snow_sds: str,
    noah_grid: GridSpec,
) -> tuple[Window, Affine, CRS, int, int]:
    with rasterio.open(snow_sds) as src:
        src_crs = CRS.from_user_input(src.crs)

        left, bottom, right, top = transform_bounds(
            noah_grid.crs,
            src_crs,
            *noah_grid.bounds,
            densify_pts=21,
        )

        # One source pixel of padding prevents boundary-rounding loss.
        pad_x = abs(float(src.transform.a))
        pad_y = abs(float(src.transform.e))
        left -= pad_x
        right += pad_x
        bottom -= pad_y
        top += pad_y

        window = from_bounds(
            left,
            bottom,
            right,
            top,
            transform=src.transform,
        )
        window = _clip_window_to_dataset(window, src.width, src.height)
        transform = src.window_transform(window)

    return (
        window,
        transform,
        src_crs,
        int(window.width),
        int(window.height),
    )


def _grid_from_window(
    transform: Affine,
    crs: CRS,
    width: int,
    height: int,
) -> GridSpec:
    x = transform.c + (np.arange(width) + 0.5) * transform.a
    y = transform.f + (np.arange(height) + 0.5) * transform.e

    return GridSpec(
        name="MODIS-MOD10A1-window",
        crs=crs,
        transform=transform,
        width=width,
        height=height,
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
    )


def _read_modis_window(
    hdf_path: Path,
    noah_grid: GridSpec,
) -> dict[str, Any]:
    snow_sds = _sds_path(hdf_path, MODIS_SNOW_SDS)
    qa_sds = _sds_path(hdf_path, MODIS_BASIC_QA_SDS)
    flags_sds = _sds_path(hdf_path, MODIS_ALGORITHM_QA_SDS)
    granule_sds = _sds_path(hdf_path, MODIS_GRANULE_SDS)

    (
        window,
        transform,
        crs,
        width,
        height,
    ) = _modis_window_for_noah_domain(snow_sds, noah_grid)

    with rasterio.open(snow_sds) as ds:
        reference_width = ds.width
        reference_height = ds.height
        reference_transform = ds.transform
        reference_crs = CRS.from_user_input(ds.crs)
        snow = ds.read(1, window=window)

    def read_aligned(sds: str, label: str) -> np.ndarray:
        with rasterio.open(sds) as ds:
            current_crs = CRS.from_user_input(ds.crs)
            if ds.width != reference_width or ds.height != reference_height:
                raise ValueError(
                    f"{hdf_path.name}: {label} raster dimensions "
                    f"{(ds.height, ds.width)} differ from NDSI_Snow_Cover "
                    f"{(reference_height, reference_width)}."
                )
            if current_crs != reference_crs:
                raise ValueError(
                    f"{hdf_path.name}: {label} CRS differs from NDSI_Snow_Cover."
                )
            if not np.allclose(
                tuple(ds.transform),
                tuple(reference_transform),
                rtol=0.0,
                atol=1.0e-9,
            ):
                raise ValueError(
                    f"{hdf_path.name}: {label} transform differs from NDSI_Snow_Cover."
                )
            return ds.read(1, window=window)

    qa = read_aligned(qa_sds, MODIS_BASIC_QA_SDS)
    flags = read_aligned(flags_sds, MODIS_ALGORITHM_QA_SDS)
    granule = read_aligned(granule_sds, MODIS_GRANULE_SDS)

    shape = (height, width)
    for name, arr in {
        "snow": snow,
        "qa": qa,
        "flags": flags,
        "granule": granule,
    }.items():
        if arr.shape != shape:
            raise ValueError(
                f"{hdf_path.name}: {name} shape {arr.shape} != expected {shape}"
            )

    if crs != reference_crs:
        raise ValueError(
            f"{hdf_path.name}: internal MODIS source CRS changed while building the read window."
        )

    return {
        "grid": _grid_from_window(transform, crs, width, height),
        "snow": snow,
        "qa": qa,
        "flags": flags,
        "granule": granule,
    }


def _reproject_fraction(
    mask: np.ndarray,
    source_grid: GridSpec,
    destination_grid: GridSpec,
) -> np.ndarray:
    """Area-average a 0/1 mask to destination area fraction."""
    source = np.asarray(mask, dtype=np.float32)
    destination = np.zeros(destination_grid.shape, dtype=np.float32)

    reproject(
        source=source,
        destination=destination,
        src_transform=source_grid.transform,
        src_crs=source_grid.crs,
        src_nodata=None,
        dst_transform=destination_grid.transform,
        dst_crs=destination_grid.crs,
        dst_nodata=0.0,
        resampling=Resampling.average,
        init_dest_nodata=True,
    )

    return np.asarray(destination, dtype=float)


# =============================================================================
# Noah grid / hourly FSNO
# =============================================================================
def _load_noah_static(
    cfg: WorkflowConfig,
    settings: ValidationSettings,
    logger: logging.Logger,
):
    noah_grid, noah_dem, noah_lat, noah_lon, noah_diag = read_noah_grid(
        cfg.data["paths"]["noahmp_geo_file"],
        cfg.section("noahmp"),
    )

    main_weight_path = cfg.output_dir / "cache" / "east_river_fraction_noahmp.npz"
    if not main_weight_path.exists():
        raise FileNotFoundError(
            "Validated East River Noah fractional-weight cache is absent:\n"
            f"  {main_weight_path}\n"
            "Run the main v0.8.2 workflow first. This script will not silently "
            "rebuild a different watershed mask."
        )

    weights = load_weights(main_weight_path, noah_grid)

    static_path = settings.cache_dir / "modis_fsno_static_noahgrid.npz"
    np.savez_compressed(
        static_path,
        watershed_fraction=np.asarray(weights, dtype=np.float32),
        noah_dem_m=np.asarray(noah_dem, dtype=np.float32),
        noah_lat=np.asarray(noah_lat, dtype=np.float32),
        noah_lon=np.asarray(noah_lon, dtype=np.float32),
        x=np.asarray(noah_grid.x, dtype=float),
        y=np.asarray(noah_grid.y, dtype=float),
    )

    logger.info("Noah grid shape: %s", noah_grid.shape)
    logger.info(
        "East River weighted Noah cells: %d",
        int(np.count_nonzero(np.asarray(weights) > 0)),
    )

    return noah_grid, weights, noah_dem, noah_lat, noah_lon, noah_diag


def _open_noah_fsno(
    cfg: WorkflowConfig,
    noah_grid_shape: tuple[int, int],
):
    source = Path(cfg.data["paths"]["noahmp_output_file"]).expanduser()

    ds_context = open_dataset_robust(source, decode_times=False)
    ds = ds_context.__enter__()

    fsno_name = "FSNO"
    if fsno_name not in ds:
        ds_context.__exit__(None, None, None)
        raise KeyError(
            f"Required Noah-MP FSNO variable {fsno_name!r} absent from {source}"
        )

    fsno_var = ds[fsno_name]
    if "Time" not in fsno_var.dims:
        ds_context.__exit__(None, None, None)
        raise ValueError(
            f"FSNO must have Time dimension; dims={fsno_var.dims}"
        )

    times = pd.DatetimeIndex(
        decode_noah_times(ds, cfg.section("noahmp"))
    )

    if times.has_duplicates:
        ds_context.__exit__(None, None, None)
        raise ValueError("Decoded Noah-MP hourly timestamps contain duplicates.")

    if not times.is_monotonic_increasing:
        ds_context.__exit__(None, None, None)
        raise ValueError("Decoded Noah-MP timestamps are not monotonic increasing.")

    if times.hasnans:
        ds_context.__exit__(None, None, None)
        raise ValueError("Decoded Noah-MP timestamps contain NaT values.")

    section = cfg.section("noahmp")
    if len(times) > 1:
        expected_seconds = float(section["expected_timestep_hours"]) * 3600.0
        diffs = np.asarray(
            (times[1:] - times[:-1]).total_seconds(),
            dtype=float,
        )
        bad = np.flatnonzero(
            ~np.isclose(diffs, expected_seconds, rtol=0.0, atol=1.0)
        )
        if len(bad):
            ds_context.__exit__(None, None, None)
            raise ValueError(
                f"Decoded Noah-MP time axis contains {len(bad)} unexpected "
                "time gaps."
            )

    expected_count = section.get("expected_record_count")
    if expected_count not in {None, ""} and len(times) != int(expected_count):
        ds_context.__exit__(None, None, None)
        raise ValueError(
            f"Noah-MP record count is {len(times)}, expected {expected_count}."
        )

    if cfg.analysis_start < times[0] or cfg.analysis_end > times[-1]:
        ds_context.__exit__(None, None, None)
        raise ValueError(
            f"Analysis window {cfg.analysis_start}..{cfg.analysis_end} is outside "
            f"decoded Noah-MP times {times[0]}..{times[-1]}."
        )

    sample_shape = np.asarray(
        fsno_var.isel(Time=0).values
    ).squeeze().shape
    if sample_shape != noah_grid_shape:
        ds_context.__exit__(None, None, None)
        raise ValueError(
            f"FSNO spatial shape {sample_shape} != Noah grid shape {noah_grid_shape}."
        )

    return ds_context, ds, fsno_var, times


def _nearest_noah_index(
    target_utc: pd.Timestamp,
    noah_times: pd.DatetimeIndex,
) -> tuple[int, pd.Timestamp, float]:
    target = pd.Timestamp(target_utc)

    pos = int(noah_times.searchsorted(target))
    candidates: list[int] = []

    if pos < len(noah_times):
        candidates.append(pos)
    if pos > 0:
        candidates.append(pos - 1)

    if not candidates:
        raise ValueError("No Noah-MP timestamps available.")

    offsets = [
        abs((noah_times[i] - target).total_seconds()) / 60.0
        for i in candidates
    ]
    best_local = int(np.argmin(offsets))
    idx = candidates[best_local]
    offset = float(offsets[best_local])

    if offset > MAX_TIME_MATCH_MINUTES:
        raise ValueError(
            f"Nearest Noah-MP hour to MODIS midpoint {target} is "
            f"{noah_times[idx]} ({offset:.2f} min away), exceeding the "
            f"{MAX_TIME_MATCH_MINUTES:.1f}-minute hourly-matching tolerance."
        )

    return idx, pd.Timestamp(noah_times[idx]), offset


def _read_fsno_slice(
    fsno_var,
    index: int,
    expected_shape: tuple[int, int],
) -> np.ndarray:
    fsno = np.asarray(fsno_var.isel(Time=int(index)).values, dtype=float).squeeze()

    if fsno.shape != expected_shape:
        raise ValueError(
            f"FSNO shape {fsno.shape} != Noah grid shape {expected_shape}"
        )

    finite = fsno[np.isfinite(fsno)]
    if finite.size:
        vmin = float(np.min(finite))
        vmax = float(np.max(finite))
        if vmin < -1.0e-6 or vmax > 1.0 + 1.0e-6:
            raise ValueError(
                f"Noah-MP FSNO outside [0,1]: min={vmin}, max={vmax}"
            )

    return np.where(
        np.isfinite(fsno),
        np.clip(fsno, 0.0, 1.0),
        np.nan,
    )


# =============================================================================
# Per-day processing
# =============================================================================
def _mix_model_granules(pointer_coverage, pointer_fsno, shape):
    """Use every contributing observation swath, or exclude the whole cell."""
    numerator = np.zeros(shape, dtype=float)
    total = np.zeros(shape, dtype=float)
    complete = np.ones(shape, dtype=bool)
    for pointer, coverage in pointer_coverage.items():
        coverage = np.asarray(coverage, dtype=float)
        fsno = np.asarray(pointer_fsno[pointer], dtype=float)
        complete &= (coverage <= 0) | np.isfinite(fsno)
        valid = (coverage > 0) & np.isfinite(fsno)
        numerator[valid] += fsno[valid] * coverage[valid]
        total[valid] += coverage[valid]
    out = np.full(shape, np.nan, dtype=float)
    good = complete & (total > 0)
    out[good] = numerator[good] / total[good]
    return out


def _process_one_day(
    date: pd.Timestamp,
    hdf_path: Path,
    cfg: WorkflowConfig,
    settings: ValidationSettings,
    noah_grid: GridSpec,
    weights: np.ndarray,
    fsno_var,
    noah_times: pd.DatetimeIndex,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, np.ndarray]]:
    modis = _read_modis_window(hdf_path, noah_grid)

    source_grid = modis["grid"]
    snow = np.asarray(modis["snow"])
    qa = np.asarray(modis["qa"])
    flags = np.asarray(modis["flags"])
    granule = np.asarray(modis["granule"])

    # Primary MODIS observation mask.
    valid_value = np.isfinite(snow) & (snow >= 0) & (snow <= 100)
    valid_qa = np.isin(qa, PRIMARY_BASIC_QA)
    valid_pointer = granule != MODIS_NODATA
    valid_obs = valid_value & valid_qa & valid_pointer

    snow_binary = np.full(snow.shape, np.nan, dtype=np.float32)
    snow_binary[valid_obs] = (snow[valid_obs] > 0).astype(np.float32)

    modis_fsca, valid_coverage = regrid_area_average(
        snow_binary,
        source_grid,
        noah_grid,
        minimum_valid_coverage=settings.cell_coverage_threshold,
    )

    # Source granules actually contributing valid primary observations.
    pointer_values = np.unique(granule[valid_obs])
    if pointer_values.size == 0:
        pointer_values = np.array([], dtype=np.uint8)
        granule_mapping: dict[int, dict[str, Any]] = {}
    else:
        granule_mapping = _granule_time_mapping(hdf_path, pointer_values)
        for pointer, info in granule_mapping.items():
            begin_date = pd.Timestamp(info["begin_utc"]).normalize()
            end_date = pd.Timestamp(info["end_utc"]).normalize()
            if begin_date != date.normalize() or end_date != date.normalize():
                raise ValueError(
                    f"{hdf_path.name}: granule_pnt={pointer} acquisition interval "
                    f"{info['begin_utc']}..{info['end_utc']} does not match the "
                    f"MOD10A1 file date {date.date()}. No cross-date time match "
                    "will be made silently."
                )

    # For each granule pointer, calculate its valid MODIS source-area fraction
    # within every Noah cell and read the corresponding hourly FSNO.
    pointer_coverage: dict[int, np.ndarray] = {}
    pointer_fsno: dict[int, np.ndarray] = {}
    time_rows: list[dict[str, Any]] = []

    for pointer in sorted(int(v) for v in pointer_values):
        info = granule_mapping[pointer]

        mask = valid_obs & (granule == pointer)
        pcoverage = _reproject_fraction(mask, source_grid, noah_grid)

        noah_idx, matched_time, offset_minutes = _nearest_noah_index(
            info["midpoint_utc"],
            noah_times,
        )
        fsno = _read_fsno_slice(
            fsno_var,
            noah_idx,
            noah_grid.shape,
        )

        pointer_coverage[pointer] = pcoverage
        pointer_fsno[pointer] = fsno

        time_rows.append(
            {
                "pointer": pointer,
                "metadata_index": int(info["metadata_index"]),
                "begin_utc": info["begin_utc"],
                "end_utc": info["end_utc"],
                "midpoint_utc": info["midpoint_utc"],
                "noah_time_utc": matched_time,
                "time_offset_minutes": offset_minutes,
            }
        )

    # Combine time-matched FSNO states within each Noah cell using the same
    # valid MODIS source-area fractions that contributed to that cell.
    fsno_matched = _mix_model_granules(pointer_coverage, pointer_fsno, noah_grid.shape)

    # Apply the exact same minimum valid MODIS source-area criterion to the
    # model counterpart.
    fsno_matched[valid_coverage < settings.cell_coverage_threshold] = np.nan

    watershed = np.asarray(weights, dtype=float) > 0
    paired = (
        watershed
        & np.isfinite(modis_fsca)
        & np.isfinite(fsno_matched)
    )

    paired_area_fraction = _weighted_fraction(paired, weights)
    watershed_weight_sum = float(
        np.sum(np.asarray(weights, dtype=float)[watershed])
    )
    modis_valid_source_watershed_fraction = (
        float(
            np.sum(
                np.clip(valid_coverage[watershed], 0.0, 1.0)
                * np.asarray(weights, dtype=float)[watershed]
            )
        )
        / watershed_weight_sum
        if watershed_weight_sum > 0
        else np.nan
    )
    modis_basin = _weighted_mean_masked(modis_fsca, weights, paired)
    fsno_basin = _weighted_mean_masked(fsno_matched, weights, paired)

    primary_usable = bool(
        np.isfinite(paired_area_fraction)
        and paired_area_fraction >= settings.basin_coverage_threshold
    )

    difference = fsno_matched - modis_fsca
    difference[~paired] = np.nan

    # Dominant pointer over valid East River Noah-cell source area.
    pointer_watershed_area: dict[int, float] = {}
    for pointer, pcov in pointer_coverage.items():
        p = np.asarray(pcov, dtype=float)
        pointer_watershed_area[pointer] = float(
            np.nansum(p[watershed] * np.asarray(weights, dtype=float)[watershed])
        )

    if pointer_watershed_area:
        dominant_pointer = max(
            pointer_watershed_area,
            key=pointer_watershed_area.get,
        )
        total_pointer_area = float(sum(pointer_watershed_area.values()))
        dominant_fraction = (
            pointer_watershed_area[dominant_pointer] / total_pointer_area
            if total_pointer_area > 0
            else np.nan
        )
        dominant_info = next(
            row for row in time_rows if row["pointer"] == dominant_pointer
        )
    else:
        dominant_pointer = np.nan
        dominant_fraction = np.nan
        dominant_info = {
            "begin_utc": pd.NaT,
            "end_utc": pd.NaT,
            "midpoint_utc": pd.NaT,
            "noah_time_utc": pd.NaT,
            "time_offset_minutes": np.nan,
        }

    daily_row = {
        "date": date.normalize(),
        "water_year": _water_year(date),
        "modis_file": hdf_path.name,
        "modis_fsca_basin": modis_basin,
        "noah_fsno_basin": fsno_basin,
        "delta_fsno_noah_minus_modis": (
            fsno_basin - modis_basin
            if np.isfinite(fsno_basin) and np.isfinite(modis_basin)
            else np.nan
        ),
        "paired_watershed_area_fraction": paired_area_fraction,
        "modis_valid_source_watershed_fraction": (
            modis_valid_source_watershed_fraction
        ),
        "primary_usable": primary_usable,
        "paired_noah_cell_count": int(np.count_nonzero(paired)),
        "dominant_granule_pointer": dominant_pointer,
        "dominant_pointer_fraction": dominant_fraction,
        "dominant_modis_begin_utc": dominant_info["begin_utc"],
        "dominant_modis_end_utc": dominant_info["end_utc"],
        "dominant_modis_midpoint_utc": dominant_info["midpoint_utc"],
        "dominant_noah_time_utc": dominant_info["noah_time_utc"],
        "dominant_time_offset_minutes": dominant_info["time_offset_minutes"],
        "granule_pointer_count": int(len(pointer_values)),
    }

    qa_row: dict[str, Any] = {
        "date": date.normalize(),
        "water_year": _water_year(date),
        "modis_file": hdf_path.name,
        "source_window_pixels": int(snow.size),
        "valid_value_pixels": int(np.count_nonzero(valid_value)),
        "basic_qa_0_pixels": int(np.count_nonzero(qa == 0)),
        "basic_qa_1_pixels": int(np.count_nonzero(qa == 1)),
        "basic_qa_other_pixels": int(
            np.count_nonzero(~np.isin(qa, PRIMARY_BASIC_QA))
        ),
        "primary_valid_pixels": int(np.count_nonzero(valid_obs)),
        "primary_valid_fraction_source_window": float(np.mean(valid_obs)),
        "algorithm_flags_unique_count": int(np.unique(flags).size),
        "algorithm_flags_fill_on_primary_valid_pixels": int(
            np.count_nonzero(valid_obs & (flags == MODIS_NODATA))
        ),
        "pointer_count_primary_valid": int(pointer_values.size),
    }

    # Record compact algorithm flag bit frequencies among primary-valid pixels.
    # Fill=255 is excluded because otherwise every bit would appear spuriously set.
    flag_diag_mask = valid_obs & (flags != MODIS_NODATA)
    if np.any(flag_diag_mask):
        flag_valid = flags[flag_diag_mask].astype(np.uint8)
        for bit in range(8):
            qa_row[f"algorithm_flag_bit{bit}_fraction"] = float(
                np.mean(((flag_valid >> bit) & 1) == 1)
            )
    else:
        for bit in range(8):
            qa_row[f"algorithm_flag_bit{bit}_fraction"] = np.nan

    # Put time mapping into a compact JSON string for the daily alignment table.
    time_mapping_json: list[dict[str, Any]] = []
    for row in time_rows:
        time_mapping_json.append(
            {
                "pointer": int(row["pointer"]),
                "metadata_index": int(row["metadata_index"]),
                "begin_utc": str(pd.Timestamp(row["begin_utc"])),
                "end_utc": str(pd.Timestamp(row["end_utc"])),
                "midpoint_utc": str(pd.Timestamp(row["midpoint_utc"])),
                "noah_time_utc": str(pd.Timestamp(row["noah_time_utc"])),
                "time_offset_minutes": float(row["time_offset_minutes"]),
                "watershed_pointer_fraction": (
                    pointer_watershed_area.get(int(row["pointer"]), 0.0)
                    / sum(pointer_watershed_area.values())
                    if sum(pointer_watershed_area.values()) > 0
                    else np.nan
                ),
            }
        )

    daily_row["granule_time_mapping_json"] = json.dumps(time_mapping_json)

    arrays = {
        "modis_fsca": np.asarray(modis_fsca, dtype=np.float32),
        "modis_valid_coverage": np.asarray(valid_coverage, dtype=np.float32),
        "noah_fsno": np.asarray(fsno_matched, dtype=np.float32),
        "delta_fsno_noah_minus_modis": np.asarray(difference, dtype=np.float32),
        "paired_mask": np.asarray(paired, dtype=np.uint8),
    }

    return daily_row, qa_row, arrays


# =============================================================================
# Cell-statistics accumulator
# =============================================================================
@dataclass
class CellAccumulator:
    n: np.ndarray
    sum_diff: np.ndarray
    sum_abs_diff: np.ndarray
    sum_sq_diff: np.ndarray
    sum_modis: np.ndarray
    sum_noah: np.ndarray
    sum_modis2: np.ndarray
    sum_noah2: np.ndarray
    sum_cross: np.ndarray

    @classmethod
    def create(cls, shape: tuple[int, int]) -> "CellAccumulator":
        zeros = lambda: np.zeros(shape, dtype=np.float64)
        return cls(
            n=np.zeros(shape, dtype=np.int64),
            sum_diff=zeros(),
            sum_abs_diff=zeros(),
            sum_sq_diff=zeros(),
            sum_modis=zeros(),
            sum_noah=zeros(),
            sum_modis2=zeros(),
            sum_noah2=zeros(),
            sum_cross=zeros(),
        )

    def update(
        self,
        modis: np.ndarray,
        noah: np.ndarray,
        paired: np.ndarray,
        *,
        use_day: bool,
    ) -> None:
        if not use_day:
            return

        x = np.asarray(modis, dtype=float)
        y = np.asarray(noah, dtype=float)
        m = (
            np.asarray(paired, dtype=bool)
            & np.isfinite(x)
            & np.isfinite(y)
        )

        d = y - x

        self.n[m] += 1
        self.sum_diff[m] += d[m]
        self.sum_abs_diff[m] += np.abs(d[m])
        self.sum_sq_diff[m] += d[m] ** 2
        self.sum_modis[m] += x[m]
        self.sum_noah[m] += y[m]
        self.sum_modis2[m] += x[m] ** 2
        self.sum_noah2[m] += y[m] ** 2
        self.sum_cross[m] += x[m] * y[m]

    def to_frame(
        self,
        weights: np.ndarray,
        dem: np.ndarray,
        lat: np.ndarray,
        lon: np.ndarray,
    ) -> pd.DataFrame:
        w = np.asarray(weights, dtype=float)
        watershed = np.isfinite(w) & (w > 0)

        n = self.n.astype(float)
        bias = np.full(n.shape, np.nan, dtype=float)
        mae = np.full(n.shape, np.nan, dtype=float)
        rmse = np.full(n.shape, np.nan, dtype=float)
        r = np.full(n.shape, np.nan, dtype=float)

        valid = n > 0
        bias[valid] = self.sum_diff[valid] / n[valid]
        mae[valid] = self.sum_abs_diff[valid] / n[valid]
        rmse[valid] = np.sqrt(self.sum_sq_diff[valid] / n[valid])

        corr_valid = n >= 2
        numerator = (
            n * self.sum_cross
            - self.sum_modis * self.sum_noah
        )
        denom_x = (
            n * self.sum_modis2
            - self.sum_modis ** 2
        )
        denom_y = (
            n * self.sum_noah2
            - self.sum_noah ** 2
        )
        denom = np.sqrt(np.maximum(denom_x, 0) * np.maximum(denom_y, 0))

        ok = corr_valid & (denom > 0)
        r[ok] = numerator[ok] / denom[ok]

        rows, cols = np.where(watershed)

        return pd.DataFrame(
            {
                "row": rows,
                "col": cols,
                "latitude": np.asarray(lat, dtype=float)[rows, cols],
                "longitude": np.asarray(lon, dtype=float)[rows, cols],
                "elevation_m": np.asarray(dem, dtype=float)[rows, cols],
                "watershed_fraction": w[rows, cols],
                "paired_n": self.n[rows, cols],
                "bias_noah_minus_modis": bias[rows, cols],
                "mae": mae[rows, cols],
                "rmse": rmse[rows, cols],
                "pearson_r": r[rows, cols],
            }
        )


# =============================================================================
# WY cache / summaries
# =============================================================================
def _save_wy_cache(
    settings: ValidationSettings,
    wy: int,
    dates: list[pd.Timestamp],
    arrays: dict[str, list[np.ndarray]],
) -> Path:
    path = settings.cache_dir / f"modis_fsno_wy{wy}.npz"

    payload: dict[str, np.ndarray] = {
        "dates": np.asarray(
            [pd.Timestamp(d).strftime("%Y-%m-%d") for d in dates],
            dtype="U10",
        )
    }

    for key, values in arrays.items():
        if len(values) == 0:
            continue
        stack = np.stack(values, axis=0)
        if stack.dtype == np.uint8:
            payload[key] = stack
        else:
            payload[key] = stack.astype(np.float32)

    np.savez_compressed(path, **payload)
    return path


def _annual_summary(daily: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for wy in sorted(daily["water_year"].dropna().astype(int).unique()):
        all_days = daily[daily["water_year"].eq(wy)].copy()
        sub = all_days[all_days["primary_usable"]].copy()

        x = sub["modis_fsca_basin"].to_numpy(dtype=float)
        y = sub["noah_fsno_basin"].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)

        if np.any(valid):
            diff = y[valid] - x[valid]
            bias = float(np.mean(diff))
            mae = float(np.mean(np.abs(diff)))
            rmse = float(np.sqrt(np.mean(diff ** 2)))
            r = _pearson_r(x[valid], y[valid])
        else:
            bias = mae = rmse = r = np.nan

        rows.append(
            {
                "water_year": int(wy),
                "modis_file_days": int(len(all_days)),
                "primary_usable_days": int(np.count_nonzero(valid)),
                "mean_paired_watershed_area_fraction": float(
                    sub["paired_watershed_area_fraction"].mean()
                ) if len(sub) else np.nan,
                "mean_modis_valid_source_watershed_fraction": float(
                    sub["modis_valid_source_watershed_fraction"].mean()
                ) if len(sub) else np.nan,
                "mean_modis_fsca_basin": float(np.nanmean(x)) if np.any(np.isfinite(x)) else np.nan,
                "mean_noah_fsno_basin": float(np.nanmean(y)) if np.any(np.isfinite(y)) else np.nan,
                "bias_noah_minus_modis": bias,
                "mae": mae,
                "rmse": rmse,
                "pearson_r": r,
            }
        )

    return pd.DataFrame(rows)


# =============================================================================
# Plotting
# =============================================================================
def _plot_yearly_comparison(
    cfg: WorkflowConfig,
    settings: ValidationSettings,
    daily: pd.DataFrame,
    years: list[int],
) -> tuple[Path, Path]:
    """
    Plot paired MODIS and Noah observations as discrete markers.

    MODIS availability is cloud/QA dependent, so connecting observations across
    missing dates would imply temporal information that does not exist.  Noah is
    plotted only at the same MODIS observation dates/times to preserve paired
    support.
    """
    fig, axes = _subplot_grid(years, sharey=True)

    for ax, wy in zip(axes.flat, years):
        sub = _water_year_plot_frame(daily, wy, cfg)

        # Distinct marker shapes make source identity readable without relying
        # on color alone.
        ax.plot(
            sub["date"],
            sub["modis_fsca_basin"],
            linestyle="none",
            marker="o",
            markersize=4.5,
            markerfacecolor="none",
            markeredgecolor=MODIS_COLOR,
            markeredgewidth=1.2,
            color=MODIS_COLOR,
        )
        ax.plot(
            sub["date"],
            sub["noah_fsno_basin"],
            linestyle="none",
            marker="x",
            markersize=5.0,
            markeredgewidth=1.2,
            color=MODEL_COLORS["Noah-MP"],
        )

        x = sub["modis_fsca_basin"].to_numpy(dtype=float)
        y = sub["noah_fsno_basin"].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)

        if np.any(valid):
            r = _pearson_r(x[valid], y[valid])
            mb = float(np.mean(y[valid] - x[valid]))
            rtxt = "NA" if not np.isfinite(r) else f"{r:.2f}"
            stat_text = (
                f"r = {rtxt}\n"
                f"MB = {mb:+.3f}\n"
                f"n = {int(np.count_nonzero(valid)):,}"
            )
        else:
            stat_text = "r = NA\nMB = NA\nn = 0"

        ax.text(
            0.97,
            0.95,
            stat_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=18,
        )

        ax.set_ylim(0, 1.05)
        ax.set_title(
            f"WY {wy}",
            fontsize=20,
            fontweight="normal",
            loc="center",
            pad=8,
        )
        _format_water_year_axis(ax, wy)

    fig.supylabel("Snow-covered area fraction", x=0.015)
    fig.supxlabel("Month", y=0.095)

    handles = [
        Line2D(
            [],
            [],
            linestyle="none",
            marker="o",
            markersize=6,
            markerfacecolor="none",
            markeredgecolor=MODIS_COLOR,
            markeredgewidth=1.2,
            color=MODIS_COLOR,
            label="MODIS-derived fSCA",
        ),
        Line2D(
            [],
            [],
            linestyle="none",
            marker="x",
            markersize=6,
            markeredgewidth=1.2,
            color=MODEL_COLORS["Noah-MP"],
            label="Noah-MP FSNO",
        ),
    ]

    # v5 standard: boxed legend in the first panel when it does not obscure
    # important data.  Upper-left is generally clear for this watershed plot.
    axes.flat[0].legend(
        handles=handles,
        loc="upper left",
        ncol=1,
        frameon=True,
        facecolor="white",
        edgecolor="0.4",
        framealpha=0.9,
        fancybox=False,
        borderpad=0.5,
        handlelength=2.2,
    )

    _apply_template_spacing(fig)

    return _save_figure(
        fig,
        cfg,
        settings,
        "figure01_modis_fsno_yearly_comparison",
    )


def _plot_yearly_difference(
    cfg: WorkflowConfig,
    settings: ValidationSettings,
    daily: pd.DataFrame,
    years: list[int],
) -> tuple[Path, Path]:
    fig, axes = _subplot_grid(years, sharey=True)

    values = daily.loc[
        daily["primary_usable"],
        "delta_fsno_noah_minus_modis",
    ].to_numpy(dtype=float)
    finite = np.abs(values[np.isfinite(values)])
    lim = (
        max(0.1, float(np.ceil(np.nanpercentile(finite, 99) / 0.1) * 0.1))
        if finite.size
        else 1.0
    )

    for ax, wy in zip(axes.flat, years):
        sub = _water_year_plot_frame(daily, wy, cfg)

        ax.plot(
            sub["date"],
            sub["delta_fsno_noah_minus_modis"],
            linestyle="none",
            marker="o",
            markersize=4.5,
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1.2,
            color="black",
        )
        ax.axhline(0.0, color="0.5", lw=1.0)
        ax.set_ylim(-lim, lim)
        ax.set_title(
            f"WY {wy}",
            fontsize=20,
            fontweight="normal",
            loc="center",
            pad=8,
        )
        _format_water_year_axis(ax, wy)

    fig.supylabel("Noah-MP FSNO - MODIS-derived fSCA", x=0.015)
    fig.supxlabel("Month", y=0.095)
    _apply_template_spacing(fig)

    return _save_figure(
        fig,
        cfg,
        settings,
        "figure02_modis_fsno_yearly_difference",
    )


def _plot_scatter(
    cfg: WorkflowConfig,
    settings: ValidationSettings,
    daily: pd.DataFrame,
    years: list[int],
) -> tuple[Path, Path]:
    fig, axes = _subplot_grid(years, sharey=True)

    for ax, wy in zip(axes.flat, years):
        sub = daily[
            daily["water_year"].eq(wy)
            & daily["primary_usable"]
        ]

        x = sub["modis_fsca_basin"].to_numpy(dtype=float)
        y = sub["noah_fsno_basin"].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)

        xx = x[valid]
        yy = y[valid]

        ax.scatter(
            xx,
            yy,
            s=24,
            facecolors="none",
            edgecolors=MODEL_COLORS["Noah-MP"],
            linewidths=1.2,
        )
        ax.plot([0, 1], [0, 1], color="0.5", lw=1.0)

        r = _pearson_r(xx, yy)
        if len(xx):
            mb = float(np.mean(yy - xx))
            rtxt = "NA" if not np.isfinite(r) else f"{r:.2f}"
            stats_text = (
                f"r = {rtxt}\n"
                f"MB = {mb:+.3f}\n"
                f"n = {len(xx):,}"
            )
        else:
            stats_text = "r = NA\nMB = NA\nn = 0"

        ax.text(
            0.97,
            0.95,
            stats_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=18,
        )
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(
            f"WY {wy}",
            fontsize=20,
            fontweight="normal",
            loc="center",
            pad=8,
        )
        _style_standard_axis(ax)

    fig.supylabel("Noah-MP FSNO", x=0.015)
    fig.supxlabel("MODIS-derived fSCA", y=0.095)
    _apply_template_spacing(fig)

    return _save_figure(
        fig,
        cfg,
        settings,
        "figure03_modis_fsno_scatter",
    )


def _plot_spatial_bias(
    cfg: WorkflowConfig,
    settings: ValidationSettings,
    cell_stats: pd.DataFrame,
    noah_lat: np.ndarray,
    noah_lon: np.ndarray,
    weights: np.ndarray,
) -> tuple[Path, Path]:
    shape = np.asarray(weights).shape
    bias = np.full(shape, np.nan, dtype=float)

    for row in cell_stats.itertuples(index=False):
        if row.paired_n > 0:
            bias[int(row.row), int(row.col)] = row.bias_noah_minus_modis

    bias[np.asarray(weights) <= 0] = np.nan

    finite = np.abs(bias[np.isfinite(bias)])
    lim = (
        max(0.05, float(np.ceil(np.nanpercentile(finite, 99) / 0.05) * 0.05))
        if finite.size
        else 1.0
    )

    fig, ax = plt.subplots(figsize=(10, 8))
    mesh = ax.pcolormesh(
        noah_lon,
        noah_lat,
        bias,
        shading="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim),
    )
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Mean Noah-MP FSNO - MODIS-derived fSCA")

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    _style_standard_axis(ax)

    return _save_figure(
        fig,
        cfg,
        settings,
        "figure04_modis_fsno_spatial_bias",
    )


# =============================================================================
# Main workflow
# =============================================================================
def run(
    config_path: Path,
    modis_dir: Path,
    output_root: Path | None,
    cell_coverage_threshold: float,
    basin_coverage_threshold: float,
    overwrite: bool,
    max_files: int | None,
    skip_figures: bool,
) -> None:
    logger = _setup_logger()

    cfg = load_config(config_path)

    if output_root is None:
        output_root = cfg.output_dir / DEFAULT_SUBDIR

    settings = ValidationSettings(
        output_root=Path(output_root).expanduser(),
        modis_dir=Path(modis_dir).expanduser(),
        cell_coverage_threshold=float(cell_coverage_threshold),
        basin_coverage_threshold=float(basin_coverage_threshold),
        overwrite=bool(overwrite),
        max_files=max_files,
        skip_figures=bool(skip_figures),
    )

    if not 0 <= settings.cell_coverage_threshold <= 1:
        raise ValueError("cell_coverage_threshold must be in [0,1].")
    if not 0 <= settings.basin_coverage_threshold <= 1:
        raise ValueError("basin_coverage_threshold must be in [0,1].")
    if settings.max_files is not None and settings.max_files <= 0:
        raise ValueError("max_files must be a positive integer when provided.")
    if not settings.modis_dir.is_dir():
        raise FileNotFoundError(
            f"MODIS directory does not exist: {settings.modis_dir}"
        )

    settings.create_tree()

    logger.info("East River workflow version: %s", ERW_VERSION)
    logger.info("Validation script version: %s", SCRIPT_VERSION)
    logger.info("Config: %s", config_path)
    logger.info("MODIS directory: %s", settings.modis_dir)
    logger.info("Output root: %s", settings.output_root)

    if not settings.skip_figures:
        _configure_publication_matplotlib()

    # -------------------------------------------------------------------------
    # Static grid / watershed
    # -------------------------------------------------------------------------
    (
        noah_grid,
        weights,
        noah_dem,
        noah_lat,
        noah_lon,
        noah_diag,
    ) = _load_noah_static(cfg, settings, logger)

    # -------------------------------------------------------------------------
    # Inventory
    # -------------------------------------------------------------------------
    inventory, alignment = _build_modis_inventory(cfg, settings)

    inventory_path = settings.diagnostics_dir / "modis_file_inventory.csv"
    alignment_path = settings.diagnostics_dir / "modis_date_alignment.csv"
    inventory.to_csv(inventory_path, index=False)
    alignment.to_csv(alignment_path, index=False)

    logger.info(
        "MODIS files in analysis window: %d",
        int(np.count_nonzero(alignment["modis_available"])),
    )
    logger.info(
        "Missing MODIS dates in analysis window: %d",
        int(np.count_nonzero(~alignment["modis_available"])),
    )
    logger.info(
        "Dates selected for this run: %d",
        int(np.count_nonzero(alignment["selected_for_run"])),
    )

    # -------------------------------------------------------------------------
    # Noah hourly FSNO
    # -------------------------------------------------------------------------
    (
        noah_context,
        noah_ds,
        fsno_var,
        noah_times,
    ) = _open_noah_fsno(cfg, noah_grid.shape)

    # -------------------------------------------------------------------------
    # Process dates
    # -------------------------------------------------------------------------
    daily_rows: list[dict[str, Any]] = []
    qa_rows: list[dict[str, Any]] = []
    cell_acc = CellAccumulator.create(noah_grid.shape)

    wy_dates: dict[int, list[pd.Timestamp]] = {}
    wy_arrays: dict[int, dict[str, list[np.ndarray]]] = {}

    try:
        selected = alignment[alignment["selected_for_run"]].copy()

        for run_number, row in enumerate(selected.itertuples(index=False), start=1):
            date = pd.Timestamp(row.date).normalize()
            hdf_path = Path(row.path)

            logger.info(
                "[%d/%d] %s | %s",
                run_number,
                len(selected),
                date.strftime("%Y-%m-%d"),
                hdf_path.name,
            )

            daily_row, qa_row, arrays = _process_one_day(
                date,
                hdf_path,
                cfg,
                settings,
                noah_grid,
                weights,
                fsno_var,
                noah_times,
            )

            daily_rows.append(daily_row)
            qa_rows.append(qa_row)

            wy = int(daily_row["water_year"])
            wy_dates.setdefault(wy, []).append(date)
            wy_arrays.setdefault(
                wy,
                {
                    "modis_fsca": [],
                    "modis_valid_coverage": [],
                    "noah_fsno": [],
                    "delta_fsno_noah_minus_modis": [],
                    "paired_mask": [],
                },
            )

            for key, value in arrays.items():
                wy_arrays[wy][key].append(value)

            cell_acc.update(
                arrays["modis_fsca"],
                arrays["noah_fsno"],
                arrays["paired_mask"].astype(bool),
                use_day=bool(daily_row["primary_usable"]),
            )

    finally:
        noah_context.__exit__(None, None, None)

    # -------------------------------------------------------------------------
    # Save caches
    # -------------------------------------------------------------------------
    cache_paths: list[str] = []
    for wy in sorted(wy_dates):
        cache_path = _save_wy_cache(
            settings,
            wy,
            wy_dates[wy],
            wy_arrays[wy],
        )
        cache_paths.append(str(cache_path))

    # -------------------------------------------------------------------------
    # Tables
    # -------------------------------------------------------------------------
    daily = pd.DataFrame(daily_rows).sort_values("date").reset_index(drop=True)
    qa = pd.DataFrame(qa_rows).sort_values("date").reset_index(drop=True)

    daily_path = settings.tables_dir / "modis_fsno_watershed_daily.csv"
    qa_path = settings.diagnostics_dir / "modis_qa_summary.csv"

    daily.to_csv(daily_path, index=False)
    qa.to_csv(qa_path, index=False)

    cell_stats = cell_acc.to_frame(
        weights,
        noah_dem,
        noah_lat,
        noah_lon,
    )
    cell_path = settings.tables_dir / "modis_fsno_cell_statistics.csv"
    cell_stats.to_csv(cell_path, index=False)

    annual = _annual_summary(daily)
    annual_path = settings.tables_dir / "modis_fsno_annual_summary.csv"
    annual.to_csv(annual_path, index=False)

    # -------------------------------------------------------------------------
    # Figures
    # -------------------------------------------------------------------------
    figure_paths: list[str] = []

    if not settings.skip_figures and len(daily):
        years = sorted(
            int(v)
            for v in cfg.section("project")["water_years"]
            if int(v) in set(daily["water_year"].astype(int))
        )

        figure_paths.extend(
            str(path)
            for path in _plot_yearly_comparison(
                cfg, settings, daily, years
            )
        )
        figure_paths.extend(
            str(path)
            for path in _plot_yearly_difference(
                cfg, settings, daily, years
            )
        )
        figure_paths.extend(
            str(path)
            for path in _plot_scatter(
                cfg, settings, daily, years
            )
        )
        figure_paths.extend(
            str(path)
            for path in _plot_spatial_bias(
                cfg,
                settings,
                cell_stats,
                noah_lat,
                noah_lon,
                weights,
            )
        )

    # -------------------------------------------------------------------------
    # Manifest
    # -------------------------------------------------------------------------
    manifest = {
        "script_version": SCRIPT_VERSION,
        "east_river_workflow_version": ERW_VERSION,
        "config_path": str(Path(config_path).expanduser().resolve()),
        "modis_directory": str(settings.modis_dir.resolve()),
        "output_root": str(settings.output_root.resolve()),
        "analysis_start": str(cfg.analysis_start),
        "analysis_end": str(cfg.analysis_end),
        "water_years": [
            int(v) for v in cfg.section("project")["water_years"]
        ],
        "modis_product": "MOD10A1 V061 Terra",
        "modis_tile_required": MODIS_TILE,
        "modis_collection_required": MODIS_COLLECTION,
        "modis_observation_variable": MODIS_SNOW_SDS,
        "modis_basic_qa_variable": MODIS_BASIC_QA_SDS,
        "modis_algorithm_qa_variable": MODIS_ALGORITHM_QA_SDS,
        "modis_granule_pointer_variable": MODIS_GRANULE_SDS,
        "primary_basic_qa_values": list(PRIMARY_BASIC_QA),
        "binary_snow_definition": (
            "NDSI_Snow_Cover > 0 => snow; "
            "NDSI_Snow_Cover == 0 => snow-free"
        ),
        "algorithm_flags_primary_filter": False,
        "cell_minimum_valid_modis_source_coverage": (
            settings.cell_coverage_threshold
        ),
        "basin_minimum_paired_watershed_area_fraction": (
            settings.basin_coverage_threshold
        ),
        "paired_watershed_area_fraction_definition": (
            "Fraction of East River fractional Noah-cell weight retained after "
            "per-cell MODIS coverage screening; accepted cells count with their "
            "full watershed fraction, consistent with the v0.8.2 regridding convention."
        ),
        "modis_valid_source_watershed_fraction_definition": (
            "East River fractional-weighted mean of the continuous MODIS valid "
            "source-area coverage field before per-cell thresholding."
        ),
        "time_basis": "UTC",
        "time_match_method": "MODIS swath midpoint -> nearest Noah-MP hourly UTC state",
        "time_match_maximum_minutes": MAX_TIME_MATCH_MINUTES,
        "missing_model_swath_policy": "Exclude whole Noah cell if any contributing MODIS swath lacks finite FSNO",
        "mixed_granule_method": (
            "Within each Noah cell, hourly FSNO states are combined using "
            "the valid MODIS source-area fraction contributed by each granule."
        ),
        "bias_sign": "Noah-MP - MODIS",
        "plotting_standard": "East_River_Plotting_Requirements_JHM_AMS_v6",
        "plotting_font": "Nimbus Sans (required; no silent fallback)",
        "figure_exports": "PNG 600 dpi + vector PDF",
        "figure01_rendering": "paired MODIS/Noah observations as discrete markers; no interpolation across missing dates",
        "missing_modis_policy": "leave missing; no interpolation",
        "noah_fsno_source": str(
            Path(cfg.data["paths"]["noahmp_output_file"]).expanduser()
        ),
        "noah_geo_source": str(
            Path(cfg.data["paths"]["noahmp_geo_file"]).expanduser()
        ),
        "watershed_weights_source": str(
            cfg.output_dir / "cache" / "east_river_fraction_noahmp.npz"
        ),
        "noah_grid_diagnostics": noah_diag,
        "inventory_file": str(inventory_path),
        "date_alignment_file": str(alignment_path),
        "daily_table": str(daily_path),
        "qa_table": str(qa_path),
        "cell_statistics_table": str(cell_path),
        "annual_summary_table": str(annual_path),
        "cache_files": cache_paths,
        "figure_files": figure_paths,
    }

    manifest_path = (
        settings.diagnostics_dir
        / "modis_validation_manifest.json"
    )
    save_json(manifest, manifest_path)

    logger.info("Done.")
    logger.info("Daily table: %s", daily_path)
    logger.info("Annual table: %s", annual_path)
    logger.info("Cell statistics: %s", cell_path)
    logger.info("Manifest: %s", manifest_path)


# =============================================================================
# CLI
# =============================================================================
def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate East River Noah-MP FSNO against time-matched "
            "MODIS Terra MOD10A1 V061 snow-cover observations."
        )
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"East River v0.8.2 YAML config (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--modis-dir",
        type=Path,
        default=DEFAULT_MODIS_DIR,
        help=f"MOD10A1 HDF directory (default: {DEFAULT_MODIS_DIR})",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "Output root. Default: "
            "<cfg.output_dir>/modis_fsno_validation_v0.1.5"
        ),
    )
    parser.add_argument(
        "--cell-coverage-threshold",
        type=float,
        default=DEFAULT_CELL_COVERAGE_THRESHOLD,
        help=(
            "Minimum valid MODIS source-area fraction required within a "
            "Noah cell (default: 0.75)."
        ),
    )
    parser.add_argument(
        "--basin-coverage-threshold",
        type=float,
        default=DEFAULT_BASIN_COVERAGE_THRESHOLD,
        help=(
            "Minimum paired East River weighted-area fraction required for "
            "a primary basin comparison day (default: 0.75)."
        ),
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help=(
            "Optional smoke-test limit on available MODIS files. "
            "Omit for the full analysis."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Reserved for compatibility with companion workflows. "
            "v0.1 currently rewrites its own output tables/caches."
        ),
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Process/cache/tables only; do not generate figures.",
    )

    return parser


def main() -> None:
    args = _parser().parse_args()

    run(
        config_path=args.config,
        modis_dir=args.modis_dir,
        output_root=args.output_root,
        cell_coverage_threshold=args.cell_coverage_threshold,
        basin_coverage_threshold=args.basin_coverage_threshold,
        overwrite=args.overwrite,
        max_files=args.max_files,
        skip_figures=args.skip_figures,
    )


if __name__ == "__main__":
    main()
