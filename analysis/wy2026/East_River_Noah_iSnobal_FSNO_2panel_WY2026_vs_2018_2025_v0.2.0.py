#!/usr/bin/env python3
"""FSNO-conditioned ablation diagnostic, v0.2.0.

WY2026 is a line of annual within-bin means (panel A) or counts (panel B).
The historical dashed curve is the median of the eight annual bin statistics;
shading is their Q25-Q75 (NumPy linear empirical quantiles), not a confidence
interval and not pooled cell-day variability. Missing bias for an empty bin
stays missing; a zero event count is real. Baseline sample counts are exported.
No envelope is drawn when fewer than two years have a defined bias statistic.
Connecting bin centers is a visual guide, not a fitted continuous response.

Eligible events: Feb 1-Apr 30, both models' signed SWE losses >1 mm/day,
on the common Noah grid, FSNO from the previous daily 23:00 state. Counts and
within-year means use unweighted eligible cell-days, as in v0.8.2 H1/H2.
Leap days are retained, matching the existing calendar-window definition.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Allow the script to live beside / outside the repository without installation.
# -----------------------------------------------------------------------------
def _bootstrap_east_river_import() -> None:
    try:
        import east_river_workflow  # noqa: F401
        return
    except ModuleNotFoundError:
        pass

    here = Path(__file__).resolve()
    candidates: list[Path] = []
    for parent in [here.parent, *here.parents]:
        candidates.extend([
            parent / "src",
            parent / "East_River_Workflow_Code_v0.8.2" / "src",
            parent / "erw_code" / "East_River_Workflow_Code_v0.8.2" / "src",
        ])
    for candidate in candidates:
        if (candidate / "east_river_workflow").is_dir():
            sys.path.insert(0, str(candidate))
            return


_bootstrap_east_river_import()

from east_river_workflow.config import WorkflowConfig, load_config  # noqa: E402
from east_river_workflow.data_access import decode_noah_times, read_noah_grid  # noqa: E402
from east_river_workflow.grids import fractional_polygon_weights, load_weights, select_watershed  # noqa: E402
from east_river_publication import configure_matplotlib, save_figure  # noqa: E402
from east_river_workflow.processing import load_spatial_cache, spatial_cache_path  # noqa: E402
from east_river_workflow.utils import open_dataset_robust, water_year_bounds  # noqa: E402
from east_river_workflow.constants import MODEL_COLORS  # noqa: E402


SCRIPT_VERSION = "0.2.0"
DEFAULT_SUBDIR = "fsno_conditioned_curves_v0.2"
DEFAULT_FOCUS_START = "02-01"
DEFAULT_FOCUS_END = "04-30"
DEFAULT_EVENT_THRESHOLD_MM_DAY = 1.0
DEFAULT_FSNO_VARIABLE = "FSNO"
DEFAULT_WY_START = 2018
DEFAULT_WY_END = 2026
FSNO_FULL_TOLERANCE = 1.0e-6
BIN_EDGES = np.linspace(0.0, 1.0, 11)  # 10 bins
BIN_LABELS = [f"{BIN_EDGES[i]:.1f}-{BIN_EDGES[i+1]:.1f}" for i in range(10)]


@dataclass(frozen=True)
class FigureSettings:
    output_root: Path
    focus_start_month_day: str = DEFAULT_FOCUS_START
    focus_end_month_day: str = DEFAULT_FOCUS_END
    event_threshold_mm_day: float = DEFAULT_EVENT_THRESHOLD_MM_DAY
    fsno_variable: str = DEFAULT_FSNO_VARIABLE
    overwrite: bool = False

    @property
    def figures_dir(self) -> Path:
        return self.output_root / "figures"

    @property
    def tables_dir(self) -> Path:
        return self.output_root / "tables"

    @property
    def diagnostics_dir(self) -> Path:
        return self.output_root / "diagnostics"

    def create_tree(self) -> None:
        for path in [self.figures_dir, self.tables_dir, self.diagnostics_dir]:
            path.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _season_date(wy: int, month_day: str) -> pd.Timestamp:
    month, day = [int(x) for x in month_day.split("-")]
    year = wy - 1 if month >= 10 else wy
    return pd.Timestamp(year=year, month=month, day=day)


def _focus_bounds(wy: int, settings: FigureSettings, cfg: WorkflowConfig) -> tuple[pd.Timestamp, pd.Timestamp]:
    start = _season_date(wy, settings.focus_start_month_day)
    end = _season_date(wy, settings.focus_end_month_day)
    wy_start, wy_end = water_year_bounds(int(wy))
    start = max(start, wy_start, cfg.analysis_start.normalize())
    end = min(end, wy_end, cfg.analysis_end.normalize())
    return start, end


def _focus_time_mask(dates: pd.DatetimeIndex, wy: int, settings: FigureSettings, cfg: WorkflowConfig) -> np.ndarray:
    start, end = _focus_bounds(wy, settings, cfg)
    return np.asarray((dates >= start) & (dates <= end), dtype=bool)


def _signed_ablation(swe_mm: np.ndarray, dates: pd.DatetimeIndex) -> np.ndarray:
    """Return t-1 minus t SWE; positive values indicate ablation."""
    swe = np.asarray(swe_mm, dtype=float)
    out = np.full(swe.shape, np.nan, dtype=float)
    if swe.shape[0] < 2:
        return out
    dt_days = np.asarray((dates[1:] - dates[:-1]).days, dtype=int)
    consecutive = dt_days == 1
    diff = swe[:-1] - swe[1:]
    valid = np.isfinite(swe[:-1]) & np.isfinite(swe[1:])
    valid &= consecutive.reshape((-1,) + (1,) * (swe.ndim - 1))
    out[1:] = np.where(valid, diff, np.nan)
    return out


def _shift_start_state(values: np.ndarray, dates: pd.DatetimeIndex) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    if arr.shape[0] < 2:
        return out
    dt_days = np.asarray((dates[1:] - dates[:-1]).days, dtype=int)
    consecutive = dt_days == 1
    shape = (-1,) + (1,) * (arr.ndim - 1)
    valid_time = consecutive.reshape(shape)
    out[1:] = np.where(valid_time, arr[:-1], np.nan)
    return out


def _load_noah_grid_and_weights(cfg: WorkflowConfig, logger: logging.Logger) -> tuple[Any, np.ndarray]:
    noah_grid, _noah_dem, _noah_lat, _noah_lon, _noah_diag = read_noah_grid(
        cfg.data["paths"]["noahmp_geo_file"], cfg.section("noahmp")
    )
    main_weight_path = cfg.output_dir / "cache" / "east_river_fraction_noahmp.npz"
    if main_weight_path.exists():
        weights = load_weights(main_weight_path, noah_grid)
        logger.info("Loaded existing Noah watershed weights: %s", main_weight_path)
    else:
        logger.warning("Main Noah watershed-weight cache is absent; rebuilding exact fractional weights.")
        watershed = select_watershed(
            cfg.data["paths"]["watershed_shapefile"],
            cfg.section("watershed")["name_field"],
            cfg.section("watershed")["name_contains"],
        )
        weights = fractional_polygon_weights(
            noah_grid,
            watershed,
            chunk_size=int(cfg.section("watershed")["fractional_chunk_size"]),
        )
    return noah_grid, np.asarray(weights, dtype=float)


def _extract_fsno_for_years(
    cfg: WorkflowConfig,
    settings: FigureSettings,
    years: list[int],
    noah_cache_dates: dict[int, pd.DatetimeIndex],
    noah_grid_shape: tuple[int, int],
    logger: logging.Logger,
) -> dict[int, dict[str, Any]]:
    source = Path(cfg.data["paths"]["noahmp_output_file"]).expanduser()
    logger.info("Extracting native Noah-MP %s from %s", settings.fsno_variable, source)
    result: dict[int, dict[str, Any]] = {}

    with open_dataset_robust(source, decode_times=False) as ds:
        if settings.fsno_variable not in ds:
            raise KeyError(
                f"Required Noah-MP FSNO variable {settings.fsno_variable!r} is absent from {source}."
            )
        fsno_var = ds[settings.fsno_variable]
        if "Time" not in fsno_var.dims:
            raise ValueError(f"{settings.fsno_variable!r} must have a Time dimension; dims={fsno_var.dims}.")

        times = decode_noah_times(ds, cfg.section("noahmp"))
        frame = pd.DataFrame({"index": np.arange(len(times), dtype=int), "time": times})
        frame["date"] = frame["time"].dt.normalize()
        frame["hour"] = frame["time"].dt.hour
        state_hour = int(cfg.section("noahmp")["daily_state_hour"])
        state = frame[frame["hour"].eq(state_hour)].copy()
        duplicate = state["date"].duplicated(keep=False)
        if duplicate.any():
            dup_dates = state.loc[duplicate, "date"].dt.strftime("%Y-%m-%d").unique().tolist()[:10]
            raise ValueError(f"More than one Noah state record at hour {state_hour} for dates such as {dup_dates}")
        state_index = state.set_index("date")["index"]

        for wy in years:
            target_dates = pd.DatetimeIndex(noah_cache_dates[wy]).normalize()
            missing = target_dates.difference(state_index.index)
            if len(missing):
                raise ValueError(
                    f"No {state_hour}:00 FSNO state record for {len(missing)} Noah cache dates in WY{wy}; "
                    f"first missing={missing[0]}."
                )
            indices = state_index.reindex(target_dates).to_numpy(dtype=int)
            fsno = np.asarray(fsno_var.isel(Time=indices).values, dtype=float)
            if fsno.shape != (len(target_dates), *noah_grid_shape):
                raise ValueError(
                    f"FSNO shape mismatch for WY{wy}: {fsno.shape}, expected {(len(target_dates), *noah_grid_shape)}"
                )
            finite = fsno[np.isfinite(fsno)]
            if finite.size and (float(np.min(finite)) < -1.0e-6 or float(np.max(finite)) > 1.0 + 1.0e-6):
                raise ValueError(
                    f"FSNO outside [0,1] in WY{wy}: min={float(np.min(finite))}, max={float(np.max(finite))}"
                )
            fsno = np.where(np.isfinite(fsno), np.clip(fsno, 0.0, 1.0), np.nan)
            result[wy] = {"dates": target_dates, "fsno": fsno}

    return result


def _align_main_caches(ic: dict[str, Any], nc: dict[str, Any]) -> tuple[pd.DatetimeIndex, np.ndarray, np.ndarray]:
    idates = pd.DatetimeIndex(ic["dates"]).normalize()
    ndates = pd.DatetimeIndex(nc["dates"]).normalize()
    common = idates.intersection(ndates).sort_values()
    if not len(common):
        raise ValueError("iSnobal and Noah caches contain no common dates.")
    ii = idates.get_indexer(common)
    ni = ndates.get_indexer(common)
    return common, ii, ni


def _fsno_bin_codes_10(values: np.ndarray) -> np.ndarray:
    """Map FSNO to integer bin codes 0..9 with 1.0 included in the last bin."""
    arr = np.asarray(values, dtype=float)
    code = np.full(arr.shape, -1, dtype=np.int8)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return code

    clipped = np.clip(arr[finite], 0.0, 1.0)
    # digitize with right=False gives [edge_i, edge_{i+1}) except last edge.
    # Multiplying by 10 and flooring is robust here; then clamp 1.0 into bin 9.
    bins = np.floor(clipped * 10.0).astype(int)
    bins = np.clip(bins, 0, 9)
    # Protect exact 1.0 and near-1.0 values inside the final bin.
    bins[np.isclose(clipped, 1.0, rtol=0.0, atol=FSNO_FULL_TOLERANCE)] = 9
    code[finite] = bins.astype(np.int8)
    return code


def _finite_quantile(x: np.ndarray, q: float) -> float:
    arr = np.asarray(x, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.quantile(arr, q)) if arr.size else np.nan


def _save_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)


# -----------------------------------------------------------------------------
# Core analysis
# -----------------------------------------------------------------------------
def _collect_year_bin_statistics(
    cfg: WorkflowConfig,
    settings: FigureSettings,
    years: list[int],
    weights: np.ndarray,
    fsno_info: dict[int, dict[str, Any]],
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-year, per-bin summary plus a cell-day table for bias diagnostics."""
    watershed_mask_2d = np.asarray(weights, dtype=float) > 0.0

    year_rows: list[dict[str, Any]] = []
    cell_rows: list[dict[str, Any]] = []

    for wy in years:
        ipath = spatial_cache_path(cfg, "isnobal", wy)
        npath = spatial_cache_path(cfg, "noahmp", wy)
        if not ipath.exists() or not npath.exists():
            missing = [str(p) for p in [ipath, npath] if not p.exists()]
            raise FileNotFoundError(
                "Required v0.8.2 spatial cache(s) are missing: " + ", ".join(missing)
            )
        logger.info("Loading WY%s caches", wy)
        ic = load_spatial_cache(ipath)
        nc = load_spatial_cache(npath)

        for source_name, cache in [("iSnobal", ic), ("Noah-MP", nc)]:
            for key in ["dates", "swe_mm"]:
                if key not in cache:
                    raise KeyError(f"{source_name} WY{wy} cache is missing required key {key!r}.")

        dates, ii, ni = _align_main_caches(ic, nc)
        start, end = _focus_bounds(wy, settings, cfg)
        required_dates = pd.date_range(start - pd.Timedelta(days=1), end, freq="D")
        if dates.has_duplicates or not required_dates.isin(dates).all():
            raise ValueError(f"WY{wy}: incomplete/duplicate model dates; require every focus day and preceding state.")
        fdates = pd.DatetimeIndex(fsno_info[wy]["dates"]).normalize()
        fi = fdates.get_indexer(dates)
        if np.any(fi < 0):
            raise ValueError(f"FSNO cache does not cover all aligned model dates for WY{wy}.")

        swe_i = np.asarray(ic["swe_mm"][ii], dtype=float)
        swe_n = np.asarray(nc["swe_mm"][ni], dtype=float)
        fsno = np.asarray(fsno_info[wy]["fsno"], dtype=float)[fi]

        ablation_i = _signed_ablation(swe_i, dates)
        ablation_n = _signed_ablation(swe_n, dates)
        bias = ablation_n - ablation_i
        fsno_start = _shift_start_state(fsno, dates)
        fsno_bin = _fsno_bin_codes_10(fsno_start)

        focus_mask_t = _focus_time_mask(dates, wy, settings, cfg)
        focus_mask = focus_mask_t[:, None, None]
        spatial_mask = np.broadcast_to(watershed_mask_2d, ablation_i.shape)
        event_mask = (
            focus_mask
            & spatial_mask
            & np.isfinite(ablation_i)
            & np.isfinite(ablation_n)
            & (ablation_i > settings.event_threshold_mm_day)
            & (ablation_n > settings.event_threshold_mm_day)
            & np.isfinite(fsno_start)
            & (fsno_bin >= 0)
        )

        for bin_code, label in enumerate(BIN_LABELS):
            keep = event_mask & (fsno_bin == bin_code)
            e = np.asarray(bias[keep], dtype=float)
            count = int(e.size)
            row = {
                "water_year": int(wy),
                "bin_code": int(bin_code),
                "bin_label": label,
                "event_count": count,
                "mean_bias_mm_day": float(np.mean(e)) if count else np.nan,
                "q25_bias_mm_day": _finite_quantile(e, 0.25),
                "median_bias_mm_day": _finite_quantile(e, 0.50),
                "q75_bias_mm_day": _finite_quantile(e, 0.75),
            }
            year_rows.append(row)

            # Save only the bias values needed for optional validation / direct checking.
            for value in e:
                cell_rows.append({
                    "water_year": int(wy),
                    "bin_code": int(bin_code),
                    "bin_label": label,
                    "bias_mm_day": float(value),
                })

    return pd.DataFrame(year_rows), pd.DataFrame(cell_rows)


def _build_comparison_table(year_bin: pd.DataFrame) -> pd.DataFrame:
    """Create one summary row per FSNO bin for the requested two-panel figure."""
    counts = year_bin["event_count"].to_numpy(dtype=float)
    if not np.isfinite(counts).all() or np.any(counts < 0) or np.any(counts != np.floor(counts)):
        raise ValueError("Annual event counts must be finite nonnegative integers; missing is not zero.")
    means = year_bin["mean_bias_mm_day"].to_numpy(dtype=float)
    if np.any((counts > 0) & ~np.isfinite(means)) or np.any((counts == 0) & np.isfinite(means)):
        raise ValueError("A bin mean is defined exactly when its event count is positive.")
    rows: list[dict[str, Any]] = []

    baseline = year_bin[year_bin["water_year"].between(2018, 2025)].copy()
    wy2026 = year_bin[year_bin["water_year"] == 2026].copy()

    for bin_code, label in enumerate(BIN_LABELS):
        b = baseline[baseline["bin_code"] == bin_code].sort_values("water_year")
        y26 = wy2026[wy2026["bin_code"] == bin_code]
        if len(y26) != 1:
            raise ValueError(f"Expected exactly one WY2026 summary row for bin {label}, found {len(y26)}.")
        y26_row = y26.iloc[0]

        if b["water_year"].tolist() != list(range(2018, 2026)):
            raise ValueError(f"Expected exactly one summary for each historical year in bin {label}.")
        baseline_means = b["mean_bias_mm_day"].to_numpy(dtype=float)
        baseline_counts = b["event_count"].to_numpy(dtype=float)

        rows.append({
            "bin_code": int(bin_code),
            "bin_label": label,
            "wy2018_2025_median_of_yearly_mean_bias_mm_day": _finite_quantile(baseline_means, 0.5),
            "wy2018_2025_median_annual_event_count": _finite_quantile(baseline_counts, 0.5),
            # Left panel: bias
            "wy2026_mean_bias_mm_day": float(y26_row["mean_bias_mm_day"]),
            "wy2026_q25_bias_mm_day": float(y26_row["q25_bias_mm_day"]),
            "wy2026_q75_bias_mm_day": float(y26_row["q75_bias_mm_day"]),
            "wy2018_2025_mean_of_yearly_mean_bias_mm_day": float(np.nanmean(baseline_means)) if np.isfinite(baseline_means).any() else np.nan,
            "wy2018_2025_q25_of_yearly_mean_bias_mm_day": _finite_quantile(baseline_means, 0.25),
            "wy2018_2025_q75_of_yearly_mean_bias_mm_day": _finite_quantile(baseline_means, 0.75),
            # Right panel: event count
            "wy2026_event_count": int(y26_row["event_count"]),
            "wy2018_2025_mean_annual_event_count": float(np.nanmean(baseline_counts)) if len(baseline_counts) else np.nan,
            "wy2018_2025_q25_annual_event_count": _finite_quantile(baseline_counts, 0.25),
            "wy2018_2025_q75_annual_event_count": _finite_quantile(baseline_counts, 0.75),
            # Diagnostics
            "baseline_years_with_bias": int(np.count_nonzero(np.isfinite(baseline_means))),
            "baseline_years_with_counts": int(np.count_nonzero(np.isfinite(baseline_counts))),
        })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def _panel_curves(ax, x, focus, median, q25, q75, n):
    # The contrast is between years, not between models: use neutral gray for
    # history, never iSnobal blue to label historical Noah-conditioned statistics.
    valid = (np.asarray(n) >= 2) & np.isfinite(q25) & np.isfinite(q75)
    ax.fill_between(x, q25, q75, where=valid, color="0.65", alpha=0.35,
                    linewidth=0, label="WY 2018–2025 IQR")
    ax.plot(x, median, color="0.35", linestyle="--", marker="s", linewidth=2.5,
            label="WY 2018–2025 median")
    ax.plot(x, focus, color=MODEL_COLORS["Noah-MP"], marker="o", linewidth=2.5,
            label="WY 2026")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_xticks(np.arange(0, 1.01, 0.2))


def _plot_two_panel_figure(cfg: WorkflowConfig, summary: pd.DataFrame, out_path: Path) -> Path:
    configure_matplotlib(cfg)
    summary = summary.sort_values("bin_code")
    if summary["bin_code"].tolist() != list(range(10)):
        raise ValueError("Expected exactly ten ordered FSNO bins.")
    x = (BIN_EDGES[:-1] + BIN_EDGES[1:]) / 2
    fig, axes = plt.subplots(1, 2, figsize=(18, 7.5))
    _panel_curves(axes[0], x, summary.wy2026_mean_bias_mm_day,
        summary.wy2018_2025_median_of_yearly_mean_bias_mm_day,
        summary.wy2018_2025_q25_of_yearly_mean_bias_mm_day,
        summary.wy2018_2025_q75_of_yearly_mean_bias_mm_day, summary.baseline_years_with_bias)
    _panel_curves(axes[1], x, summary.wy2026_event_count,
        summary.wy2018_2025_median_annual_event_count,
        summary.wy2018_2025_q25_annual_event_count,
        summary.wy2018_2025_q75_annual_event_count, summary.baseline_years_with_counts)
    axes[0].axhline(0, color="0.6", linewidth=0.8)
    axes[0].set_ylabel("Mean ablation bias (mm/day)\nNoah-MP − iSnobal")
    axes[1].set_ylabel("Eligible cell-days per water year")
    axes[1].set_ylim(bottom=0)
    axes[0].set_title("Ablation bias", pad=8)
    axes[1].set_title("Event count", pad=8)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend([handles[i] for i in [2,1,0]], [labels[i] for i in [2,1,0]],
                   loc="lower left", borderpad=0.5, handlelength=2.2)
    fig.supxlabel("Noah-MP FSNO at start of daily ablation interval", y=0.04)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.92, bottom=0.16, wspace=0.28)
    save_figure(fig, out_path)
    return out_path.with_suffix(".png")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to East River workflow config YAML.")
    parser.add_argument("--output-subdir", default=DEFAULT_SUBDIR, help="Output subdirectory under config output_dir.")
    parser.add_argument("--focus-start", default=DEFAULT_FOCUS_START, help="Focus-window start as MM-DD (default: 02-01).")
    parser.add_argument("--focus-end", default=DEFAULT_FOCUS_END, help="Focus-window end as MM-DD (default: 04-30).")
    parser.add_argument("--event-threshold", type=float, default=DEFAULT_EVENT_THRESHOLD_MM_DAY,
                        help="Ablation threshold in mm/day applied to BOTH models with an AND rule (default: 1.0).")
    parser.add_argument("--fsno-variable", default=DEFAULT_FSNO_VARIABLE, help="FSNO variable name in LDASOUT (default: FSNO).")
    parser.add_argument("--wy-start", type=int, default=DEFAULT_WY_START, help="First water year to include (default: 2018).")
    parser.add_argument("--wy-end", type=int, default=DEFAULT_WY_END, help="Last water year to include (default: 2026).")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger("fsno_2panel")

    cfg = load_config(args.config)
    years = list(range(int(args.wy_start), int(args.wy_end) + 1))
    if years != list(range(2018, 2027)):
        raise ValueError("This diagnostic requires exactly WY2018-WY2026.")
    if 2018 not in years or 2025 not in years:
        raise ValueError("This script expects WY2018-WY2025 to exist for the baseline comparison.")

    settings = FigureSettings(
        output_root=cfg.output_dir / args.output_subdir,
        focus_start_month_day=args.focus_start,
        focus_end_month_day=args.focus_end,
        event_threshold_mm_day=float(args.event_threshold),
        fsno_variable=args.fsno_variable,
        overwrite=bool(args.overwrite),
    )
    settings.create_tree()

    noah_grid, weights = _load_noah_grid_and_weights(cfg, logger)

    # Get Noah cache dates first so FSNO can be aligned exactly to those dates.
    noah_cache_dates: dict[int, pd.DatetimeIndex] = {}
    for wy in years:
        npath = spatial_cache_path(cfg, "noahmp", wy)
        if not npath.exists():
            raise FileNotFoundError(f"Required Noah-MP cache missing: {npath}")
        nc = load_spatial_cache(npath)
        if "dates" not in nc:
            raise KeyError(f"Noah-MP WY{wy} cache missing required key 'dates'.")
        noah_cache_dates[wy] = pd.DatetimeIndex(nc["dates"])

    fsno_info = _extract_fsno_for_years(cfg, settings, years, noah_cache_dates, tuple(noah_grid.shape), logger)
    year_bin, cell_bias = _collect_year_bin_statistics(cfg, settings, years, weights, fsno_info, logger)
    summary = _build_comparison_table(year_bin)

    summary_csv = settings.tables_dir / "fsno_two_panel_wy2026_vs_wy2018_2025_summary.csv"
    yearbin_csv = settings.tables_dir / "fsno_two_panel_per_year_bin_summary.csv"
    cell_csv = settings.tables_dir / "fsno_two_panel_bias_values_long.csv"
    fig_path = settings.figures_dir / "fsno_two_panel_wy2026_vs_wy2018_2025.png"
    manifest_path = settings.diagnostics_dir / "fsno_two_panel_manifest.json"

    _save_csv(summary, summary_csv)
    _save_csv(year_bin, yearbin_csv)
    _save_csv(cell_bias, cell_csv)
    _plot_two_panel_figure(cfg, summary, fig_path)

    manifest = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "config": str(Path(args.config).resolve()),
        "output_root": str(settings.output_root),
        "years": years,
        "focus_window": {"start": settings.focus_start_month_day, "end": settings.focus_end_month_day},
        "event_rule": "keep only cell-days where A_I > threshold AND A_N > threshold",
        "event_threshold_mm_day": float(settings.event_threshold_mm_day),
        "fsno_binning": {
            "bin_edges": BIN_EDGES.tolist(),
            "labels": BIN_LABELS,
            "include_fsno_1_in_last_bin": True,
            "use_fsno_start_of_interval": True,
        },
        "figure": str(fig_path),
        "figure_pdf": str(fig_path.with_suffix(".pdf")),
        "historical_curve": "Median of annual bin statistics; equal weight per defined year",
        "envelope": "Empirical Q25-Q75 of annual bin statistics; not a confidence interval; n>=2",
        "event_weighting": "Unweighted eligible watershed cell-days, not fractional-area weighting",
        "leap_days": "Retained; annual Feb-Apr exposure is 89 or 90 days",
        "plotting_standard": "East_River_Plotting_Requirements_JHM_AMS_v6",
        "tables": [str(summary_csv), str(yearbin_csv), str(cell_csv)],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))

    print(json.dumps({
        "status": "complete",
        "figure": str(fig_path),
        "figure_pdf": str(fig_path.with_suffix(".pdf")),
        "historical_curve": "Median of annual bin statistics; equal weight per defined year",
        "envelope": "Empirical Q25-Q75 of annual bin statistics; not a confidence interval; n>=2",
        "event_weighting": "Unweighted eligible watershed cell-days, not fractional-area weighting",
        "leap_days": "Retained; annual Feb-Apr exposure is 89 or 90 days",
        "plotting_standard": "East_River_Plotting_Requirements_JHM_AMS_v6",
        "summary_csv": str(summary_csv),
        "year_bin_csv": str(yearbin_csv),
        "manifest": str(manifest_path),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
