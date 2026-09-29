#!/usr/bin/env python3
"""East River WY2021/WY2026 diagnostics with estimated daily net longwave.

Version 0.3.1. No model rerun. Feb 1-Apr 30, midnight-to-midnight UTC.
Rows: (1) FSNO and delta SWE; (2) delta albedo and absorbed SW;
(3) Noah-MP meteorological forcing temperature and delta estimated net LW (Noah-MP minus iSnobal).

Longwave uses common snow-present cells, not the whole-basin support used
by the other rows. iSnobal = (23*Rn[D,23]+Rn[D+1,00])/24 minus solar integrated
from 25 hourly endpoints. Noah = daily integral of -FIRA hourly samples.
Gray shading covers ONLY solar quadrature uncertainty conditional on
uninterrupted snow; it is not total uncertainty or a confidence interval.

Snow screening is at available checkpoints, not proof of continuous snow.
Noah-MP canopy/partial-snow radiation and iSnobal snow-surface radiation remain
physically different diagnostics. Reported common area may change each day.
Nimbus Sans; PNG 600 dpi and vector PDF. See README.md for assumptions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Allow the script to live near the repository root, examples/, or elsewhere.
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
from east_river_workflow.constants import MODEL_COLORS  # noqa: E402
from east_river_publication import configure_matplotlib, save_figure  # noqa: E402
from east_river_daily_longwave import build_daily_longwave  # noqa: E402


SCRIPT_VERSION = "0.3.1"
DEFAULT_OUTPUT_SUBDIR = "wy2021_wy2026_comparison_v0.3"
DEFAULT_SEPARATION_SUBDIR = "separation_diagnostic_v0.1"
DEFAULT_FORCING_SUBDIR = "ldasin_wy2026_forcing_v0.1.1"
TARGET_WATER_YEARS = (2021, 2026)
FOCUS_START = "02-01"
FOCUS_END = "04-30"
SECONDS_PER_DAY = 86400.0
MJ_TO_WM2_DAY = 1.0e6 / SECONDS_PER_DAY
DEFAULT_PRECIP_EVENT_THRESHOLD_MM_DAY = 1.0


@dataclass(frozen=True)
class Settings:
    output_root: Path
    separation_dir: Path
    forcing_dir: Path
    focus_start_month_day: str
    focus_end_month_day: str
    precip_event_threshold_mm_day: float

    @property
    def tables_dir(self) -> Path:
        return self.output_root / "tables"

    @property
    def figures_dir(self) -> Path:
        return self.output_root / "figures"

    def create_tree(self) -> None:
        for path in [self.output_root, self.tables_dir, self.figures_dir]:
            path.mkdir(parents=True, exist_ok=True)


def _season_date(wy: int, month_day: str) -> pd.Timestamp:
    month, day = [int(x) for x in month_day.split("-")]
    year = wy - 1 if month >= 10 else wy
    return pd.Timestamp(year=year, month=month, day=day)


def _focus_bounds(wy: int, settings: Settings, cfg: WorkflowConfig) -> tuple[pd.Timestamp, pd.Timestamp]:
    start = _season_date(wy, settings.focus_start_month_day)
    end = _season_date(wy, settings.focus_end_month_day)
    start = max(start, cfg.analysis_start.normalize())
    end = min(end, cfg.analysis_end.normalize())
    return start, end


def _save_figure(fig, cfg, settings, stem) -> Path:
    png, _ = save_figure(fig, settings.figures_dir / stem,
                        close=bool(cfg.section("plotting").get("close_after_save", True)))
    return png


def _read_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if "date" in frame:
        frame["date"] = pd.to_datetime(frame["date"])
    return frame


def _load_requested_daily_products(
    cfg: WorkflowConfig,
    settings: Settings,
    years: tuple[int, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sep_path = settings.separation_dir / "tables" / "separation_watershed_daily_full_period.csv"
    forcing_path = settings.forcing_dir / "tables" / "ldasin_watershed_daily_oct1_jun30.csv"
    if not sep_path.exists():
        raise FileNotFoundError(
            f"Required separation table is missing: {sep_path}\n"
            "Run East_River_Noah_iSnobal_Separation_Diagnostic_v0.1.py first."
        )
    if not forcing_path.exists():
        raise FileNotFoundError(
            f"Required forcing table is missing: {forcing_path}\n"
            "Run East_River_LDASIN_WY2026_Forcing_Climatology_v0.1.py first."
        )

    sep = _read_table(sep_path)
    forcing = _read_table(forcing_path)

    keep_sep = [
        "date", "water_year", "fsno_n_watershed_mean", "delta_swe_n_minus_i_mm",
        "delta_effective_albedo_n_minus_i", "delta_absorbed_sw_n_minus_i_mj_m2_day",
    ]
    keep_forcing = ["date", "water_year", "precip_mm_day", "temperature_mean_c",
                    "hour_count"]

    missing_sep = sorted(set(keep_sep) - set(sep.columns))
    missing_forcing = sorted(set(keep_forcing) - set(forcing.columns))
    if missing_sep:
        raise KeyError(f"Separation table is missing required columns: {missing_sep}")
    if missing_forcing:
        raise KeyError(f"Forcing table is missing required columns: {missing_forcing}")

    sep = sep[keep_sep].copy()
    forcing = forcing[keep_forcing].copy()

    focus_parts_sep: list[pd.DataFrame] = []
    focus_parts_forcing: list[pd.DataFrame] = []
    for wy in years:
        start, end = _focus_bounds(wy, settings, cfg)
        focus_parts_sep.append(sep[(sep["water_year"] == wy) & (sep["date"] >= start) & (sep["date"] <= end)].copy())
        focus_parts_forcing.append(
            forcing[(forcing["water_year"] == wy) & (forcing["date"] >= start) & (forcing["date"] <= end)].copy()
        )

    sep_focus = pd.concat(focus_parts_sep, ignore_index=True)
    forcing_focus = pd.concat(focus_parts_forcing, ignore_index=True)

    for label, frame in [("separation", sep_focus), ("forcing", forcing_focus)]:
        if frame.duplicated(["date", "water_year"]).any():
            raise ValueError(f"Duplicate {label} date/year keys.")
        expected = pd.DatetimeIndex([])
        for wy in years:
            start, end = _focus_bounds(wy, settings, cfg)
            expected = expected.append(pd.date_range(start, end, freq="D"))
        actual = pd.DatetimeIndex(frame["date"]).sort_values()
        if not actual.equals(expected.sort_values()):
            raise ValueError(f"{label} does not contain every requested comparison date.")
        if not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
            raise ValueError(f"Non-finite {label} comparison values.")

    if not forcing_focus["hour_count"].eq(24).all():
        raise ValueError("Comparison forcing must contain 24-hour daily means.")
    if not sep_focus["fsno_n_watershed_mean"].between(0, 1).all():
        raise ValueError("FSNO must lie in [0,1].")
    sep_focus["delta_net_sw_n_minus_i_w_m2"] = sep_focus["delta_absorbed_sw_n_minus_i_mj_m2_day"] * MJ_TO_WM2_DAY
    forcing_focus["precip_event"] = forcing_focus["precip_mm_day"] >= settings.precip_event_threshold_mm_day
    return sep_focus, forcing_focus


def _merge_all_products(sep, forcing) -> pd.DataFrame:
    keys = ["date", "water_year"]
    reference = pd.MultiIndex.from_frame(sep[keys])
    for name, frame in [("separation", sep), ("forcing", forcing)]:
        index = pd.MultiIndex.from_frame(frame[keys])
        if index.has_duplicates or len(index) != len(reference) or len(index.difference(reference)):
            raise ValueError(f"{name} date keys differ; no silent inner-join date loss is allowed.")
    return sep.merge(forcing, on=keys, how="left", validate="one_to_one").sort_values(keys).reset_index(drop=True)


def _symmetric_limits(values: np.ndarray, pad_fraction: float = 0.08, minimum_half_range: float = 1.0) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        half = minimum_half_range
    else:
        half = max(float(np.nanmax(np.abs(arr))), minimum_half_range)
    half *= (1.0 + pad_fraction)
    return -half, half


def _padded_limits(values: np.ndarray, pad_fraction: float = 0.08) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0, 1.0
    lo = float(np.nanmin(arr))
    hi = float(np.nanmax(arr))
    if math.isclose(lo, hi, rel_tol=0.0, abs_tol=1e-12):
        span = 1.0 if math.isclose(lo, 0.0, rel_tol=0.0, abs_tol=1e-12) else abs(lo) * 0.2
        return lo - span, hi + span
    pad = (hi - lo) * pad_fraction
    return lo - pad, hi + pad


def _format_date_axis(ax: plt.Axes, dates: pd.Series) -> None:
    ax.set_xlim(dates.min(), dates.max())
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.spines[["top"]].set_visible(False)


def _add_precip_short_ticks(ax: plt.Axes, event_dates: pd.Series) -> None:
    if len(event_dates) == 0:
        return
    ymin, ymax = ax.get_ylim()
    tick_top = ymin + 0.08 * (ymax - ymin)
    ax.vlines(event_dates, ymin=ymin, ymax=tick_top, linewidth=1.0, color=MODEL_COLORS["Noah-MP"])
    ax.set_ylim(ymin, ymax)


def make_figure(frame: pd.DataFrame, cfg: WorkflowConfig, settings: Settings) -> plt.Figure:
    configure_matplotlib(cfg)
    plt.rcParams.update({
        "font.size": 18,
        "axes.labelsize": 18,
        "axes.titlesize": 20,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
    })

    c_noah = MODEL_COLORS.get("Noah-MP", "tab:blue")
    # Noah-MP states/forcing use Noah red; model differences use neutral gray.
    color_primary = c_noah
    color_secondary = "0.25"  # Model differences are not iSnobal-only quantities.

    fig, axes = plt.subplots(3, 2, figsize=(21, 16.5), squeeze=False)
    fig.subplots_adjust(top=0.95, bottom=0.15, left=0.07, right=0.93, hspace=0.30, wspace=0.30)

    yr_data = {wy: frame[frame["water_year"] == wy].copy() for wy in TARGET_WATER_YEARS}
    for wy, sub in yr_data.items():
        if sub.empty:
            raise ValueError(f"No merged rows available for WY{wy} in requested focus window.")

    fsno_lim = (0.0, 1.0)
    dswe_lim = _symmetric_limits(frame["delta_swe_n_minus_i_mm"].to_numpy(), minimum_half_range=10.0)
    dalb_lim = _symmetric_limits(frame["delta_effective_albedo_n_minus_i"].to_numpy(), minimum_half_range=0.02)
    dnsw_lim = _symmetric_limits(frame["delta_net_sw_n_minus_i_w_m2"].to_numpy(), minimum_half_range=5.0)
    temp_lim = _padded_limits(frame["temperature_mean_c"].to_numpy())
    inlw_lim = _symmetric_limits(frame[["delta_lw_solar_lower_w_m2", "delta_lw_solar_upper_w_m2"]].to_numpy().ravel(), minimum_half_range=5.0)

    row_specs = [
        {
            "left_col": "fsno_n_watershed_mean",
            "right_col": "delta_swe_n_minus_i_mm",
            "left_label": "Noah-MP FSNO",
            "right_label": "Δ SWE (mm)",
            "left_lim": fsno_lim,
            "right_lim": dswe_lim,
            "legend_labels": ("FSNO", "Δ SWE"),
            "left_is_delta": False,
            "right_is_delta": True,
            "annotate_precip": False,
        },
        {
            "left_col": "delta_effective_albedo_n_minus_i",
            "right_col": "delta_net_sw_n_minus_i_w_m2",
            "left_label": "Δ Albedo",
            "right_label": "Δ Absorbed SW (W/m²)",
            "left_lim": dalb_lim,
            "right_lim": dnsw_lim,
            "legend_labels": ("Δ Albedo", "Δ Absorbed SW"),
            "left_is_delta": True,
            "right_is_delta": True,
            "annotate_precip": True,
        },
        {
            "left_col": "temperature_mean_c",
            "right_col": "delta_net_lw_n_minus_i_w_m2",
            "left_label": "Noah-MP forcing temperature (°C)",
            "right_label": "Δ Estimated net LW (W/m²)",
            "left_lim": temp_lim,
            "right_lim": inlw_lim,
            "legend_labels": ("Noah-MP forcing", "Δ Estimated net LW"),
            "left_is_delta": False,
            "right_is_delta": True,
            "annotate_precip": False,
        },
    ]

    for j, wy in enumerate(TARGET_WATER_YEARS):
        axes[0, j].set_title(f"WY {wy}", fontweight="normal", pad=8)

    for row_idx, spec in enumerate(row_specs):
        for col_idx, wy in enumerate(TARGET_WATER_YEARS):
            sub = yr_data[wy]
            ax = axes[row_idx, col_idx]
            ax_r = ax.twinx()

            line_left, = ax.plot(
                sub["date"], sub[spec["left_col"]],
                color=(color_secondary if spec["left_is_delta"] else color_primary), linewidth=2.5, linestyle="-",
                label=spec["legend_labels"][0],
            )
            line_right, = ax_r.plot(
                sub["date"], sub[spec["right_col"]],
                color=(color_secondary if spec["right_is_delta"] else color_primary), linewidth=2.5, linestyle="--",
                label=spec["legend_labels"][1],
            )

            if row_idx == 2:
                ax_r.fill_between(sub["date"], sub["delta_lw_solar_lower_w_m2"],
                                  sub["delta_lw_solar_upper_w_m2"], color="0.35", alpha=0.18,
                                  linewidth=0, label="Solar integration range")
                ax.text(0.98, 0.04,
                        f"Common area: {sub.common_area_fraction.min():.0%}–{sub.common_area_fraction.max():.0%}",
                        transform=ax.transAxes, ha="right", va="bottom", fontsize=14)
            if col_idx == 0:
                ax.set_ylabel(spec["left_label"])
            if col_idx == 1:
                ax_r.set_ylabel(spec["right_label"])
            ax.set_ylim(*spec["left_lim"])
            ax_r.set_ylim(*spec["right_lim"])

            _format_date_axis(ax, sub["date"])
            ax.spines["right"].set_visible(False)
            ax_r.spines["top"].set_visible(False)

            if spec["left_is_delta"]:
                ax.axhline(0.0, color=color_primary, linewidth=0.8, linestyle=":")
            if spec["right_is_delta"]:
                ax_r.axhline(0.0, color=color_secondary, linewidth=0.8, linestyle=":")

            if spec["annotate_precip"]:
                _add_precip_short_ticks(ax, sub.loc[sub["precip_event"], "date"])
                if col_idx == 0:
                    ax.text(
                        0.98, 0.06,
                        f"Short ticks: precip ≥ {settings.precip_event_threshold_mm_day:g} mm/day",
                        transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=14,
                    )

            if col_idx == 0:
                handles = [line_left, line_right]
                if row_idx == 2:
                    handles.append(ax_r.collections[-1])
                labels = [h.get_label() for h in handles]
                ax.legend(handles, labels, loc="upper left", frameon=True)

    fig.supxlabel("Month", y=0.085, fontsize=18)
    fig.text(0.5, 0.035, "Δ = Noah-MP − iSnobal. Net LW: common snow-present cells; shading: solar integration only.",
             ha="center", va="center", fontsize=14)
    return fig


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to east_river_config.yaml")
    parser.add_argument("--output-subdir", default=DEFAULT_OUTPUT_SUBDIR)
    parser.add_argument("--separation-subdir", default=DEFAULT_SEPARATION_SUBDIR)
    parser.add_argument("--forcing-subdir", default=DEFAULT_FORCING_SUBDIR)
    parser.add_argument("--focus-start", default=FOCUS_START)
    parser.add_argument("--focus-end", default=FOCUS_END)
    parser.add_argument("--precip-threshold-mm-day", type=float, default=DEFAULT_PRECIP_EVENT_THRESHOLD_MM_DAY)
    parser.add_argument("--noah-lw-time-support", choices=["sampled-trapezoid", "hourly-end-mean"],
                        default="sampled-trapezoid",
                        help="Default integrates hourly FIRA samples. Use hourly-end-mean only if your writer explicitly saves preceding-hour means.")
    parser.add_argument("--isnobal-run-config", default="awsm_config_backup.ini",
                        help="Existing config filename inside each run folder; validates hourly daily restart/output intervals.")
    parser.add_argument("--minimum-noah-fsno", type=float, default=0.0,
                        help="Optional stricter snow-fraction screen; default still requires FSNO > 0 at all 25 checkpoints.")
    return parser


def run(args: argparse.Namespace) -> dict[str, Any]:
    cfg = load_config(args.config)
    settings = Settings(
        output_root=cfg.output_dir / str(args.output_subdir),
        separation_dir=cfg.output_dir / str(args.separation_subdir),
        forcing_dir=cfg.output_dir / str(args.forcing_subdir),
        focus_start_month_day=str(args.focus_start),
        focus_end_month_day=str(args.focus_end),
        precip_event_threshold_mm_day=float(args.precip_threshold_mm_day),
    )
    sep, forcing = _load_requested_daily_products(cfg, settings, TARGET_WATER_YEARS)
    merged = _merge_all_products(sep, forcing)
    longwave, raw_lineage = build_daily_longwave(
        cfg, merged.date, args.noah_lw_time_support, args.isnobal_run_config, args.minimum_noah_fsno)
    keys = ["date", "water_year"]
    if (longwave.duplicated(keys).any() or
            set(map(tuple, longwave[keys].to_numpy())) != set(map(tuple, merged[keys].to_numpy()))):
        raise ValueError("Longwave date keys differ from the requested comparison dates.")
    merged = merged.merge(longwave, on=keys, how="left", validate="one_to_one")
    settings.create_tree()
    merged_path = settings.tables_dir / "wy2021_wy2026_daily_merged.csv"
    merged.to_csv(merged_path, index=False)
    longwave.to_csv(settings.tables_dir / "daily_estimated_net_longwave.csv", index=False)
    fig = make_figure(merged, cfg, settings)
    fig_path = _save_figure(fig, cfg, settings, "figure_wy2021_wy2026_daily_comparison")
    input_paths = [settings.separation_dir / "tables/separation_watershed_daily_full_period.csv",
                   settings.forcing_dir / "tables/ldasin_watershed_daily_oct1_jun30.csv"]
    result = {
        "status": "complete", "script_version": SCRIPT_VERSION,
        "output_root": str(settings.output_root), "merged_table": str(merged_path),
        "figure": str(fig_path), "figure_pdf": str(fig_path.with_suffix(".pdf")),
        "config": str(Path(args.config).resolve()),
        "config_sha256": hashlib.sha256(Path(args.config).read_bytes()).hexdigest(),
        "code_sha256": {name: hashlib.sha256((Path(__file__).parent / name).read_bytes()).hexdigest()
                        for name in [Path(__file__).name, "east_river_daily_longwave.py", "east_river_publication.py"]},
        "source_tables": [{"path": str(p.resolve()), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in input_paths],
        "longwave": "Estimated daily net LW on common snow-present support, positive toward surface",
        "noah_lw_time_support": args.noah_lw_time_support,
        "minimum_noah_fsno": args.minimum_noah_fsno,
        "isnobal_run_config": args.isnobal_run_config,
        "raw_input_files": raw_lineage,
        "interval": "UTC midnight-to-midnight; iSnobal weights 23h at D23 + 1h at next D00",
        "solar": "Trapezoidal integral of 25 hourly endpoints, with conditional per-hour left-substep bounds",
        "snow_screen": "iSnobal SWE>0 at D00,D23,nextD00; Noah SWE>0 and FSNO>0 at all 25 hourly samples; all contributing native cells pass",
        "limitations": ["Snow between saved checkpoints is unknown", "Solar range excludes snow-transition, Noah-sampling, canopy, remapping and model errors", "Noah FIRA includes canopy/partial-snow effects; not identical to snow-surface net LW", "Daily selected area varies; inspect common_area_fraction", "Other rows retain their original full-watershed support"],
        "no_common_snow_days": int(longwave.status.eq("no_common_snow_support").sum()),
        "retired_calculation": "iSnobal net_rad minus mean net_solar; unequal averaging intervals and snow-only energy support",
        "difference_sign": "Noah-MP minus iSnobal", "focus_years": list(TARGET_WATER_YEARS),
        "shortwave_conversion": "MJ/m2/day * 1e6/86400 -> W/m2",
        "plotting_standard": "East_River_Plotting_Requirements_JHM_AMS_v6",
    }
    (settings.output_root / "comparison_manifest.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    result = run(args)
    print(pd.Series(result).to_json(indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
