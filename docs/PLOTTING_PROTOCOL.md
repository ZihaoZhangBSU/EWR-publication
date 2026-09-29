# East River Plotting Requirements — JHM/AMS Publication Standard

**Status:** Current project-wide plotting standard  
**Version:** 6  
**Date:** 24 September 2026  
**Target journal:** *Journal of Hydrometeorology* (American Meteorological Society)  
**Primary internal visual/code reference:** `Watershed_comparison (1).ipynb`  
**Default plotting font:** Nimbus Sans

## 1. Purpose and authority

This document defines the general visual and publication conventions for East River project figures. It is intended to be source-friendly: requirements have stable IDs, project decisions are separated from external AMS guidance, and source provenance is recorded explicitly.

When requirements differ, use this authority order:

1. Explicit instructions from the user for the current figure.
2. This plotting-requirements document.
3. Figure-specific scientific needs.
4. Existing East River plotting code/templates.
5. AMS/JHM publication guidance.

Plotting requirements govern visual presentation and export. Scientific definitions, filtering, masking, units, sign conventions, spatial support, temporal aggregation, and statistical calculations must come from the validated workflow or from an explicitly approved analysis.

## 2. Hard project requirements

**ER-PLOT-001 — Default font.** Use **Nimbus Sans** as the intended font for all East River figure text.

**ER-PLOT-002 — Font verification.** For final publication figures, verify that Nimbus Sans is actually available and being used. Do not silently accept a Matplotlib fallback.

**ER-PLOT-003 — Base font size.** Use **18 pt** as the standard base figure-text size for the established East River development template unless a figure-specific layout requires adjustment.

**ER-PLOT-004 — Panel-title size.** Water-year panel titles should use **20 pt**, normal weight, centered, with `pad=8`.

**ER-PLOT-005 — Font weight.** Use normal-weight typography by default. Avoid bold text unless emphasis is genuinely needed.

**ER-PLOT-006 — Primary line width.** Use **2.5 pt** for primary iSnobal and Noah-MP comparison curves in standard time-series figures.

**ER-PLOT-007 — Output formats.** Every final figure must be written in **both PNG and PDF** formats using the same filename stem.

**ER-PLOT-008 — PNG resolution.** Save the PNG at **600 dpi** by default.

**ER-PLOT-009 — PDF quality.** Save the PDF as vector output where possible and embed fonts.

**ER-PLOT-010 — Main figure title.** Do not add a main figure title automatically. Add one only when explicitly requested or scientifically necessary.

### AMS font compatibility note

Nimbus Sans is a **project choice**. Current AMS figure guidance says to use the same font across figures and lists Arial, Helvetica, Times, Symbol, Mathematical Pi, and European Pi as standard fonts. Nimbus Sans is **not explicitly named** in that AMS list. This is not a stated prohibition, but the distinction must not be hidden. For submission, keep the font embedded in the PDF and verify that the rendered figure remains correct in the AMS production workflow. See `SRC-AMS-02`.

## 3. Statistical annotations

**ER-PLOT-011 — Statistical annotations are optional.** Statistical annotations are **not a general requirement**. Include them only when they materially help interpret the figure or when the user explicitly requests them.

**ER-PLOT-012 — Metric abbreviations.** When statistical annotations are included:

- use `r` for **Pearson correlation**;
- use `MB` for **mean bias**.

Do not routinely spell out “Pearson correlation” or “mean bias” inside a panel when `r` and `MB` are sufficient.

Preferred example when both statistics are appropriate:

```text
r = 0.87
MB = +12.4 mm
```

**ER-PLOT-013 — Standard statistics position.** When a compact statistics block is used, place it near the upper-right of the panel by default:

```python
ax.text(
    0.97,
    0.95,
    text,
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=18,
)
```

This position is adaptable. Move the annotation when the default placement would overlap important data.

## 4. Typography implementation

Preferred Matplotlib publication setup:

```python
import matplotlib as mpl
import matplotlib.font_manager as fm

FONT_NAME = "Nimbus Sans"

available_names = {f.name for f in fm.fontManager.ttflist}
if FONT_NAME not in available_names:
    raise RuntimeError(
        "Nimbus Sans is required for East River publication figures "
        "but is not available in the current environment."
    )

mpl.rcParams.update({
    "font.family": FONT_NAME,
    "font.size": 18,
    "axes.titlesize": 20,
    "axes.titleweight": "normal",
    "axes.labelweight": "normal",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
```

The installation method for Nimbus Sans is environment-specific and is intentionally outside this plotting standard.

## 5. Standard multipanel geometry

For comparable multi-water-year figures, use three columns:

```python
ncols = 3
nrows = math.ceil(len(years) / ncols)

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(21, 5.5 * nrows),
    sharey=True,
    squeeze=False,
)
```

For nine water years, the established development canvas is `figsize=(21, 16.5)`.

Default spacing for comparable 3-column water-year layouts:

```python
fig.subplots_adjust(
    top=0.95,
    bottom=0.15,
    left=0.07,
    hspace=0.30,
    wspace=0.16,
)
```

Hide unused axes:

```python
for ax in axes.flat[len(years):]:
    ax.set_visible(False)
```

These values are defaults for comparable figures, not constraints for maps, 2×2 summaries, mechanism-chain figures, distribution figures, or other scientifically different layouts.

## 6. Shared axis labels

When panels share the same axis meaning, use figure-level shared labels rather than repeating labels in every panel.

Established positions:

```python
fig.supxlabel("Month", y=0.095)
fig.supylabel(ylabel, x=0.015)
```

The authoritative shared y-label position for the watershed template is **`x=0.015`**.

## 7. Panel titles

For water-year panels:

```python
ax.set_title(
    f"WY {wy}",
    fontsize=20,
    fontweight="normal",
    loc="center",
    pad=8,
)
```

Panel titles are distinct from a main figure title.

### Multipanel identifiers

**ER-PLOT-027 — Required panel labels.** Every figure with two or more distinct visible data panels must label them `(a)`, `(b)`, `(c)`, … in reading order: left to right across a row, then top to bottom. Restart at `(a)` for each figure. Continue `(aa)`, `(ab)`, … after `(z)` if needed.

- Use Nimbus Sans, 18 pt, normal weight, in black, consistently across all panels.
- Place the label just above the upper-left corner of its panel, by default 8 points above the top border. Keep centered water-year/descriptive titles separate. Adjust spacing when necessary to prevent overlap with titles, adjacent panels, legends, statistics, or data.
- A twin-y or secondary axis belongs to the same panel and receives no additional label. A colorbar is not a data panel. Hidden/unused subplot slots receive no label and consume no letter.
- A single-panel figure does not need `(a)`, even if it has a colorbar or secondary axis.
- Include each label exactly once. Remove older inline `(a)`/`(b)` title prefixes when adding separate corner labels; preserve the descriptive title text.
- Preserve labels in both PNG and vector PDF exports and check for clipping at the final figure size.

This is an explicit project decision requested by the user on 24 September 2026; it is not presented as a new external AMS requirement. The shared `east_river_publication.save_figure` helper in the comparison package adds the labels automatically. Other project plotting code must apply the same rule before export.

Example when axes have already been ordered and filtered to actual data panels:

```python
for label, ax in zip(["(a)", "(b)", "(c)", "(d)"], data_axes):
    ax.annotate(
        label, xy=(0, 1), xycoords="axes fraction",
        xytext=(0, 8), textcoords="offset points",
        ha="left", va="bottom", fontsize=18,
        fontfamily="Nimbus Sans", fontweight="normal",
        color="black", annotation_clip=False, clip_on=False,
    )
```

## 8. Legend

For standard multipanel comparison figures, use the established East River legend treatment: a boxed legend inside the first subplot at the upper left when it does not obscure important data.

```python
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
```

Use the base figure font size unless a figure-specific layout requires a modest reduction.

Legend placement is adaptable. If the upper-left location obscures scientifically important content, move the legend while preserving the same clean, boxed treatment when practical.

**ER-PLOT-014 — Reference-image legends do not override the project legend.** When a visual reference is used for boxplot/bar construction, copy the relevant plotting geometry or marker treatment only. Keep the East River legend convention unless the user explicitly requests a different legend.

## 9. Axes and ticks

For standard line/scatter figures:

```python
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
```

Avoid unnecessary gridlines. Add reference lines such as `y=0` only when they aid scientific interpretation.

## 10. Water-year time axes

For full water-year hydrographs:

- show October through September;
- use abbreviated month names;
- use major month ticks every two months by default;
- remove the top and right spines.

Reference implementation:

```python
start, end = water_year_bounds(int(wy))
ax.set_xlim(start, end)

ax.xaxis.set_major_locator(
    mdates.MonthLocator(interval=2)
)
ax.xaxis.set_major_formatter(
    mdates.DateFormatter("%b")
)
```

For shorter focus windows, use a denser locator when needed for readability.

## 11. Axis limits

For physically nonnegative snow-state variables such as SWE and snow depth, normally start the y-axis at zero. Use clean upper bounds with modest headroom rather than arbitrary maxima.

Established rounding for the watershed comparison template:

- SWE: 50 mm increments;
- snow depth: 100 mm increments.

For signed differences, anomalies, errors, or biases, use a visible zero reference and use symmetric limits around zero when that improves scientific interpretation.

## 12. Color and source identity

**ER-PLOT-015 — Source identity overrides reference-image palettes.** Use the established East River source/model colors even when a visual reference uses another palette. Do not copy arbitrary blue/salmon/green colors from an example figure when they conflict with East River source identity.

Default source identity:

- observations / ASO: **black**;
- iSnobal: **`#1F77B4`**;
- Noah-MP: **`#D62728`**.

Use the established `MODEL_COLORS` definitions when available rather than redefining colors locally.

Maintain source identity consistently between figures.

Do not rely on color alone when line style, marker shape, hatching, or another redundant cue would make an important comparison more accessible.

AMS requires color image files to be submitted in RGB and strongly encourages color choices that remain interpretable for readers with color-vision impairment. See `SRC-AMS-04`.

## 13. Grouped distribution and boxplot protocol

This section governs grouped distribution plots such as elevation-bin comparisons. The approved visual language is a **notched boxplot** when the notch has statistical meaning.

**ER-PLOT-016 — Data-driven notches.** When notched boxplots are used, the notch geometry must be computed from the actual plotted distribution or from a documented statistical confidence interval. Do **not** manually force all boxes to have the same notch shape. Different notch shapes are expected when the distributions or sample sizes differ.

If a notch is not statistically meaningful for a particular figure, use a regular rectangular boxplot instead of decorative notching.

**ER-PLOT-017 — Group order.** When observations/ASO, iSnobal, and Noah-MP are shown together within an elevation bin or other category, use the default left-to-right order:

1. observations / ASO;
2. iSnobal;
3. Noah-MP.

Keep the same order across all bins/panels in a figure.

**ER-PLOT-018 — Box, whisker, and cap styling.** Use black box outlines, whiskers, and caps. A standard visual starting point is approximately **1.5 pt** for these strokes.

The **statistical definition of whisker extent and outlier handling is not a style rule**. It must come from the scientific definition for the figure and must not be silently changed to match a visual reference.

**ER-PLOT-019 — Median.** When a median is shown, use a **black horizontal median line**, approximately **2.0 pt** in the established grouped-distribution style.

**ER-PLOT-020 — Mean.** When a mean is scientifically appropriate to display, show it as a **yellow diamond with a black edge**. The mean marker is optional; do not add it when the mean is not relevant to the figure.

**ER-PLOT-021 — Observation fill.** Preserve observations/ASO as the black source identity. For filled boxes, black may be shown with reduced opacity so the box appears light gray while retaining a black outline. The approved demo used approximately `alpha=0.30`.

**ER-PLOT-022 — Model box fill.** Use the East River model colors with black outlines. For filled distribution boxes, a partially transparent fill around `alpha=0.70–0.75` is the preferred starting range.

**ER-PLOT-023 — Background and grid.** For grouped boxplot/bar-style summary figures, use a light neutral-gray plotting background and subtle horizontal dashed gridlines behind the data. The approved demo uses approximately:

```python
ax.set_facecolor("0.94")
ax.yaxis.grid(
    True,
    linestyle="--",
    linewidth=0.8,
    alpha=0.45,
)
ax.xaxis.grid(False)
ax.set_axisbelow(True)
```

Do not use vertical gridlines by default for these grouped categorical figures.

**ER-PLOT-024 — Boxplot styling template.** The following snippet captures the visual protocol. `WHISKER_RULE` and `SHOW_FLIERS` are intentionally left as scientific choices rather than hard-coded style choices.

```python
bp = ax.boxplot(
    values,
    positions=positions,
    widths=0.20,
    notch=True,
    patch_artist=True,
    whis=WHISKER_RULE,      # define from the scientific specification
    showfliers=SHOW_FLIERS, # define from the scientific specification
    manage_ticks=False,
    medianprops={"color": "black", "linewidth": 2.0},
    whiskerprops={"color": "black", "linewidth": 1.5},
    capprops={"color": "black", "linewidth": 1.5},
    boxprops={"edgecolor": "black", "linewidth": 1.5},
)
```

For the mean marker when used:

```python
ax.scatter(
    x,
    mean_value,
    marker="D",
    s=88,
    facecolor="yellow",
    edgecolor="black",
    linewidth=1.0,
    zorder=5,
)
```

**ER-PLOT-025 — No automatic statistical text.** Do not add `r`, `MB`, sample sizes, p-values, or other text annotations merely because the figure is a statistical distribution plot. Statistical annotations remain optional under `ER-PLOT-011`.

## 14. True bar-chart adaptation

The approved boxplot reference should not be interpreted to mean that every “bar-style” categorical figure must be a boxplot.

**ER-PLOT-026 — True bars.** When the scientific quantity is a single aggregate per category and a true bar chart is appropriate:

- use the established East River source/model colors;
- use black bar edges;
- retain the light-gray background and subtle horizontal dashed gridlines when appropriate;
- retain the existing East River legend convention;
- do not add notches, median lines, or mean diamonds unless those elements have an explicit scientific meaning for that figure.

## 15. Mandatory dual export

Every final figure must write both PNG and PDF:

```python
png_path = out_dir / f"{stem}.png"
pdf_path = out_dir / f"{stem}.pdf"

fig.savefig(
    png_path,
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.05,
)

fig.savefig(
    pdf_path,
    bbox_inches="tight",
    pad_inches=0.05,
)
```

Neither format is optional.

For a multipanel figure, export the complete multipanel figure as one PNG and one PDF rather than exporting each panel separately unless the user explicitly requests separate panel files.

## 16. JHM/AMS publication requirements relevant to this project

The *Journal of Hydrometeorology* is an American Meteorological Society journal. Use current AMS figure guidance for publication-format requirements unless JHM publishes a more specific requirement for a particular submission situation.

**AMS-001 — Figure completeness.** Submit one complete figure per figure file. Multiple-panel figures that appear as one figure should be combined into a single file whenever possible. See `SRC-AMS-03`.

**AMS-002 — Consistent fonts.** Use the same font across figures. AMS lists Arial, Helvetica, Times, Symbol, Mathematical Pi, and European Pi as standard fonts. See `SRC-AMS-02`.

**AMS-003 — Font embedding.** Embed fonts in figure documents/vector outputs where applicable. See `SRC-AMS-02`.

**AMS-004 — Raster resolution.** AMS guidance gives approximately 300 dpi for tones/halftones, 600–900 dpi for combination images, and 1000–1200 dpi for line art. The East River 600-dpi PNG default is appropriate for many line-plus-color scientific figures, while the companion PDF preserves vector content. See `SRC-AMS-02`.

**AMS-005 — Color mode.** Submit color image files in original RGB color. See `SRC-AMS-04`.

**AMS-006 — Accessibility.** Prepare colors so that readers with color-vision impairment can distinguish the scientific content. See `SRC-AMS-04`.

**AMS-007 — Final figure widths.** AMS lists standard final widths of 3.2 in for one column, 5.5 in for two columns, and 6.5 in for more than two columns. In-figure text and detail can affect final size. See `SRC-AMS-01`.

**AMS-008 — Publication resizing.** AMS technical editors may resize figures for layout and clarity. See `SRC-AMS-01` and `SRC-AMS-03`.

**AMS-009 — Readability.** AMS author guidance emphasizes that figure labels, contours, colorbars, and other text must remain readable after typesetting; a font may need to look large in the original figure because the published figure will be smaller. Sans-serif fonts are recommended where possible. See `SRC-AMS-05`.

## 17. Publication-size readability check

The 21-in-wide East River canvas is a **development geometry**, not the expected printed width.

Before finalizing a figure, inspect it at its intended AMS publication width. For a complex 3×3 figure, check readability at approximately **6.5 in wide**.

Do not reduce a figure until labels, tick text, legends, annotations, or symbols become difficult to read.

The 18-pt base and 20-pt panel-title values describe the established East River development template. If a figure is authored directly at final journal dimensions, adjust source font sizes as needed to preserve equivalent final readability.

## 18. Scientific plotting guardrails

Plotting code must not silently:

- redefine a scientific metric;
- change a bias sign convention;
- change temporal aggregation;
- change spatial support or weighting;
- substitute a variable;
- change unit conversions;
- create a new data-filtering rule;
- infer missing intermediate data from a screenshot or filename;
- change the statistical definition of a boxplot whisker, outlier, notch, mean, or median merely to reproduce a visual shape.

Use validated workflow outputs and exact metric definitions whenever available.

If a figure requires a new calculation, implement or document that calculation explicitly before treating it as plotting-only logic.

## 19. Final pre-export checklist

- Multipanel figures have exactly one `(a)`, `(b)`, … label per visible data panel in reading order; no colorbar/twin-axis labels, duplicate title prefixes, overlaps, or clipped letters.

- [ ] Nimbus Sans is the intended figure font.
- [ ] Nimbus Sans availability has been verified; no silent fallback occurred.
- [ ] Typography is normal weight unless emphasis is necessary.
- [ ] Base text follows the 18-pt development standard when applicable.
- [ ] Water-year panel titles are 20 pt when applicable.
- [ ] Primary iSnobal/Noah-MP comparison lines are 2.5 pt when applicable.
- [ ] Statistical annotations are included only when useful or requested.
- [ ] When statistics are shown, `r` means Pearson correlation and `MB` means mean bias.
- [ ] Statistics do not obscure important data.
- [ ] Shared labels are used when appropriate.
- [ ] Established shared y-label position is `x=0.015` for the watershed template.
- [ ] No unnecessary main figure title is present.
- [ ] Legend placement follows the East River convention unless explicitly changed.
- [ ] Top/right spines are removed where appropriate.
- [ ] Axis limits are scientifically interpretable.
- [ ] Observation/ASO, iSnobal, and Noah-MP colors follow East River source identity.
- [ ] Reference-image palettes have not overridden project color conventions.
- [ ] For notched boxplots, notch geometry is data-driven and has not been manually standardized.
- [ ] Boxplot whisker/outlier definitions come from the scientific specification, not from visual imitation.
- [ ] Mean diamonds and median lines are shown only when scientifically appropriate.
- [ ] Grouped boxplot/bar figures use the approved background/grid treatment when appropriate.
- [ ] Color is not the sole cue when accessibility requires redundancy.
- [ ] Figure remains readable at expected AMS/JHM publication width.
- [ ] One complete multipanel figure is exported per figure file.
- [ ] PNG output exists at 600 dpi.
- [ ] PDF output exists.
- [ ] Fonts are embedded in the PDF.
- [ ] Color output is RGB.
- [ ] Plotting code has not redefined scientific calculations.

## 20. Source provenance

### SRC-ERW-01 — Internal project plotting reference

`Watershed_comparison (1).ipynb`

Used as the primary internal source for the approved visual template, including:

- 3-column dynamic water-year layout;
- `figsize=(21, 5.5 * nrows)`;
- 20-pt water-year titles;
- 2.5-pt primary model lines;
- shared x-label at `y=0.095`;
- shared y-label at `x=0.015`;
- upper-right statistics position when statistics are included;
- boxed upper-left first-panel legend;
- 1.0-pt left/bottom axes;
- 1.0-pt major-tick width and 5-point tick length;
- 600-dpi PNG output;
- companion PDF output.

The following are explicit later project decisions and therefore supersede older template wording where necessary:

- Nimbus Sans is the default plotting font.
- Statistical annotations are optional.
- When shown, `r` means Pearson correlation and `MB` means mean bias.
- Both PNG and PDF are mandatory outputs.
- Reference-image palettes do not override East River source colors.
- Observations/ASO use black source identity, iSnobal uses `#1F77B4`, and Noah-MP uses `#D62728`.
- The established East River legend convention remains in force for grouped boxplot/bar figures.
- Notched boxplots use **data-driven notch geometry** when notches are statistically meaningful.
- In grouped distribution figures, black outlines/whiskers/caps, black medians, optional yellow-diamond means, a light-gray background, and subtle dashed horizontal gridlines are the approved visual treatment.

### SRC-ERW-02 — Workflow elevation-band definitions

The East River workflow defines elevation bands as:

- Lower: `< 2896 m`;
- Middle: `>= 2896 m and < 3353 m`;
- Upper: `>= 3353 m`;
- All: all valid watershed cells.

These elevation definitions are scientific/data definitions from the workflow and are not styling choices.

### SRC-AMS-01 — AMS Figure Sizing

American Meteorological Society, “Figure Sizing.”  
https://www.ametsoc.org/ams/publications/author-information/figure-information-for-authors/figure-sizing/

Supports the standard final widths of 3.2, 5.5, and 6.5 in and the note that in-figure detail may alter the final production size.

### SRC-AMS-02 — AMS Microsoft Office Images / general figure preparation

American Meteorological Society, “Microsoft Office Images.”  
https://www.ametsoc.org/ams/publications/author-information/figure-information-for-authors/microsoft-office-images/

Supports the same-font-across-figures guidance, the AMS list of standard fonts, font embedding, and the 300 / 600–900 / 1000–1200 dpi resolution guidance.

### SRC-AMS-03 — AMS Formatting and Manuscript Components

American Meteorological Society, “Formatting and Manuscript Components.”  
https://www.ametsoc.org/ams/publications/author-information/formatting-and-manuscript-components/

Supports combining multipanel figures into one complete figure file, providing source files at revision, and allowing technical editors to resize figures.

### SRC-AMS-04 — AMS Color Mode (RGB and CMYK)

American Meteorological Society, “Color Mode (RGB and CMYK).”  
https://www.ametsoc.org/ams/publications/author-information/figure-information-for-authors/color-mode-rgb-and-cmyk/

Supports RGB submission for color figures and color-vision accessibility guidance.

### SRC-AMS-05 — AMS author guidance on readable figures

American Meteorological Society, “How to Be a More Effective Author.” *Monthly Weather Review* 150(11).  
https://journals.ametsoc.org/view/journals/mwre/150/11/MWR-D-22-0277.1.xml

Supports large readable in-figure text, sufficient quality to avoid pixelation, use of sans-serif fonts where possible, and avoiding inaccessible color choices.

**External-source verification date:** 10 September 2026.

## 21. Change log

**v6 — 24 September 2026**

- Added ER-PLOT-027: required lowercase parenthesized identifiers on every multipanel figure.
- Defined reading order, 18-pt Nimbus Sans styling, default placement, and handling of twins, colorbars, hidden panels, and repeated exports.
- Preserved all v5 scientific-scope boundaries, colors, typography, legend, statistical, and export requirements.


**v5 — 11 September 2026**

- Added the approved grouped distribution / notched-boxplot protocol.
- Added the rule that notch geometry must remain data-driven rather than manually standardized.
- Preserved East River source colors: observations/ASO black, iSnobal blue, Noah-MP red.
- Added black box/whisker/cap outlines, black median lines, optional yellow-diamond means, light-gray plotting background, and subtle dashed horizontal gridlines as the preferred grouped-distribution visual treatment.
- Explicitly preserved the existing East River legend convention rather than copying legend placement from the external boxplot reference image.
- Separated boxplot visual styling from scientific choices such as whisker extent and outlier handling.
- Added a true-bar-chart adaptation so the boxplot protocol is not incorrectly applied to aggregate bar charts.
- Preserved Nimbus Sans, 18-pt base font, 20-pt water-year titles, 2.5-pt primary time-series lines, `x=0.015` shared y-label position, optional statistical annotations, `r`/`MB` notation, and mandatory PNG + PDF export.

**v4 — 10 September 2026**

- Changed the project default plotting font from Helvetica to **Nimbus Sans**.
- Added an explicit AMS compatibility note because Nimbus Sans is not named in the AMS standard-font list.
- Preserved the existing 18-pt base font, 20-pt water-year title, 2.5-pt primary line width, `x=0.015` shared y-label position, optional statistical annotations, `r`/`MB` notation, and mandatory PNG + PDF export.
