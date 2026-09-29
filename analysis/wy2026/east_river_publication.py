"""Shared, verified East River plotting v6 settings for companion scripts."""
from pathlib import Path
import os
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

MODEL_COLORS = {"Noah-MP": "#D62728", "iSnobal": "#1F77B4", "observations": "black"}


def configure_matplotlib(cfg=None):
    """Register the bundled licensed font and reject fallback, including math text."""
    font = Path(os.environ.get("ERW_NIMBUS_FONT", Path(__file__).parent / "fonts" / "NimbusSans-Regular.ttf"))
    if font.is_file():
        fm.fontManager.addfont(str(font))
    resolved = fm.findfont(fm.FontProperties(family="Nimbus Sans"), fallback_to_default=False)
    if fm.FontProperties(fname=resolved).get_name() != "Nimbus Sans":
        raise RuntimeError(f"Nimbus Sans font resolution failed: {resolved}")
    mpl.rcParams.update({
        "font.family": "Nimbus Sans", "font.size": 18,
        "axes.labelsize": 18, "axes.titlesize": 20,
        "axes.titleweight": "normal", "axes.labelweight": "normal",
        "xtick.labelsize": 18, "ytick.labelsize": 18, "legend.fontsize": 18,
        "lines.linewidth": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans",
        "mathtext.it": "Nimbus Sans:italic", "mathtext.bf": "Nimbus Sans",
        "legend.frameon": True, "legend.facecolor": "white",
        "legend.edgecolor": "0.4", "legend.framealpha": 0.9,
        "legend.fancybox": False, "axes.grid": False,
        "axes.linewidth": 1.0, "xtick.major.width": 1.0,
        "ytick.major.width": 1.0, "xtick.major.size": 5, "ytick.major.size": 5,
    })
    return resolved


def add_panel_labels(fig):
    """Label distinct visible data panels in reading order (ER-PLOT-027).

    Twin axes occupy one panel; colorbars and hidden/empty disabled axes are
    not panels. Labels sit just above each panel's upper-left corner so they
    do not compete with inside-panel legends. Repeated exports are idempotent.
    """
    marker = "east-river-panel-label"
    for ax in fig.axes:
        for text in list(ax.texts):
            if text.get_gid() == marker:
                text.remove()
    panels = []
    for ax in fig.axes:
        if (not ax.get_visible() or ax.get_label() == "<colorbar>"
                or getattr(ax, "_colorbar", None) is not None
                or (not ax.axison and not ax.has_data())):
            continue
        # Both exact panel positions and Matplotlib's explicit twin grouping
        # are checked, including figures whose colorbars modify GridSpecs.
        bounds = tuple(round(v, 8) for v in ax.get_position(original=True).bounds)
        twins = ax._twinned_axes.get_siblings(ax)
        if any(other in twins or other_bounds == bounds for other, other_bounds in panels):
            continue
        panels.append((ax, bounds))
    panels.sort(key=lambda item: (-round(item[1][1]+item[1][3], 6), item[1][0]))
    if len(panels) < 2:
        return []

    def letters(index):
        value = index + 1
        result = ""
        while value:
            value, remainder = divmod(value-1, 26)
            result = chr(97+remainder) + result
        return result

    labels = []
    for index, (ax, _) in enumerate(panels):
        text = ax.annotate(f"({letters(index)})", xy=(0, 1), xycoords="axes fraction",
                           xytext=(0, 8), textcoords="offset points",
                           ha="left", va="bottom", fontsize=18,
                           fontfamily="Nimbus Sans", fontweight="normal", color="black",
                           annotation_clip=False, clip_on=False, zorder=20)
        text.set_gid(marker)
        text.set_in_layout(True)
        labels.append(text)
    return labels


def save_figure(fig, path, *, close=True):
    """Export matching 600-dpi RGB PNG and embedded-font vector PDF."""
    from PIL import Image
    add_panel_labels(fig)
    path = Path(path).with_suffix("")
    path.parent.mkdir(parents=True, exist_ok=True)
    png, pdf = path.with_suffix(".png"), path.with_suffix(".pdf")
    fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    with Image.open(png) as im:
        rgb = im.convert("RGB")
    rgb.save(png, dpi=(600, 600))
    rgb.close()
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    if close:
        plt.close(fig)
    return png, pdf
