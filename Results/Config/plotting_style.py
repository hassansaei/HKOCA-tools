"""
Shared plotting configuration for all figure-generating scripts in the atlas paper.

Import this at the top of every fig_*.py script:

    from common.plotting_style import save_figure, FIGURE_DPI
"""

from pathlib import Path
import matplotlib.pyplot as plt

# Single source of truth for output resolution across the whole paper
FIGURE_DPI = 300

# Optional: keep a consistent look across figures. Extend as needed.
plt.rcParams.update({
    "figure.dpi": 100,       # screen/preview resolution while developing
    "savefig.dpi": FIGURE_DPI,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
})


def save_figure(fig, out_path: Path, dpi: int = FIGURE_DPI, tight: bool = True):
    """
    Save a matplotlib figure to disk with consistent settings.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    out_path : Path
        Full path, including filename and extension (.png / .pdf).
    dpi : int
        Resolution, defaults to the paper-wide standard (300).
    tight : bool
        Whether to apply tight bounding box on save.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_path,
        dpi=dpi,
        bbox_inches="tight" if tight else None,
    )
    print(f"[saved] {out_path}  (dpi={dpi})")
