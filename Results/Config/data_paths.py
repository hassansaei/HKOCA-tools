"""
Central place for all input data locations and output roots.

Edit this file, not the individual fig_*.py scripts, when a data path changes
(e.g. moving to a new server, updating the harmonized data release).
"""

from pathlib import Path

# --- Raw / processed data roots ---
HARMONIZED_DATA_DIR = Path(
    "/data-master/pure-workspace/labss/hmami/new_data/Updated_data/updated_harmonized_files"
)

# --- Paper output root (each step folder has its own outputs/ subfolder) ---
PAPER_ROOT = Path(__file__).resolve().parent.parent  # .../paper/
