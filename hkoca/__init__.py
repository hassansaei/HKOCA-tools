"""HKOCA modular single-cell toolkit.

Public modules
--------------
pipeline     End-to-end A–Z analysis (orchestrates stages as they are added)
cellbender   Ambient RNA removal only (CellBender)
qc_filter    Harmonization → doublet detection → QC filtering
annotation   Cell-type annotation
integration  Batch integration / atlas projection
"""

from hkoca._version import __version__

__all__ = ["__version__"]
