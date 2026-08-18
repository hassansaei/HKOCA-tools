"""Runtime checks for scPoli / PyTorch projection dependencies."""

from __future__ import annotations

PROJECTION_ENV_HINT = (
    "This step runs in the hkoca_projection conda env (PyTorch + scArches)."
)


def _import_projection_stack() -> None:
    """Import torch/scPoli after restoring anndata APIs that scarches 0.6 needs."""
    from hkoca.projection.scpoli import patch_anndata_for_scarches

    patch_anndata_for_scarches()
    import torch  # noqa: F401
    from scarches.models.scpoli import scPoli  # noqa: F401


def projection_stack_available() -> bool:
    try:
        _import_projection_stack()
    except Exception:
        return False
    return True


def require_projection_stack() -> None:
    try:
        _import_projection_stack()
    except Exception as exc:
        raise ImportError(
            f"Failed to import PyTorch/scArches in hkoca_projection: {exc}. "
            + PROJECTION_ENV_HINT
        ) from exc
