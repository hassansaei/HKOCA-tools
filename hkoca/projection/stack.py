"""Runtime checks for scPoli / PyTorch projection dependencies."""

from __future__ import annotations

PROJECTION_ENV_HINT = (
    "The hkoca_projection conda env must be built ahead of time from "
    "conda/environment_projection.yaml (PyTorch + scArches). "
    "hkoca is shared from the harmonize install via PYTHONPATH; "
    "do not pip install packages when running the pipeline."
)


def projection_stack_available() -> bool:
    try:
        import torch  # noqa: F401
        from scarches.models.scpoli import scPoli  # noqa: F401
    except ImportError:
        return False
    return True


def require_projection_stack() -> None:
    try:
        import torch  # noqa: F401
        from scarches.models.scpoli import scPoli  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "scPoli projection requires PyTorch and scArches in hkoca_projection. "
            + PROJECTION_ENV_HINT
        ) from exc
