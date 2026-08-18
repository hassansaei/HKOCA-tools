"""Runtime checks for scPoli / PyTorch projection dependencies."""

from __future__ import annotations

PROJECTION_ENV_HINT = (
    "This step runs in the hkoca_projection conda env (PyTorch + scArches)."
)


def projection_stack_available() -> bool:
    try:
        import torch  # noqa: F401
        from scarches.models.scpoli import scPoli  # noqa: F401
    except Exception:
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
