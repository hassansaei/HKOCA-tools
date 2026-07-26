"""Cell-type annotation stage.

Will wrap marker-based / model-based annotation against the kidney
organoid atlas reference (SnapSeed and related helpers).
"""

from __future__ import annotations


def status_message() -> str:
    return (
        "Annotation module — not implemented yet.\n"
        "Will support marker-based and atlas-guided cell-type labeling."
    )


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="hkoca annotation",
        description="Cell-type annotation against the kidney organoid atlas",
    )
    parser.parse_args(argv)
    print(status_message())
    return 1


__all__ = ["main", "status_message"]
