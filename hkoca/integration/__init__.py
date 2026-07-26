"""Batch integration and atlas projection stage.

Will wrap integration methods (e.g. scPoli / scVI family) and projection
of query datasets onto the kidney organoid atlas.
"""

from __future__ import annotations


def status_message() -> str:
    return (
        "Integration module — not implemented yet.\n"
        "Will support batch integration and projection onto the organoid atlas."
    )


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="hkoca integration",
        description="Batch integration and atlas projection",
    )
    parser.parse_args(argv)
    print(status_message())
    return 1


__all__ = ["main", "status_message"]
