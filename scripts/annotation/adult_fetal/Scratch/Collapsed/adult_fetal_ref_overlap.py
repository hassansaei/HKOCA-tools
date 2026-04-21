#!/usr/bin/env python3
"""
pipeline.py
-----------
Unified Adult-Fetal Kidney Reference Pipeline.

Replaces:
  kidney_reference_stage1_clean.ipynb
  hpa_marker_validation.ipynb
  overlap.ipynb

Three explicit phases, run in sequence or individually:
  Phase 1 — Load, preprocess, marker extraction, fetal rescue
    Phase 2 — Adult HPA validation + ranking
  Phase 3 — Overlap analysis + SnapSeed YAML export

Usage
-----
  # Validate inputs without processing:
  python pipeline.py --config config.yaml --dry-run

    # Full run (default runs all phases 1 2 3):
  python pipeline.py --config config.yaml

    # Selected phases only:
    python pipeline.py --config config.yaml --only-phases 1
    python pipeline.py --config config.yaml --only-phases 2 3
    python pipeline.py --config config.yaml --skip-phases 3

    # Force rerun (disable resume via hidden .meta.json):
    python pipeline.py --config config.yaml --no-resume

  # Verify output contract after a completed run:
  python pipeline.py --config config.yaml --verify
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import sys
import textwrap
import time
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pandas_types

    PandasDataFrame = pandas_types.DataFrame
    PandasSeries = pandas_types.Series
else:
    PandasDataFrame = Any
    PandasSeries = Any

import yaml

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None  # type: ignore[assignment]

try:
    import numpy as np
except ModuleNotFoundError:
    np = None  # type: ignore[assignment]

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None  # type: ignore[assignment]

try:
    import seaborn as sns
except ModuleNotFoundError:
    sns = None  # type: ignore[assignment]


# =============================================================================
# SECTION 1 — Logging & Config
# =============================================================================

def setup_logging() -> logging.Logger:
    logger = logging.getLogger("kidney_pipeline")
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(
        fmt="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logger.addHandler(handler)
    return logger


log = setup_logging()


def add_file_logging(log_file: Path) -> None:
    """Attach a file handler once to preserve per-run logs for reproducibility."""
    resolved = log_file.resolve()
    for handler in log.handlers:
        if isinstance(handler, logging.FileHandler) and Path(handler.baseFilename).resolve() == resolved:
            return

    resolved.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(resolved, encoding="utf-8")
    handler.setLevel(log.level)
    handler.setFormatter(logging.Formatter(
        fmt="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    log.addHandler(handler)
    log.info("File logging enabled: %s", resolved)


def _require_dict(cfg: dict[str, Any], key: str) -> dict[str, Any]:
    section = cfg.get(key)
    if not isinstance(section, dict):
        raise ValueError(f"Config section '{key}' must be a mapping/dict.")
    return section


def _require_keys(section: dict[str, Any], section_name: str, keys: list[str]) -> None:
    missing = [k for k in keys if k not in section]
    if missing:
        raise KeyError(f"Config section '{section_name}' missing required key(s): {missing}")


def _require_positive_int(value: Any, key_name: str) -> int:
    if not isinstance(value, int) or value <= 0:
        raise ValueError(f"Config key '{key_name}' must be a positive integer (got: {value!r}).")
    return value


def _get_package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for dist_name in ["anndata", "scanpy", "numpy", "pandas", "matplotlib", "seaborn", "PyYAML"]:
        try:
            versions[dist_name] = importlib_metadata.version(dist_name)
        except importlib_metadata.PackageNotFoundError:
            versions[dist_name] = None
    return versions


def _validate_manual_mapping_seed(seed: Any) -> None:
    if not isinstance(seed, list):
        raise ValueError("Config key 'overlap.manual_mapping_seed' must be a list when provided.")

    required = [
        "fetal_cell_type",
        "adult_cell_type_candidates",
        "mapping_type",
        "confidence",
        "is_broad_generic",
    ]
    for i, row in enumerate(seed):
        if not isinstance(row, dict):
            raise ValueError(f"overlap.manual_mapping_seed[{i}] must be a mapping/dict.")
        _require_keys(row, f"overlap.manual_mapping_seed[{i}]", required)
        if not isinstance(row["adult_cell_type_candidates"], list) or not row["adult_cell_type_candidates"]:
            raise ValueError(
                f"overlap.manual_mapping_seed[{i}].adult_cell_type_candidates "
                "must be a non-empty list."
            )


def _validate_overlap_acceptance_criteria(criteria: Any) -> None:
    if criteria is None:
        return
    if not isinstance(criteria, dict):
        raise ValueError("Config key 'overlap.acceptance_criteria' must be a mapping/dict when provided.")

    allowed = {
        "min_best_mapped_jaccard",
        "min_best_mapped_shared_genes",
        "expected_shared_mapped_min",
        "expected_shared_mapped_max",
        "expected_adult_only_min",
        "expected_adult_only_max",
        "expected_fetal_only_min",
        "expected_fetal_only_max",
    }
    unknown = sorted(set(criteria.keys()) - allowed)
    if unknown:
        raise ValueError(
            "Unknown keys in overlap.acceptance_criteria: "
            f"{unknown}. Allowed keys: {sorted(allowed)}"
        )

    for key in [
        "min_best_mapped_jaccard",
        "min_best_mapped_shared_genes",
        "expected_shared_mapped_min",
        "expected_shared_mapped_max",
        "expected_adult_only_min",
        "expected_adult_only_max",
        "expected_fetal_only_min",
        "expected_fetal_only_max",
    ]:
        if key not in criteria or criteria[key] is None:
            continue
        val = criteria[key]
        if not isinstance(val, (int, float)):
            raise ValueError(f"Config key 'overlap.acceptance_criteria.{key}' must be numeric (got: {val!r}).")
        if float(val) < 0:
            raise ValueError(f"Config key 'overlap.acceptance_criteria.{key}' must be >= 0 (got: {val!r}).")

    for lo_key, hi_key in [
        ("expected_shared_mapped_min", "expected_shared_mapped_max"),
        ("expected_adult_only_min", "expected_adult_only_max"),
        ("expected_fetal_only_min", "expected_fetal_only_max"),
    ]:
        lo = criteria.get(lo_key)
        hi = criteria.get(hi_key)
        if lo is not None and hi is not None and float(lo) > float(hi):
            raise ValueError(
                f"Config keys 'overlap.acceptance_criteria.{lo_key}' and "
                f"'overlap.acceptance_criteria.{hi_key}' are inconsistent: min > max."
            )


def _validate_dataset_manifest_config(manifest_cfg: Any) -> None:
    if manifest_cfg is None:
        return
    if not isinstance(manifest_cfg, dict):
        raise ValueError("Config key 'dataset_manifest' must be a mapping/dict when provided.")
    for stage in ["adult", "fetal"]:
        if stage in manifest_cfg and not isinstance(manifest_cfg[stage], dict):
            raise ValueError(f"Config key 'dataset_manifest.{stage}' must be a mapping/dict when provided.")


def _normalize_cfg_label(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value).strip().lower())


def _validate_cross_section_consistency(
    hpa_validation: dict[str, Any],
    overlap: dict[str, Any],
) -> None:
    manual_override_raw = hpa_validation.get("manual_override_cell_types", [])
    manual_override_labels = [
        _normalize_cfg_label(x)
        for x in manual_override_raw
        if isinstance(x, str) and str(x).strip()
    ]
    if not manual_override_labels:
        return

    hpa_map = hpa_validation.get("hpa_cell_type_mapping", {})
    overlap_map = overlap.get("harmonization_map", {})
    hpa_keys_norm = {_normalize_cfg_label(k) for k in hpa_map.keys()}
    overlap_keys_norm = {_normalize_cfg_label(k) for k in overlap_map.keys()}

    missing_hpa = sorted(set(manual_override_labels) - hpa_keys_norm)
    if missing_hpa:
        raise ValueError(
            "Config inconsistency: all hpa_validation.manual_override_cell_types must be present in "
            f"hpa_validation.hpa_cell_type_mapping. Missing: {missing_hpa}"
        )

    missing_overlap = sorted(set(manual_override_labels) - overlap_keys_norm)
    if missing_overlap:
        raise ValueError(
            "Config inconsistency: all hpa_validation.manual_override_cell_types must be present in "
            f"overlap.harmonization_map to avoid silent Phase 3 drop-through. Missing: {missing_overlap}"
        )


def _validate_marker_extraction_optionals(marker_extraction: dict[str, Any]) -> None:
    min_cells = marker_extraction.get("min_cells_per_cell_type_for_marker_calling")
    if min_cells is not None:
        _require_positive_int(min_cells, "marker_extraction.min_cells_per_cell_type_for_marker_calling")

    min_pct_cells = marker_extraction.get("min_pct_cells")
    if min_pct_cells is not None:
        if not isinstance(min_pct_cells, (int, float)) or not (0.0 <= float(min_pct_cells) <= 1.0):
            raise ValueError(
                "Config key 'marker_extraction.min_pct_cells' must be numeric in [0, 1]."
            )

    max_pct_rest = marker_extraction.get("max_pct_rest")
    if max_pct_rest is not None:
        if not isinstance(max_pct_rest, (int, float)) or not (0.0 <= float(max_pct_rest) <= 1.0):
            raise ValueError(
                "Config key 'marker_extraction.max_pct_rest' must be numeric in [0, 1]."
            )

    stress_genes = marker_extraction.get("adult_stress_signature_genes")
    if stress_genes is not None and not isinstance(stress_genes, list):
        raise ValueError("Config key 'marker_extraction.adult_stress_signature_genes' must be a list when provided.")
    if isinstance(stress_genes, list) and any(not isinstance(g, str) for g in stress_genes):
        raise ValueError("Config key 'marker_extraction.adult_stress_signature_genes' must contain only strings.")

    stress_mode = str(marker_extraction.get("adult_stress_mitigation_mode", "none")).strip().lower()
    if stress_mode not in {"none", "exclude_high_stress", "regress_out"}:
        raise ValueError(
            "Config key 'marker_extraction.adult_stress_mitigation_mode' must be one of: "
            "none, exclude_high_stress, regress_out."
        )

    stress_q = marker_extraction.get("adult_stress_exclude_quantile", 0.99)
    if not isinstance(stress_q, (int, float)) or not (0.0 < float(stress_q) < 1.0):
        raise ValueError(
            "Config key 'marker_extraction.adult_stress_exclude_quantile' must be numeric in (0, 1)."
        )
    
    rescue_score = marker_extraction.get("fetal_rescue_min_score", 0.0)
    if not isinstance(rescue_score, (int, float)):
        raise ValueError(
            "Config key 'marker_extraction.fetal_rescue_min_score' must be numeric."
        )
    
    rescue_margin = marker_extraction.get("fetal_rescue_min_margin", 0.0)
    if not isinstance(rescue_margin, (int, float)) or float(rescue_margin) < 0.0:
        raise ValueError(
            "Config key 'marker_extraction.fetal_rescue_min_margin' must be numeric and >= 0."
        )

    rescue_strategy = str(marker_extraction.get("fetal_rescue_ambiguous_strategy", "keep_original")).strip().lower()
    if rescue_strategy not in {"keep_original", "assign_label"}:
        raise ValueError(
            "Config key 'marker_extraction.fetal_rescue_ambiguous_strategy' must be one of: "
            "keep_original, assign_label."
        )

    rescue_label = marker_extraction.get("fetal_rescue_ambiguous_label", "ambiguous_tubule")
    if not isinstance(rescue_label, str) or not rescue_label.strip():
        raise ValueError(
            "Config key 'marker_extraction.fetal_rescue_ambiguous_label' must be a non-empty string."
        )


def validate_config_schema(cfg: dict[str, Any]) -> None:
    if not isinstance(cfg, dict):
        raise ValueError("Config root must be a mapping/dict.")

    paths = _require_dict(cfg, "paths")
    preprocessing = _require_dict(cfg, "preprocessing")
    marker_extraction = _require_dict(cfg, "marker_extraction")
    hpa_validation = _require_dict(cfg, "hpa_validation")
    overlap = _require_dict(cfg, "overlap")
    _require_dict(cfg, "umap")

    _require_keys(paths, "paths", ["adult_h5ad", "fetal_h5ad", "hpa_tsv", "output_dir"])
    _require_keys(
        preprocessing,
        "preprocessing",
        [
            "min_genes_per_cell",
            "min_cells_per_gene",
            "n_highly_variable_genes",
            "noise_gene_prefixes",
            "noise_gene_exact",
            "apply_min_celltype_filter",
            "min_cells_per_type",
            "adult_labels_to_remove",
            "fetal_labels_to_remove",
        ],
    )
    _require_keys(
        marker_extraction,
        "marker_extraction",
        [
            "top_n_markers",
            "pval_epsilon",
            "pval_adj_cutoff",
            "require_positive_logfc",
            "fetal_generic_kidney_labels",
            "fetal_rescue_distal_markers",
            "fetal_rescue_loh_markers",
            "fetal_rescue_label_distal",
            "fetal_rescue_label_loh",
        ],
    )
    _require_keys(
        hpa_validation,
        "hpa_validation",
        [
            "validation_ncpm_threshold",
            "top_n_validated",
            "pval_threshold",
            "logfc_threshold",
            "weight_extraction_strength",
            "weight_hpa_strength",
            "weight_hpa_specificity",
            "manual_override_cell_types",
            "manual_validated_markers",
            "hpa_cell_type_mapping",
        ],
    )
    _require_keys(overlap, "overlap", ["harmonization_map"])

    _require_positive_int(marker_extraction["top_n_markers"], "marker_extraction.top_n_markers")
    _require_positive_int(hpa_validation["top_n_validated"], "hpa_validation.top_n_validated")
    _validate_marker_extraction_optionals(marker_extraction)

    n_jobs = preprocessing.get("n_jobs", 1)
    if not isinstance(n_jobs, int) or n_jobs == 0 or n_jobs < -1:
        raise ValueError("Config key 'preprocessing.n_jobs' must be -1 or a positive integer.")

    if not isinstance(overlap["harmonization_map"], dict):
        raise ValueError("Config key 'overlap.harmonization_map' must be a mapping/dict.")

    if "manual_mapping_seed" in overlap:
        _validate_manual_mapping_seed(overlap["manual_mapping_seed"])

    _validate_overlap_acceptance_criteria(overlap.get("acceptance_criteria"))
    _validate_dataset_manifest_config(cfg.get("dataset_manifest"))
    _validate_cross_section_consistency(hpa_validation, overlap)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _new_run_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def write_run_manifest(
    manifest_path: Path,
    run_id: str,
    config_path: Path,
    cfg: dict[str, Any],
    paths: dict[str, Path],
    phases: list[int],
    args: argparse.Namespace,
) -> None:
    manifest = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "cwd": str(Path.cwd()),
        "command": " ".join(sys.argv),
        "python_version": sys.version,
        "selected_phases": phases,
        "dry_run": bool(args.dry_run),
        "verify_only": bool(args.verify),
        "config": {
            "path": str(config_path),
            "sha256": _sha256_file(config_path),
            "seed": cfg.get("seed", 42),
            "top_n_markers": cfg["marker_extraction"]["top_n_markers"],
            "top_n_validated": cfg["hpa_validation"]["top_n_validated"],
        },
        "resolved_paths": {k: str(v) for k, v in paths.items()},
        "dependency_versions": _get_package_versions(),
    }
    save_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", manifest_path, "run manifest")


def _json_sha256(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _path_signature(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.exists():
        return {
            "path": str(resolved),
            "exists": False,
            "size": None,
            "mtime_ns": None,
            "type": None,
        }

    stat = resolved.stat()
    if resolved.is_dir():
        n_children = len(list(resolved.iterdir()))
        path_type = "dir"
    else:
        n_children = None
        path_type = "file"

    return {
        "path": str(resolved),
        "exists": True,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "type": path_type,
        "n_children": n_children,
    }


def _is_artifact_complete(path: Path) -> bool:
    if not path.exists():
        return False
    if path.is_dir():
        return any(path.iterdir())
    return path.stat().st_size > 0


def _missing_artifacts(artifacts: dict[str, Path]) -> list[str]:
    return [name for name, path in artifacts.items() if not _is_artifact_complete(path)]


def _has_any_existing_outputs(artifacts: dict[str, Path]) -> bool:
    return any(path.exists() for path in artifacts.values())


def _new_resume_meta() -> dict[str, Any]:
    return {
        "meta_version": 1,
        "phase_state": {},
        "updated_utc": None,
    }


def _load_resume_meta(meta_path: Path) -> dict[str, Any]:
    if not meta_path.exists():
        return _new_resume_meta()

    try:
        data = json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception as e:
        log.warning("Resume metadata unreadable (%s). Starting with empty resume state.", e)
        return _new_resume_meta()

    if not isinstance(data, dict):
        log.warning("Resume metadata malformed (root is not an object). Resetting state.")
        return _new_resume_meta()

    if data.get("meta_version") != 1:
        log.warning(
            "Resume metadata version mismatch (got=%r, expected=1). Resetting state.",
            data.get("meta_version"),
        )
        return _new_resume_meta()

    if not isinstance(data.get("phase_state"), dict):
        log.warning("Resume metadata malformed (missing phase_state map). Resetting state.")
        return _new_resume_meta()

    return data


def _save_resume_meta(meta_path: Path, meta: dict[str, Any]) -> None:
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = meta_path.with_name(meta_path.name + ".tmp")
    tmp_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp_path.replace(meta_path)


def _record_phase_resume_state(
    meta: dict[str, Any],
    phase_id: int,
    fingerprint: str,
    artifacts: dict[str, Path],
    run_id: str,
) -> None:
    phase_map = meta.setdefault("phase_state", {})
    phase_map[str(phase_id)] = {
        "fingerprint": fingerprint,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "run_id": run_id,
        "artifacts": {name: str(path.resolve()) for name, path in artifacts.items()},
    }
    meta["updated_utc"] = datetime.now(timezone.utc).isoformat()


def _can_resume_phase(
    phase_id: int,
    fingerprint: str,
    expected_artifacts: dict[str, Path],
    meta: dict[str, Any],
) -> tuple[bool, str]:
    phase_state = meta.get("phase_state", {}).get(str(phase_id))
    if not isinstance(phase_state, dict):
        return False, "no previous resume metadata"

    previous = phase_state.get("fingerprint")
    if previous != fingerprint:
        return False, "related configuration/code fingerprint changed"

    missing = _missing_artifacts(expected_artifacts)
    if missing:
        return False, f"missing/empty artifacts: {missing}"

    return True, "fingerprint unchanged and artifacts are complete"


def _ensure_upstream_compatibility(
    meta: dict[str, Any],
    upstream_phase_id: int,
    expected_fingerprint: str,
    current_phase_label: str,
) -> None:
    phase_state = meta.get("phase_state", {}).get(str(upstream_phase_id))
    if not isinstance(phase_state, dict):
        log.warning(
            "%s: no resume metadata found for upstream Phase %d; cannot verify compatibility.",
            current_phase_label,
            upstream_phase_id,
        )
        return

    recorded = phase_state.get("fingerprint")
    if recorded != expected_fingerprint:
        raise RuntimeError(
            f"{current_phase_label} depends on upstream Phase {upstream_phase_id}, but related parameters "
            "changed since that upstream phase was last completed. "
            f"Re-run including Phase {upstream_phase_id} to regenerate compatible artifacts."
        )


def _build_phase1_fingerprint(cfg: dict[str, Any], paths: dict[str, Path], script_sha256: str) -> str:
    payload = {
        "phase": 1,
        "script_sha256": script_sha256,
        "seed": cfg.get("seed", 42),
        "preprocessing": cfg["preprocessing"],
        "marker_extraction": cfg["marker_extraction"],
        "fetal_main_output_top_n": cfg["hpa_validation"]["top_n_validated"],
        "umap": cfg["umap"],
        "inputs": {
            "adult_h5ad": _path_signature(paths["adult_h5ad"]),
            "fetal_h5ad": _path_signature(paths["fetal_h5ad"]),
        },
    }
    return _json_sha256(payload)


def _build_phase2_fingerprint(
    cfg: dict[str, Any],
    paths: dict[str, Path],
    phase1_fingerprint: str,
    script_sha256: str,
) -> str:
    payload = {
        "phase": 2,
        "script_sha256": script_sha256,
        "hpa_validation": cfg["hpa_validation"],
        "marker_extraction_pval_epsilon": cfg["marker_extraction"].get("pval_epsilon", 1e-300),
        "inputs": {
            "hpa_tsv": _path_signature(paths["hpa_tsv"]),
        },
        "depends_on_phase1": phase1_fingerprint,
    }
    return _json_sha256(payload)


def _build_phase3_fingerprint(
    cfg: dict[str, Any],
    phase1_fingerprint: str,
    phase2_fingerprint: str,
    script_sha256: str,
) -> str:
    payload = {
        "phase": 3,
        "script_sha256": script_sha256,
        "overlap": cfg["overlap"],
        "hpa_top_n_validated": cfg["hpa_validation"]["top_n_validated"],
        "umap_fig_dpi": cfg.get("umap", {}).get("fig_dpi", 300),
        "depends_on_phase1": phase1_fingerprint,
        "depends_on_phase2": phase2_fingerprint,
    }
    return _json_sha256(payload)


def load_config(config_path: Path) -> dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with config_path.open() as fh:
        cfg = yaml.safe_load(fh)
    validate_config_schema(cfg)
    log.info("Loaded config: %s", config_path)
    return cfg


def _normalize_token(value: str) -> str:
    """Normalize file-name tokens for tolerant matching (spaces, unicode spaces, punctuation)."""
    return re.sub(r"[^a-z0-9]+", "", str(value).casefold())


def _resolve_output_path(path_value: str, config_dir: Path) -> Path:
    p = Path(path_value)
    return p.resolve() if p.is_absolute() else (config_dir / p).resolve()


def _resolve_input_path(path_value: str, config_dir: Path) -> Path:
    """
    Resolve input paths robustly.

    Strategy:
    1) Absolute path as-is.
    2) Config-relative path.
    3) Ancestor-relative path (for repo moves).
    4) Filename-normalized fallback inside candidate parent dirs
       (handles cases like Fetal_17962.h5ad vs Fetal_17\u202f962.h5ad).
    """
    raw = Path(path_value)

    if raw.is_absolute():
        return raw.resolve()

    primary = (config_dir / raw).resolve()
    if primary.exists():
        return primary

    for anchor in config_dir.parents:
        candidate = (anchor / raw).resolve()
        if candidate.exists():
            log.warning("Resolved input via ancestor path: %s -> %s", path_value, candidate)
            return candidate

    normalized_target = _normalize_token(raw.name)
    search_dirs: list[Path] = []
    for anchor in [config_dir, *config_dir.parents]:
        candidate_dir = (anchor / raw.parent).resolve()
        if candidate_dir not in search_dirs:
            search_dirs.append(candidate_dir)

    matches: list[Path] = []
    seen: set[Path] = set()

    for directory in search_dirs:
        if not directory.exists() or not directory.is_dir():
            continue
        for child in sorted(directory.iterdir(), key=lambda p: p.name.casefold()):
            if child.is_file() and _normalize_token(child.name) == normalized_target:
                resolved = child.resolve()
                if resolved not in seen:
                    matches.append(resolved)
                    seen.add(resolved)

    if len(matches) == 1:
        log.warning(
            "Resolved input via normalized filename match: requested=%s matched=%s",
            path_value,
            matches[0],
        )
        return matches[0]

    if len(matches) > 1:
        candidates = "\n".join(f"  - {p}" for p in matches)
        raise FileNotFoundError(
            "Ambiguous normalized filename match for input path "
            f"'{path_value}'. Multiple candidates found:\n{candidates}\n"
            "Please set an explicit absolute or config-relative path."
        )

    return primary


def build_paths(cfg: dict, config_dir: Path) -> dict[str, Path]:
    """Resolve all pipeline paths with robust input handling and deterministic output roots."""
    p = cfg["paths"]
    output_dir = _resolve_output_path(p["output_dir"], config_dir)
    aux_dir = output_dir / p.get("aux_subdir", "aux")
    # Keep support figures under aux so the main output folder stays minimal.
    figures_dir = aux_dir / p.get("figures_subdir", "figures")
    diagnostics_dir = aux_dir / p.get("diagnostics_subdir", "diagnostics")

    return {
        "adult_h5ad":   _resolve_input_path(p["adult_h5ad"], config_dir),
        "fetal_h5ad":   _resolve_input_path(p["fetal_h5ad"], config_dir),
        "hpa_tsv":      _resolve_input_path(p["hpa_tsv"], config_dir),
        "output_dir":   output_dir,
        "aux_dir":      aux_dir.resolve(),
        "figures_dir":  figures_dir.resolve(),
        "diagnostics_dir": diagnostics_dir.resolve(),
    }


# =============================================================================
# SECTION 2 — Shared Utilities
# =============================================================================

def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input not found [{label}]: {path}")


def require_columns(df: PandasDataFrame, required: list[str], source: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Schema error [{source}]: missing columns {missing}. Present: {list(df.columns)}")


def save_tsv(df: PandasDataFrame, path: Path, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    log.info("Saved [%s]: %s  (%d rows)", label, path, len(df))


def save_text(text: str, path: Path, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    log.info("Saved [%s]: %s", label, path)


def save_figure(fig: Any, path: Path, label: str, dpi: int = 300) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved [%s]: %s", label, path)


def gate_phase(artifacts: dict[str, Path], phase_name: str) -> None:
    """Fail fast if any upstream artifact is missing or empty."""
    missing = [f"  [{k}]  {v}" for k, v in artifacts.items() if not v.exists()]
    empty   = [f"  [{k}]  {v}" for k, v in artifacts.items() if v.exists() and v.stat().st_size == 0]
    errors = []
    if missing:
        errors.append("Missing upstream artifacts required by " + phase_name + ":\n" + "\n".join(missing))
    if empty:
        errors.append("Empty upstream artifacts required by " + phase_name + ":\n" + "\n".join(empty))
    if errors:
        raise RuntimeError("\n\n".join(errors) + f"\n\nRun all preceding phases before {phase_name}.")
    log.info("Gate passed for %s — %d upstream artifacts present.", phase_name, len(artifacts))


def require_runtime_packages(required: dict[str, Any], context: str) -> None:
    missing = [name for name, module in required.items() if module is None]
    if missing:
        raise ModuleNotFoundError(
            f"{context} requires missing Python package(s): {', '.join(sorted(missing))}. "
            "Install them in the active environment and retry."
        )


def build_harmonization_lookup(raw_map: dict[str, str]) -> dict[str, str]:
    def _norm(name: str) -> str:
        t = str(name).strip().lower()
        t = re.sub(r"[_\-]+", " ", t)
        return re.sub(r"\s+", " ", t).strip()
    return {_norm(k): v for k, v in raw_map.items()}


def harmonize_cell_type(name: str, lookup: dict[str, str]) -> str:
    def _norm(n: str) -> str:
        t = str(n).strip().lower()
        t = re.sub(r"[_\-]+", " ", t)
        return re.sub(r"\s+", " ", t).strip()
    return lookup.get(_norm(name), _norm(name))


def build_top_markers(
    markers: PandasDataFrame,
    top_n: int,
    pval_eps: float = 1e-300,
    min_pct_cells: float = 0.0,
    max_pct_rest: float = 1.0,
) -> PandasDataFrame:
    """Rank markers within each cell type by combined p-value + logFC rank score."""
    require_columns(markers, ["cell_type", "gene", "logfoldchange", "pval_adj", "pct_cells", "pct_rest"], "markers")
    top = markers.copy()

    before_quality = len(top)
    if min_pct_cells > 0.0:
        top = top[top["pct_cells"] >= float(min_pct_cells)].copy()
    if max_pct_rest < 1.0:
        top = top[top["pct_rest"] <= float(max_pct_rest)].copy()
    removed_quality = before_quality - len(top)
    if removed_quality > 0:
        log.info(
            "Marker quality gates removed %d rows (min_pct_cells>=%.3f, max_pct_rest<=%.3f).",
            removed_quality,
            min_pct_cells,
            max_pct_rest,
        )

    top["pval_adj_safe"] = top["pval_adj"].clip(lower=pval_eps)
    top["p_rank"]   = top.groupby("cell_type")["pval_adj_safe"].rank(method="min", ascending=True)
    top["lfc_rank"] = top.groupby("cell_type")["logfoldchange"].rank(method="min", ascending=False)
    top["combined_rank_score"] = top["p_rank"] + top["lfc_rank"]
    top = top.sort_values(
        ["cell_type", "combined_rank_score", "pval_adj_safe", "logfoldchange"],
        ascending=[True, True, True, False],
    ).reset_index(drop=True)
    top["rank_top"] = top.groupby("cell_type").cumcount() + 1
    return top[top["rank_top"] <= top_n].copy()


def resolve_marker_calling_min_cells(cfg: dict[str, Any]) -> int:
    """
    Resolve minimum cells per cell type required for marker calling acceptance.
    """
    marker_cfg = cfg["marker_extraction"]
    pre_cfg = cfg["preprocessing"]
    explicit = marker_cfg.get("min_cells_per_cell_type_for_marker_calling")
    if explicit is not None:
        return int(explicit)
    if pre_cfg.get("apply_min_celltype_filter", False):
        return int(pre_cfg["min_cells_per_type"])
    return 1


def enforce_marker_calling_min_cells(adata: Any, min_cells: int, label: str) -> None:
    """Acceptance gate: every cell type used for marker calling must meet minimum size."""
    counts = adata.obs["cell_type"].astype(str).value_counts()
    violating = counts[counts < min_cells]
    if len(violating) > 0:
        details = ", ".join(f"{ct}:{int(n)}" for ct, n in violating.items())
        raise ValueError(
            f"[{label}] Marker-calling acceptance failed: {len(violating)} cell type(s) below "
            f"minimum cells ({min_cells}). Offenders: {details}"
        )
    log.info(
        "[%s] Marker-calling acceptance passed: all %d cell types have >= %d cells.",
        label,
        len(counts),
        min_cells,
    )


def exclude_adult_stress_signature_markers(markers_adult: PandasDataFrame, cfg_me: dict[str, Any]) -> PandasDataFrame:
    """Optional mitigation for stress-signature confounding in adult marker tables."""
    stress_genes_raw = cfg_me.get("adult_stress_signature_genes", [])
    stress_genes = {str(g).strip().upper() for g in stress_genes_raw if str(g).strip()}
    if not stress_genes:
        return markers_adult

    out = markers_adult.copy()
    out["gene_upper"] = out["gene"].astype(str).str.upper()
    before = len(out)
    out = out[~out["gene_upper"].isin(stress_genes)].copy()
    removed = before - len(out)
    out = out.drop(columns=["gene_upper"])

    log.info(
        "[adult] Stress-signature filter removed %d marker rows using %d configured genes.",
        removed,
        len(stress_genes),
    )
    if out.empty:
        raise ValueError(
            "[adult] Stress-signature filtering removed all adult marker rows. "
            "Adjust marker_extraction.adult_stress_signature_genes."
        )
    return out


def write_dataset_manifest(cfg: dict[str, Any], paths: dict[str, Path]) -> Path | None:
    """Write dataset_manifest.tsv with provenance fields for adult/fetal references."""
    if pd is None:
        log.warning("Skipping dataset_manifest.tsv: pandas is unavailable.")
        return None

    manifest_cfg = cfg.get("dataset_manifest", {})
    if not isinstance(manifest_cfg, dict):
        manifest_cfg = {}

    rows: list[dict[str, Any]] = []
    for stage, path_key in [("adult", "adult_h5ad"), ("fetal", "fetal_h5ad")]:
        stage_cfg = manifest_cfg.get(stage, {})
        if not isinstance(stage_cfg, dict):
            stage_cfg = {}

        input_path = paths[path_key]
        row = {
            "stage": stage,
            "source": stage_cfg.get("source", "unspecified"),
            "accession": stage_cfg.get("accession", "unspecified"),
            "tissue_region": stage_cfg.get("tissue_region", "kidney"),
            "donor_count": stage_cfg.get("donor_count", "unspecified"),
            "platform": stage_cfg.get("platform", "unspecified"),
            "download_date": stage_cfg.get("download_date", "unspecified"),
            "input_path": str(input_path),
            "input_exists": bool(input_path.exists()),
        }
        rows.append(row)

    out_df = pd.DataFrame(rows)
    out_path = paths["aux_dir"] / "dataset_manifest.tsv"
    save_tsv(out_df, out_path, "dataset manifest")
    return out_path


def write_marker_method_params(cfg: dict[str, Any], paths: dict[str, Path]) -> Path:
    """Write marker_method_params.yaml to make statistical settings explicit and reproducible."""
    payload = {
        "schema_version": 1,
        "seed": cfg.get("seed", 42),
        "preprocessing": {
            "min_genes_per_cell": cfg["preprocessing"]["min_genes_per_cell"],
            "min_cells_per_gene": cfg["preprocessing"]["min_cells_per_gene"],
            "n_highly_variable_genes": cfg["preprocessing"]["n_highly_variable_genes"],
            "apply_min_celltype_filter": cfg["preprocessing"]["apply_min_celltype_filter"],
            "min_cells_per_type": cfg["preprocessing"]["min_cells_per_type"],
        },
        "marker_extraction": {
            "method": "wilcoxon",
            "corr_method": "benjamini-hochberg",
            "use_raw": True,
            "top_n_markers": cfg["marker_extraction"]["top_n_markers"],
            "min_cells_per_cell_type_for_marker_calling": resolve_marker_calling_min_cells(cfg),
            "min_pct_cells": cfg["marker_extraction"].get("min_pct_cells", 0.0),
            "max_pct_rest": cfg["marker_extraction"].get("max_pct_rest", 1.0),
            "pval_adj_cutoff": cfg["marker_extraction"]["pval_adj_cutoff"],
            "require_positive_logfc": cfg["marker_extraction"]["require_positive_logfc"],
            "pval_epsilon": cfg["marker_extraction"]["pval_epsilon"],
            "adult_stress_signature_genes": cfg["marker_extraction"].get("adult_stress_signature_genes", []),
            "adult_stress_mitigation_mode": cfg["marker_extraction"].get("adult_stress_mitigation_mode", "none"),
            "adult_stress_exclude_quantile": cfg["marker_extraction"].get("adult_stress_exclude_quantile", 0.99),
            "fetal_rescue_min_score": cfg["marker_extraction"].get("fetal_rescue_min_score", 0.0),
            "fetal_rescue_min_margin": cfg["marker_extraction"].get("fetal_rescue_min_margin", 0.0),
            "fetal_rescue_ambiguous_strategy": cfg["marker_extraction"].get("fetal_rescue_ambiguous_strategy", "keep_original"),
            "fetal_rescue_ambiguous_label": cfg["marker_extraction"].get("fetal_rescue_ambiguous_label", "ambiguous_tubule"),
        },
        "hpa_validation": {
            "validation_ncpm_threshold": cfg["hpa_validation"]["validation_ncpm_threshold"],
            "top_n_validated": cfg["hpa_validation"]["top_n_validated"],
            "pval_threshold": cfg["hpa_validation"]["pval_threshold"],
            "logfc_threshold": cfg["hpa_validation"]["logfc_threshold"],
            "weight_extraction_strength": cfg["hpa_validation"]["weight_extraction_strength"],
            "weight_hpa_strength": cfg["hpa_validation"]["weight_hpa_strength"],
            "weight_hpa_specificity": cfg["hpa_validation"]["weight_hpa_specificity"],
        },
        "overlap_yaml_policy": {
            "adult_top_n_per_harmonized_cell_type": 25,
            "fetal_top_n_per_harmonized_cell_type": None,
            "gene_ordering": [
                "rank (if available)",
                "final_joint_score (desc)",
                "combined_rank_score (asc)",
                "pval_adj (asc)",
                "logfoldchange (desc)",
            ],
            "acceptance_criteria": cfg["overlap"].get("acceptance_criteria", {}),
        },
    }

    out_path = paths["aux_dir"] / "marker_method_params.yaml"
    text = yaml.safe_dump(payload, sort_keys=False)
    save_text(text, out_path, "marker method params")
    return out_path


# =============================================================================
# SECTION 3 — Phase 1: Preprocess + Marker Extraction + Fetal Rescue
# =============================================================================

def _ratio_ensembl_like(values: list[str]) -> float:
    if not values:
        return 0.0
    pat = re.compile(r"^ENS[A-Z0-9]*\d+(?:\.\d+)?$", flags=re.IGNORECASE)
    n = sum(1 for v in values if pat.match(str(v).strip()))
    return n / len(values)


def _ratio_numeric_like(values: list[str]) -> float:
    if not values:
        return 0.0
    n = sum(1 for v in values if str(v).strip().isdigit())
    return n / len(values)


def _ensure_gene_symbols(adata, label: str):
    current_names = adata.var_names.astype(str).tolist()
    source_used = "var_names"

    if "feature_name" in adata.var.columns:
        feature_names = adata.var["feature_name"].astype(str).tolist()
        if len(feature_names) != len(current_names):
            raise ValueError(
                f"[{label}] feature_name length mismatch with var_names "
                f"({len(feature_names)} vs {len(current_names)})."
            )

        if feature_names != current_names:
            current_ens = _ratio_ensembl_like(current_names)
            feature_ens = _ratio_ensembl_like(feature_names)
            current_num = _ratio_numeric_like(current_names)
            feature_num = _ratio_numeric_like(feature_names)
            feature_invalid_ratio = (
                sum(1 for g in feature_names if str(g).strip() == "" or str(g).strip().lower() == "nan")
                / len(feature_names)
            )

            should_switch = (
                (current_ens >= 0.70 and feature_ens <= 0.30)
                or (current_num >= 0.80 and feature_num <= 0.20)
            )

            if feature_invalid_ratio > 0.01:
                should_switch = False
                log.warning(
                    "[%s] Gene naming: kept existing var_names; feature_name has %.2f%% blank/NaN-like values.",
                    label,
                    feature_invalid_ratio * 100.0,
                )

            if should_switch:
                changed_pairs = [
                    f"{old}->{new}"
                    for old, new in zip(current_names, feature_names)
                    if old != new
                ][:5]
                adata.var_names = pd.Index(feature_names)
                source_used = "feature_name"
                log.info(
                    "[%s] Gene naming: switched var_names to feature_name (%d changes; examples: %s).",
                    label,
                    sum(1 for old, new in zip(current_names, feature_names) if old != new),
                    changed_pairs,
                )
            elif current_ens <= 0.30 and feature_ens >= 0.70:
                log.warning(
                    "[%s] Gene naming: kept existing var_names; feature_name appears Ensembl-like and may corrupt marker symbols.",
                    label,
                )
            else:
                log.info(
                    "[%s] Gene naming: kept existing var_names; feature_name differs but no clear symbol-quality gain.",
                    label,
                )

    before_names = adata.var_names.astype(str).tolist()
    dup_counts = pd.Index(before_names).value_counts()
    duplicate_name_count = int((dup_counts > 1).sum())
    duplicate_entry_count = int((dup_counts - 1).clip(lower=0).sum())
    if duplicate_entry_count > 0:
        examples = dup_counts[dup_counts > 1].index.astype(str).tolist()[:8]
        log.warning(
            "[%s] Gene naming: found %d duplicate entries across %d duplicated symbols before make_unique (examples: %s).",
            label,
            duplicate_entry_count,
            duplicate_name_count,
            examples,
        )

    adata.var_names_make_unique()
    after_names = adata.var_names.astype(str).tolist()
    renamed = sum(1 for old, new in zip(before_names, after_names) if old != new)
    if renamed > 0:
        rename_examples = [
            f"{old}->{new}"
            for old, new in zip(before_names, after_names)
            if old != new
        ][:5]
        log.warning(
            "[%s] Gene naming: var_names_make_unique renamed %d entries (examples: %s).",
            label,
            renamed,
            rename_examples,
        )

    log.info("[%s] Gene symbols ready using %s (%d genes).", label, source_used, adata.n_vars)
    return adata


def _apply_adult_stress_mitigation(adata, cfg_me: dict[str, Any], label: str):
    if label != "adult":
        return adata

    mode = str(cfg_me.get("adult_stress_mitigation_mode", "none")).strip().lower()
    if mode == "none":
        return adata

    stress_genes_raw = cfg_me.get("adult_stress_signature_genes", [])
    stress_genes = [str(g).strip() for g in stress_genes_raw if str(g).strip()]
    if not stress_genes:
        log.warning("[adult] Stress mitigation requested but no stress genes configured; skipping.")
        return adata

    var_name_lookup = {str(g).upper(): str(g) for g in adata.var_names.astype(str)}
    available = []
    for g in stress_genes:
        resolved = var_name_lookup.get(g.upper())
        if resolved is not None and resolved not in available:
            available.append(resolved)
    if not available:
        log.warning("[adult] Stress mitigation requested but none of the stress genes are present; skipping.")
        return adata

    import scanpy as sc

    sc.tl.score_genes(adata, gene_list=available, score_name="_stress_score", use_raw=False)

    if mode == "exclude_high_stress":
        q = float(cfg_me.get("adult_stress_exclude_quantile", 0.99))
        threshold = float(adata.obs["_stress_score"].quantile(q))
        keep_mask = adata.obs["_stress_score"] <= threshold
        removed = int((~keep_mask).sum())
        if removed >= int(adata.n_obs):
            raise ValueError(
                "[adult] Stress mitigation would remove all cells. "
                "Relax adult_stress_exclude_quantile or disable mitigation."
            )
        if removed > 0:
            adata = adata[keep_mask].copy()
        log.info(
            "[adult] Stress mitigation: excluded %d high-stress cells at quantile %.3f (threshold=%.4f, genes_used=%d).",
            removed,
            q,
            threshold,
            len(available),
        )
    elif mode == "regress_out":
        sc.pp.regress_out(adata, keys=["_stress_score"])
        log.info(
            "[adult] Stress mitigation: regressed out stress score using %d present stress genes.",
            len(available),
        )
    else:
        raise ValueError(f"Unsupported adult stress mitigation mode: {mode}")

    if "_stress_score" in adata.obs.columns:
        adata.obs.drop(columns=["_stress_score"], inplace=True)
    return adata


def _remove_labels(adata, labels_to_remove: set[str], label: str):
    ct_norm = adata.obs["cell_type"].astype(str).str.strip().str.lower()
    mask = ct_norm.isin(labels_to_remove)
    found = sorted(set(ct_norm.unique()) & labels_to_remove)
    if found:
        log.info("[%s] Removing %d cells with labels: %s", label, int(mask.sum()), found)
    adata = adata[~mask].copy()
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    return adata


def _collapse_adult_endothelial_cell_types(adata):
    """Collapse adult endothelial subtypes into a single label for marker calling."""
    ct = adata.obs["cell_type"].astype(str).str.strip()
    ct_lc = ct.str.lower()
    endothelial_mask = ct_lc.str.contains(r"\bendothelial\b", regex=True, na=False)
    if int(endothelial_mask.sum()) == 0:
        return adata

    original_labels = sorted(ct[endothelial_mask].unique().tolist())
    ct.loc[endothelial_mask] = "endothelial cells"
    adata.obs["cell_type"] = ct.astype("category")
    log.info(
        "[adult] Collapsed %d endothelial cells from %d labels into 'endothelial cells': %s",
        int(endothelial_mask.sum()),
        len(original_labels),
        original_labels,
    )
    return adata


def _drop_fetal_post_rescue_generic_labels(adata_fetal, labels_to_drop: set[str]):
    """Drop residual generic fetal labels after rescue assignment."""
    ct_norm = adata_fetal.obs["cell_type"].astype(str).str.strip().str.lower()
    mask = ct_norm.isin(labels_to_drop)
    found = sorted(set(ct_norm.unique()) & labels_to_drop)
    removed = int(mask.sum())
    if removed > 0:
        log.info(
            "[fetal] Post-rescue cleanup: dropping %d cells with labels: %s",
            removed,
            found,
        )
        adata_fetal = adata_fetal[~mask].copy()
    else:
        log.info("[fetal] Post-rescue cleanup: no residual generic labels to drop.")

    adata_fetal.obs["cell_type"] = adata_fetal.obs["cell_type"].astype("category")
    return adata_fetal


def _filter_noise_genes(adata, noise_prefixes: tuple, noise_exact: set, label: str):
    mask_x = ~(adata.var_names.str.startswith(noise_prefixes) | adata.var_names.isin(noise_exact))
    before_x = adata.n_vars
    adata = adata[:, mask_x].copy()
    if adata.raw is not None:
        raw_ad = adata.raw.to_adata()
        mask_raw = ~(raw_ad.var_names.str.startswith(noise_prefixes) | raw_ad.var_names.isin(noise_exact))
        before_raw = raw_ad.n_vars
        adata.raw = raw_ad[:, mask_raw]
        log.info("[%s] Noise filter: .X %d→%d, .raw %d→%d", label, before_x, adata.n_vars, before_raw, adata.raw.n_vars)
    else:
        log.info("[%s] Noise filter: .X %d→%d", label, before_x, adata.n_vars)
    return adata


def _preprocess_reference(adata, label: str, cfg_pre: dict, cfg_me: dict, seed: int):
    import scanpy as sc
    adata = adata.copy()
    log.info("[%s] Preprocessing: %d cells × %d genes", label, adata.n_obs, adata.n_vars)
    sc.pp.filter_cells(adata, min_genes=cfg_pre["min_genes_per_cell"])
    sc.pp.filter_genes(adata, min_cells=cfg_pre["min_cells_per_gene"])
    log.info("[%s] After QC: %d cells × %d genes", label, adata.n_obs, adata.n_vars)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata = _apply_adult_stress_mitigation(adata, cfg_me, label)

    adata.raw = adata
    batch_key = "dataset_id" if "dataset_id" in adata.obs.columns else None
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg_pre["n_highly_variable_genes"], batch_key=batch_key)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50, random_state=seed)
    n_pcs = min(30, adata.obsm["X_pca"].shape[1])
    if n_pcs < 2:
        raise ValueError(f"[{label}] Not enough PCs ({n_pcs}) to build neighbour graph.")
    sc.pp.neighbors(adata, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    log.info("[%s] Preprocessing done: %d cells × %d HVGs", label, adata.n_obs, adata.n_vars)
    return adata


def _enforce_min_cells_per_type(adata, min_cells: int, label: str):
    counts = adata.obs["cell_type"].value_counts()
    keep = counts[counts >= min_cells].index
    removed = counts[counts < min_cells]
    if len(removed):
        log.info("[%s] Dropped %d low-count cell types (<%d cells): %s", label, len(removed), min_cells, removed.index.tolist())
    filtered = adata[adata.obs["cell_type"].isin(keep)].copy()
    filtered.obs["cell_type"] = filtered.obs["cell_type"].astype("category")
    return filtered


def _plot_umap_two_panel(adata_adult, adata_fetal, out_path: Path, cfg_umap: dict) -> None:
    import scanpy as sc
    import matplotlib
    import matplotlib.patches as mpatches
    import matplotlib.patheffects as pe

    def _get_palette(adata, key):
        cats = adata.obs[key].cat.categories.tolist()
        stored = adata.uns.get(f"{key}_colors")
        if stored is not None and len(stored) >= len(cats):
            return dict(zip(cats, stored))
        cmap = matplotlib.colormaps.get_cmap("tab20")
        return {c: cmap(i / max(len(cats) - 1, 1)) for i, c in enumerate(cats)}

    def _compute_cluster_anchors(adata, key):
        coords = adata.obsm["X_umap"]
        anchors = {}
        for cat in adata.obs[key].cat.categories:
            mask = np.asarray(adata.obs[key] == cat)
            pts = coords[mask]
            if pts.shape[0] == 0:
                continue
            med = np.median(pts, axis=0)
            d2 = np.sum((pts - med) ** 2, axis=1)
            anchors[cat] = pts[np.argmin(d2)]
        return anchors

    def _repel_points(points, xlim, ylim, min_dist=0.06, n_iter=220, pull=0.028, max_shift=0.07):
        if len(points) <= 1:
            return points.copy()

        x0, x1 = xlim
        y0, y1 = ylim
        xr = max(abs(x1 - x0), 1e-9)
        yr = max(abs(y1 - y0), 1e-9)

        p = points.copy()
        p[:, 0] = (p[:, 0] - x0) / xr
        p[:, 1] = (p[:, 1] - y0) / yr
        a = p.copy()

        for _ in range(n_iter):
            p += (a - p) * pull

            for i in range(len(p)):
                for j in range(i + 1, len(p)):
                    d = p[j] - p[i]
                    dist = np.hypot(d[0], d[1]) + 1e-9
                    if dist < min_dist:
                        step = 0.5 * (min_dist - dist) / dist
                        shift = d * step
                        p[i] -= shift
                        p[j] += shift

            delta = p - a
            delta_norm = np.linalg.norm(delta, axis=1, keepdims=True) + 1e-9
            too_far = delta_norm[:, 0] > max_shift
            if np.any(too_far):
                p[too_far] = a[too_far] + delta[too_far] * (max_shift / delta_norm[too_far])

            p[:, 0] = np.clip(p[:, 0], 0.03, 0.97)
            p[:, 1] = np.clip(p[:, 1], 0.03, 0.97)

        out = p.copy()
        out[:, 0] = out[:, 0] * xr + x0
        out[:, 1] = out[:, 1] * yr + y0
        return out

    def _plot_umap_panel(
        adata,
        ax,
        title,
        color_key="cell_type",
        tissue_key="tissue",
        panel_label=None,
        label_min_dist=0.065,
        label_max_shift=0.07,
    ):
        draw_tissue_bg = False
        if tissue_key in adata.obs.columns:
            n_tissues = adata.obs[tissue_key].astype("category").cat.categories.size
            draw_tissue_bg = n_tissues > 1

        if draw_tissue_bg:
            sc.pl.umap(
                adata,
                color=tissue_key,
                ax=ax,
                show=False,
                size=12,
                alpha=0.28,
                legend_loc="none",
                title="",
                frameon=False,
            )

        sc.pl.umap(
            adata,
            color=color_key,
            ax=ax,
            show=False,
            size=3,
            alpha=0.86,
            legend_loc="none",
            title="",
            frameon=False,
        )

        palette = _get_palette(adata, color_key)
        cell_types = adata.obs[color_key].cat.categories.tolist()
        anchors = _compute_cluster_anchors(adata, color_key)

        present = [ct for ct in cell_types if ct in anchors]
        anchor_pts = np.array([anchors[ct] for ct in present])

        moved_pts = _repel_points(
            anchor_pts,
            ax.get_xlim(),
            ax.get_ylim(),
            min_dist=label_min_dist,
            max_shift=label_max_shift,
        )
        pt_map = {ct: moved_pts[i] for i, ct in enumerate(present)}

        for idx, ct in enumerate(cell_types, start=1):
            if ct not in anchors:
                continue
            cx, cy = anchors[ct]
            tx, ty = pt_map[ct]

            ax.plot([cx, tx], [cy, ty], color="#4d4d4d", linewidth=0.45, alpha=0.74, zorder=9)
            ax.scatter([cx], [cy], s=8, c="#1f1f1f", alpha=0.75, linewidths=0, zorder=9)

            txt = ax.text(
                tx,
                ty,
                str(idx),
                fontsize=6.8,
                fontweight="bold",
                ha="center",
                va="center",
                color="#111111",
                zorder=10,
                bbox=dict(boxstyle="square,pad=0.18", fc="white", ec="#202020", lw=0.7),
            )
            txt.set_path_effects([pe.withStroke(linewidth=1.8, foreground="white")])

        ax.set_title(title, fontsize=9.8, fontweight="bold", pad=5)
        ax.set_xlabel("UMAP1", fontsize=8)
        ax.set_ylabel("UMAP2", fontsize=8)
        ax.tick_params(labelsize=7)
        for spine in ax.spines.values():
            spine.set_visible(False)

        if panel_label:
            ax.text(-0.08, 1.04, panel_label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")

        return cell_types, palette

    def _build_legend(
        fig,
        ax_ref,
        adata,
        color_key,
        tissue_key,
        cell_types,
        palette,
        base_fontsize=6.0,
        legend_width=0.15,
        gap=0.02,
    ):
        bbox = ax_ref.get_position()
        left = bbox.x1 + gap
        width = legend_width
        leg_ax = fig.add_axes([left, bbox.y0, width, bbox.height], facecolor="white", zorder=30)
        leg_ax.patch.set_alpha(1.0)
        leg_ax.axis("off")
        leg_ax.set_xlim(0, 1)
        leg_ax.set_ylim(0, 1)

        t = leg_ax.transAxes
        y = 0.99

        if tissue_key in adata.obs.columns:
            tissues = adata.obs[tissue_key].cat.categories.tolist()
            t_colors = adata.uns.get(
                f"{tissue_key}_colors",
                [matplotlib.colormaps["Set2"](i) for i in range(len(tissues))],
            )

            leg_ax.text(0.00, y, "Tissue", fontsize=base_fontsize + 1.4, fontweight="bold", va="top", transform=t, color="#202020")
            y -= 0.045
            leg_ax.plot([0, 0.97], [y, y], color="#cfcfcf", linewidth=0.6, transform=t, clip_on=False)
            y -= 0.020

            row_h = 0.042
            for i, (tissue, color) in enumerate(zip(tissues, t_colors)):
                leg_ax.add_patch(mpatches.Rectangle((0.00, y - row_h + 0.004), 0.94, row_h - 0.006, facecolor="#f8f8f8" if i % 2 == 0 else "white", edgecolor="none", transform=t, clip_on=False))
                leg_ax.add_patch(mpatches.Rectangle((0.012, y - 0.026), 0.022, 0.022, facecolor=color, edgecolor="#d2d2d2", linewidth=0.5, transform=t, clip_on=False))
                leg_ax.text(0.045, y - 0.015, tissue, fontsize=base_fontsize, va="center", transform=t, color="#303030")
                y -= row_h

        y -= 0.010
        leg_ax.text(0.00, y, "Cell type", fontsize=base_fontsize + 1.4, fontweight="bold", va="top", transform=t, color="#202020")
        y -= 0.045
        leg_ax.plot([0, 0.97], [y, y], color="#cfcfcf", linewidth=0.6, transform=t, clip_on=False)
        y -= 0.018

        n_items = len(cell_types)
        available = max(y + 0.005, 0.10)

        two_col = False
        n_rows = n_items
        row_h = available / max(n_rows, 1)
        if row_h < 0.029 and n_items > 10:
            two_col = True
            n_rows = int(np.ceil(n_items / 2.0))
            row_h = available / max(n_rows, 1)

        row_h = min(0.038, max(0.024, row_h))
        col_w = 0.50 if two_col else 1.0
        max_char = 23 if two_col else 46
        effective_font = base_fontsize if row_h >= 0.030 else max(5.4, base_fontsize - 0.5)

        for i, ct in enumerate(cell_types, start=1):
            col = (i - 1) % 2 if two_col else 0
            row = (i - 1) // 2 if two_col else (i - 1)
            xoff = col * col_w
            ypos = y - row * row_h

            if ypos < -0.01:
                continue

            leg_ax.add_patch(mpatches.Rectangle((xoff + 0.00, ypos - row_h + 0.003), col_w - 0.03, row_h - 0.004, facecolor="#f8f8f8" if row % 2 == 0 else "white", edgecolor="none", transform=t, clip_on=False))

            box_h = min(0.022, row_h * 0.70)
            y_box = ypos - (row_h * 0.50 + box_h * 0.45)

            leg_ax.add_patch(mpatches.Rectangle((xoff + 0.008, y_box), 0.022, box_h, facecolor="white", edgecolor="#333333", linewidth=0.6, transform=t, clip_on=False))
            leg_ax.text(xoff + 0.019, y_box + box_h / 2.0, str(i), ha="center", va="center", fontsize=effective_font - 0.2, fontweight="bold", transform=t, color="#202020")

            leg_ax.add_patch(mpatches.Rectangle((xoff + 0.036, y_box), 0.022, box_h, facecolor=palette.get(ct, "#888888"), edgecolor="#d2d2d2", linewidth=0.5, transform=t, clip_on=False))

            label = ct[:max_char] + ("..." if len(ct) > max_char else "")
            leg_ax.text(xoff + 0.066, y_box + box_h / 2.0, label, fontsize=effective_font, va="center", transform=t, color="#303030")

    fig = plt.figure(figsize=(28, 10), facecolor="white")
    plot_w = 0.29
    leg_w = 0.15
    gap = 0.02
    bottom = 0.06
    height = 0.88

    ax_adult = fig.add_axes([0.02, bottom, plot_w, height])
    ax_fetal = fig.add_axes([0.52, bottom, plot_w, height])

    ct_adult, pal_adult = _plot_umap_panel(
        adata_adult,
        ax_adult,
        title="The Tabula Sapiens Consortium et al., 2022",
        color_key="cell_type",
        tissue_key="tissue",
        panel_label="a",
        label_min_dist=0.068,
        label_max_shift=0.065,
    )
    _build_legend(
        fig,
        ax_adult,
        adata_adult,
        color_key="cell_type",
        tissue_key="tissue",
        cell_types=ct_adult,
        palette=pal_adult,
        legend_width=leg_w,
        gap=gap,
    )

    ct_fetal, pal_fetal = _plot_umap_panel(
        adata_fetal,
        ax_fetal,
        title="Yu et al., 2021",
        color_key="cell_type",
        tissue_key="tissue",
        panel_label="b",
        label_min_dist=0.068,
        label_max_shift=0.070,
    )
    _build_legend(
        fig,
        ax_fetal,
        adata_fetal,
        color_key="cell_type",
        tissue_key="tissue",
        cell_types=ct_fetal,
        palette=pal_fetal,
        legend_width=leg_w,
        gap=gap,
    )

    save_figure(fig, out_path, "two-panel UMAP", dpi=cfg_umap.get("fig_dpi", 300))


def _extract_markers_full(adata, stage_label: str, cfg_me: dict) -> PandasDataFrame:
    import scanpy as sc
    adata = adata.copy()
    if adata.obs["cell_type"].nunique() < 2:
        raise ValueError(f"[{stage_label}] Need ≥2 cell types for marker testing.")
    n_genes = adata.raw.n_vars if adata.raw is not None else adata.n_vars
    log.info("[%s] rank_genes_groups (n_genes=%d) ...", stage_label, n_genes)
    sc.tl.rank_genes_groups(
        adata, groupby="cell_type", method="wilcoxon",
        use_raw=True, pts=True, n_genes=n_genes, corr_method="benjamini-hochberg",
    )
    rgg = adata.uns["rank_genes_groups"]
    rows = []
    for ct in rgg["names"].dtype.names:
        for rank, (gene, lf, p_adj, pt, pt_r) in enumerate(zip(
            rgg["names"][ct], rgg["logfoldchanges"][ct],
            rgg["pvals_adj"][ct], rgg["pts"][ct], rgg["pts_rest"][ct],
        ), start=1):
            rows.append({"stage": stage_label, "cell_type": ct, "rank_raw": rank,
                         "gene": str(gene), "logfoldchange": float(lf),
                         "pval_adj": float(p_adj), "pct_cells": float(pt), "pct_rest": float(pt_r)})
    markers = pd.DataFrame(rows)
    markers = markers[markers["pval_adj"] < cfg_me["pval_adj_cutoff"]].copy()
    if cfg_me["require_positive_logfc"]:
        markers = markers[markers["logfoldchange"] > 0].copy()

    min_pct_cells = float(cfg_me.get("min_pct_cells", 0.0))
    max_pct_rest = float(cfg_me.get("max_pct_rest", 1.0))
    before_pct_gate = len(markers)
    if min_pct_cells > 0.0:
        markers = markers[markers["pct_cells"] >= min_pct_cells].copy()
    if max_pct_rest < 1.0:
        markers = markers[markers["pct_rest"] <= max_pct_rest].copy()
    removed_pct_gate = before_pct_gate - len(markers)
    if removed_pct_gate > 0:
        log.info(
            "[%s] Marker pct-gates removed %d rows (min_pct_cells>=%.3f, max_pct_rest<=%.3f)",
            stage_label,
            removed_pct_gate,
            min_pct_cells,
            max_pct_rest,
        )

    if markers.empty:
        raise ValueError(
            f"[{stage_label}] No marker rows remain after DE + pct filters. "
            "Relax marker_extraction.min_pct_cells / max_pct_rest thresholds."
        )

    markers = markers.sort_values(["cell_type", "pval_adj", "logfoldchange"], ascending=[True, True, False]).reset_index(drop=True)
    markers["rank_filtered"] = markers.groupby("cell_type").cumcount() + 1
    log.info("[%s] Extracted %d marker rows across %d cell types", stage_label, len(markers), markers["cell_type"].nunique())
    return markers


def _fetal_rescue(adata_fetal, cfg_me: dict, rescue_tsv_path: Path):
    import scanpy as sc
    generic_labels = {s.lower() for s in cfg_me["fetal_generic_kidney_labels"]}
    distal_markers: list = cfg_me["fetal_rescue_distal_markers"]
    loh_markers: list    = cfg_me["fetal_rescue_loh_markers"]
    label_distal: str    = cfg_me["fetal_rescue_label_distal"]
    label_loh: str       = cfg_me["fetal_rescue_label_loh"]
    min_margin: float = float(cfg_me.get("fetal_rescue_min_margin", 0.0))
    min_score: float = float(cfg_me.get("fetal_rescue_min_score", 0.0))

    ambiguous_strategy: str = str(cfg_me.get("fetal_rescue_ambiguous_strategy", "keep_original")).strip().lower()
    ambiguous_label: str = str(cfg_me.get("fetal_rescue_ambiguous_label", "ambiguous_tubule")).strip()

    adata_fetal.obs["cell_type"] = adata_fetal.obs["cell_type"].astype(str).str.strip()
    generic_mask = adata_fetal.obs["cell_type"].str.lower().isin(generic_labels)
    n_generic = int(generic_mask.sum())
    log.info("Fetal rescue: %d generic-label cells to classify", n_generic)
    if n_generic == 0:
        log.info("No generic fetal labels found — rescue skipped.")
        return adata_fetal

    available = set(adata_fetal.raw.var_names.astype(str) if adata_fetal.raw is not None else adata_fetal.var_names.astype(str))
    distal_present = [g for g in distal_markers if g in available]
    loh_present    = [g for g in loh_markers    if g in available]
    log.info("Rescue — distal genes present: %s", distal_present)
    log.info("Rescue — LoH genes present: %s", loh_present)
    if not distal_present or not loh_present:
        raise ValueError(f"Cannot rescue: zero marker genes present (distal={distal_present}, loh={loh_present}).")

    score_kw = {"use_raw": adata_fetal.raw is not None}
    sc.tl.score_genes(adata_fetal, gene_list=distal_present, score_name="_rescue_distal", **score_kw)
    sc.tl.score_genes(adata_fetal, gene_list=loh_present,    score_name="_rescue_loh",    **score_kw)

    rescue_df = adata_fetal.obs.loc[generic_mask, ["cell_type", "_rescue_distal", "_rescue_loh"]].copy()
    rescue_df.columns = ["original_cell_type", "rescue_distal_score", "rescue_loh_score"]
    base_assignment = np.where(
        rescue_df["rescue_distal_score"] >= rescue_df["rescue_loh_score"], label_distal, label_loh
    )
    rescue_df["rescued_cell_type"] = base_assignment
    rescue_df["score_margin"] = (rescue_df["rescue_distal_score"] - rescue_df["rescue_loh_score"]).abs()

    rescue_df["max_score"] = rescue_df[["rescue_distal_score", "rescue_loh_score"]].max(axis=1)
    
    # Flag as ambiguous if the margin is too tight OR if neither score is high enough
    ambiguous_mask = (rescue_df["score_margin"] < min_margin) | (rescue_df["max_score"] < min_score)
    # ----------------------------------

    if (min_margin > 0.0) or (min_score > -np.inf): # Update this if-statement
        if ambiguous_strategy == "keep_original":
            rescue_df.loc[ambiguous_mask, "rescued_cell_type"] = rescue_df.loc[ambiguous_mask, "original_cell_type"]
            rescue_df["rescue_decision"] = np.where(ambiguous_mask, "ambiguous_kept_original", "assigned")
        else:
            rescue_df.loc[ambiguous_mask, "rescued_cell_type"] = ambiguous_label
            rescue_df["rescue_decision"] = np.where(
                ambiguous_mask,
                f"ambiguous_assigned_label:{ambiguous_label}",
                "assigned",
            )
    else:
        rescue_df["rescue_decision"] = "assigned"

    log.info(
        "Fetal rescue confidence gate: %d/%d cells flagged ambiguous (min_margin=%.4f, min_score=%.4f, strategy=%s).",
        int(ambiguous_mask.sum()),
        len(rescue_df),
        min_margin,
        min_score,
        ambiguous_strategy,
    )

    adata_fetal.obs.loc[generic_mask, "cell_type"] = rescue_df["rescued_cell_type"].values
    adata_fetal.obs["cell_type"] = adata_fetal.obs["cell_type"].astype("category")
    for col in ["_rescue_distal", "_rescue_loh"]:
        if col in adata_fetal.obs.columns:
            adata_fetal.obs.drop(columns=[col], inplace=True)

    save_tsv(rescue_df.reset_index(), rescue_tsv_path, "fetal rescue assignments")
    log.info("Rescue counts: %s", rescue_df["rescued_cell_type"].value_counts().to_dict())
    return adata_fetal


def run_phase1(cfg: dict, paths: dict[str, Path]) -> dict[str, Path]:
    require_runtime_packages(
        {
            "matplotlib": plt,
            "numpy": np,
            "pandas": pd,
        },
        "Phase 1",
    )

    import anndata as ad
    import scanpy as sc
    sc.settings.verbosity = 1
    n_jobs = cfg["preprocessing"].get("n_jobs", 1)
    if not isinstance(n_jobs, int) or n_jobs == 0 or n_jobs < -1:
        raise ValueError("Config key 'preprocessing.n_jobs' must be -1 or a positive integer.")
    sc.settings.n_jobs = n_jobs
    ad.settings.allow_write_nullable_strings = True

    seed: int      = cfg.get("seed", 42)
    cfg_pre: dict  = cfg["preprocessing"]
    cfg_me: dict   = cfg["marker_extraction"]
    cfg_umap: dict = cfg["umap"]
    fetal_top_n_main: int = int(cfg["hpa_validation"]["top_n_validated"])
    marker_min_cells: int = resolve_marker_calling_min_cells(cfg)
    pval_eps: float = float(cfg_me["pval_epsilon"])
    noise_prefixes  = tuple(cfg_pre["noise_gene_prefixes"])
    noise_exact     = set(cfg_pre["noise_gene_exact"])

    output_dir  = paths["output_dir"]
    aux_dir     = paths["aux_dir"]
    figures_dir = paths["figures_dir"]
    ensure_dirs(output_dir, aux_dir, figures_dir)
    log.info("Phase 1 runtime settings: scanpy.n_jobs=%d", sc.settings.n_jobs)

    # Load
    require_file(paths["adult_h5ad"], "adult h5ad")
    require_file(paths["fetal_h5ad"], "fetal h5ad")
    adata_adult = _ensure_gene_symbols(sc.read_h5ad(str(paths["adult_h5ad"])), "adult")
    adata_fetal = _ensure_gene_symbols(sc.read_h5ad(str(paths["fetal_h5ad"])), "fetal")
    log.info("Adult raw: %d cells × %d genes", adata_adult.n_obs, adata_adult.n_vars)
    log.info("Fetal raw: %d cells × %d genes", adata_fetal.n_obs, adata_fetal.n_vars)
    for adata, lbl in [(adata_adult, "adult"), (adata_fetal, "fetal")]:
        if "cell_type" not in adata.obs.columns:
            raise KeyError(f"Reference [{lbl}] missing required obs column 'cell_type'.")

    # Filter labels
    adata_adult = _remove_labels(adata_adult, {s.lower() for s in cfg_pre["adult_labels_to_remove"]}, "adult")
    adata_fetal = _remove_labels(adata_fetal, {s.lower() for s in cfg_pre["fetal_labels_to_remove"]}, "fetal")
    adata_adult = _collapse_adult_endothelial_cell_types(adata_adult)

    # Preprocess
    adata_adult = _preprocess_reference(adata_adult, "adult", cfg_pre, cfg_me, seed)
    adata_fetal = _preprocess_reference(adata_fetal, "fetal", cfg_pre, cfg_me, seed)
    adata_adult = _filter_noise_genes(adata_adult, noise_prefixes, noise_exact, "adult")
    adata_fetal = _filter_noise_genes(adata_fetal, noise_prefixes, noise_exact, "fetal")

    if cfg_pre["apply_min_celltype_filter"]:
        min_c = cfg_pre["min_cells_per_type"]
        adata_adult = _enforce_min_cells_per_type(adata_adult, min_c, "adult")
        adata_fetal = _enforce_min_cells_per_type(adata_fetal, min_c, "fetal")

    enforce_marker_calling_min_cells(adata_adult, marker_min_cells, "adult")

    # Checkpoint h5ad
    adult_h5ad_out = aux_dir / "adult_preprocessed.h5ad"
    fetal_h5ad_out = aux_dir / "fetal_preprocessed.h5ad"
    adata_adult.write_h5ad(str(adult_h5ad_out))
    log.info("Checkpoint saved: %s", adult_h5ad_out)
    adata_fetal.write_h5ad(str(fetal_h5ad_out))
    log.info("Checkpoint saved: %s", fetal_h5ad_out)

    # UMAP figure immediately after checkpoint
    _plot_umap_two_panel(adata_adult, adata_fetal, figures_dir / "umap_adult_fetal_two_panel.png", cfg_umap)

    # Marker extraction — adult
    markers_adult = _extract_markers_full(adata_adult, "adult", cfg_me)
    markers_adult = exclude_adult_stress_signature_markers(markers_adult, cfg_me)

    # Fetal rescue → recompute fetal markers
    rescue_tsv = aux_dir / "fetal_generic_kidney_rescue_assignments.tsv"
    adata_fetal = _fetal_rescue(adata_fetal, cfg_me, rescue_tsv)
    adata_fetal = _drop_fetal_post_rescue_generic_labels(
        adata_fetal,
        {"kidney cell", "kidney epithelia cell", "kidney epithelial cell"},
    )
    enforce_marker_calling_min_cells(adata_fetal, marker_min_cells, "fetal")
    markers_fetal = _extract_markers_full(adata_fetal, "fetal", cfg_me)

    # Main output top table (fetal)
    min_pct_cells = float(cfg_me.get("min_pct_cells", 0.0))
    max_pct_rest = float(cfg_me.get("max_pct_rest", 1.0))
    markers_fetal_top = build_top_markers(
        markers_fetal,
        fetal_top_n_main,
        pval_eps,
        min_pct_cells=min_pct_cells,
        max_pct_rest=max_pct_rest,
    )

    fetal_top_tsv = output_dir / f"fetal_kidney_markers_top{fetal_top_n_main}.tsv"
    save_tsv(markers_fetal_top, fetal_top_tsv, f"fetal top{fetal_top_n_main} markers")

    # All-markers → aux outputs
    adult_all_tsv = aux_dir / "adult_kidney_markers_all.tsv"
    fetal_all_tsv = aux_dir / "fetal_kidney_markers_all.tsv"
    save_tsv(markers_adult, adult_all_tsv, "adult all markers")
    save_tsv(markers_fetal, fetal_all_tsv, "fetal all markers")

    log.info("Phase 1 complete.")
    return {
        "adult_markers_all_tsv": adult_all_tsv,
        "fetal_markers_all_tsv": fetal_all_tsv,
        "fetal_markers_top_tsv": fetal_top_tsv,
        "fetal_rescue_tsv":      rescue_tsv,
        "adult_h5ad":            adult_h5ad_out,
        "fetal_h5ad":            fetal_h5ad_out,
        "umap_figure":           figures_dir / "umap_adult_fetal_two_panel.png",
    }


# =============================================================================
# SECTION 4 — Phase 2: Adult HPA Validation + Ranking
# =============================================================================

_MARKER_COLS = ["stage", "cell_type", "gene", "logfoldchange", "pct_cells", "pct_rest", "pval_adj"]
_HPA_COLS    = ["Gene", "Gene name", "Cell type", "nCPM"]


def _load_hpa(hpa_path: Path) -> dict[str, PandasDataFrame]:
    """Load HPA TSV and return a gene-name → subset lookup dict."""
    require_file(hpa_path, "HPA reference TSV")
    hpa_raw = pd.read_csv(hpa_path, sep="\t", low_memory=False)
    require_columns(hpa_raw, _HPA_COLS, "HPA TSV")
    hpa = hpa_raw[_HPA_COLS].copy()
    hpa["Gene name"] = hpa["Gene name"].astype(str).str.strip()
    hpa["Cell type"]  = hpa["Cell type"].astype(str).str.strip()
    hpa["nCPM"]       = pd.to_numeric(hpa["nCPM"], errors="coerce")
    hpa = hpa.dropna(subset=["Gene name", "Cell type", "nCPM"])
    hpa["gene_name_lower"]  = hpa["Gene name"].str.lower()
    hpa["cell_type_lower"]  = hpa["Cell type"].str.lower()
    lookup = {g: sub[["Cell type", "cell_type_lower", "nCPM"]]
              for g, sub in hpa.groupby("gene_name_lower", sort=False)}
    log.info("HPA loaded: %d rows, %d unique genes", len(hpa), len(lookup))
    return lookup


def _validate_row(row: PandasSeries, gene_lookup: dict, hpa_mapping: dict,
                  ncpm_threshold: float, manual_override_cts: set, manual_markers: set) -> dict:
    my_ct    = str(row["cell_type"]).strip()
    my_ct_lc = my_ct.lower()
    gene     = str(row["gene"]).strip()

    def _lookup_hpa_mapping(label: str) -> str | None:
        direct = hpa_mapping.get(label)
        if direct is not None:
            return direct
        label_lc = str(label).strip().lower()
        for k, v in hpa_mapping.items():
            if str(k).strip().lower() == label_lc:
                return v
        return None

    result: dict[str, Any] = {
        "my_cell_type": my_ct, "gene": gene,
        "logfoldchange": row.get("logfoldchange", np.nan),
        "pct_cells":     row.get("pct_cells", np.nan),
        "pct_rest":      row.get("pct_rest", np.nan),
        "pval_adj":      row.get("pval_adj", np.nan),
        "hpa_mapped_cell_type": np.nan, "hpa_nCPM": np.nan,
        "hpa_specificity": np.nan, "is_validated": False, "validation_note": "Unmapped",
    }

    if my_ct_lc in manual_override_cts:
        mapped = _lookup_hpa_mapping(my_ct)
        if mapped is None:
            raise KeyError(
                f"Manual override cell type '{my_ct}' is missing in "
                "hpa_validation.hpa_cell_type_mapping. Add an explicit mapping in config.yaml."
            )
            
        # If the gene is on the manual list, validate it immediately and return.
        if gene.upper() in manual_markers:
            result["hpa_mapped_cell_type"] = mapped
            result.update(is_validated=True, hpa_nCPM=float(ncpm_threshold),
                          hpa_specificity=1.0, validation_note="Validated (manual override)")
            return result

    mapped = _lookup_hpa_mapping(my_ct)
    if mapped is None:
        result["validation_note"] = "Unmapped cell type"
        return result

    result["hpa_mapped_cell_type"] = mapped
    gene_subset = gene_lookup.get(gene.lower())
    if gene_subset is None or gene_subset.empty:
        result["validation_note"] = "Gene not found in HPA"
        return result

    pattern = re.escape(str(mapped).strip().lower())
    matched = gene_subset[gene_subset["cell_type_lower"].str.contains(pattern, case=False, na=False, regex=True)]
    if matched.empty:
        result["validation_note"] = "No matching HPA cell type for mapped label"
        return result

    ncpm_val = float(matched["nCPM"].max())
    sum_maxes = float(gene_subset.groupby("cell_type_lower")["nCPM"].max().sum())
    specificity = (ncpm_val / sum_maxes) if sum_maxes > 0 else 0.0
    result.update(
        hpa_nCPM=ncpm_val, hpa_specificity=float(specificity),
        is_validated=bool(ncpm_val >= ncpm_threshold),
        validation_note="Validated" if ncpm_val >= ncpm_threshold else f"Below threshold (< {ncpm_threshold})",
    )
    return result


def run_phase2(cfg: dict, paths: dict[str, Path], phase1_artifacts: dict[str, Path]) -> dict[str, Path]:
    require_runtime_packages({"numpy": np, "pandas": pd}, "Phase 2")

    gate_phase({"adult_markers_all_tsv": phase1_artifacts["adult_markers_all_tsv"]}, "Phase 2")

    cfg_hpa: dict    = cfg["hpa_validation"]
    aux_dir: Path    = paths["aux_dir"]
    diagnostics_dir: Path = paths["diagnostics_dir"] / "adult_hpa_validation"
    ensure_dirs(aux_dir, diagnostics_dir)

    require_file(paths["hpa_tsv"], "HPA reference TSV")
    gene_lookup = _load_hpa(paths["hpa_tsv"])

    ncpm_threshold: float = cfg_hpa["validation_ncpm_threshold"]
    top_n: int            = cfg_hpa["top_n_validated"]
    pval_eps: float       = float(cfg_hpa.get("pval_epsilon", cfg["marker_extraction"].get("pval_epsilon", 1e-300)))
    manual_override_cts   = {s.lower() for s in cfg_hpa["manual_override_cell_types"]}
    manual_markers        = {s.upper() for s in cfg_hpa["manual_validated_markers"]}
    hpa_mapping: dict     = cfg_hpa["hpa_cell_type_mapping"]

    all_markers_raw = pd.read_csv(phase1_artifacts["adult_markers_all_tsv"], sep="\t")
    require_columns(all_markers_raw, _MARKER_COLS, "adult_kidney_markers_all.tsv")
    all_markers_df = all_markers_raw[_MARKER_COLS].copy()
    all_markers_df["cell_type"] = all_markers_df["cell_type"].astype(str).str.strip()
    all_markers_df["gene"]      = all_markers_df["gene"].astype(str).str.strip()

    log.info("Adult HPA validation: validating %d rows ...", len(all_markers_df))
    records = [
        _validate_row(row, gene_lookup, hpa_mapping, ncpm_threshold, manual_override_cts, manual_markers)
        for _, row in all_markers_df.iterrows()
    ]
    validated_df = pd.DataFrame(records)
    save_tsv(validated_df, diagnostics_dir / "adult_hpa_validation_audit.tsv", "adult HPA validation audit")

    # Post-validation filter
    manual_mask = validated_df["validation_note"].str.contains("manual override", case=False, na=False)
    rank_src = validated_df[
        validated_df["is_validated"]
        & (((validated_df["pval_adj"] < cfg_hpa["pval_threshold"]) &
             (validated_df["logfoldchange"] > cfg_hpa["logfc_threshold"])) | manual_mask)
    ].copy()
    if rank_src.empty:
        raise ValueError("Adult HPA validation: no markers passed p-value/logFC + HPA validation filters.")

    # Fill manual-override rows with sentinel HPA metrics
    rank_manual_mask = rank_src["validation_note"].str.contains("manual override", case=False, na=False)
    for col, fill in [("hpa_nCPM", ncpm_threshold), ("hpa_specificity", 1.0)]:
        rank_src.loc[rank_manual_mask & rank_src[col].isna(), col] = fill
    rank_src = rank_src.dropna(subset=["pval_adj", "logfoldchange", "hpa_nCPM"]).copy()

    # Joint ranking
    rank_src["pval_adj_safe"]    = rank_src["pval_adj"].clip(lower=pval_eps)
    rank_src["neg_log10_padj"]   = -np.log10(rank_src["pval_adj_safe"])
    rank_src["p_rank"]   = rank_src.groupby("my_cell_type")["pval_adj_safe"].rank(method="min", ascending=True)
    rank_src["lfc_rank"] = rank_src.groupby("my_cell_type")["logfoldchange"].rank(method="min", ascending=False)
    rank_src["combined_rank_score"] = rank_src["p_rank"] + rank_src["lfc_rank"]

    def _scale(series: PandasSeries) -> PandasSeries:
        mn, mx = series.min(), series.max()
        return (series - mn) / (mx - mn) if mx > mn else pd.Series(1.0, index=series.index)

    rank_src["extraction_strength_scaled"] = _scale(1.0 / rank_src["combined_rank_score"])
    rank_src["hpa_strength_scaled"]        = _scale(np.log1p(rank_src["hpa_nCPM"].clip(lower=0.0)))

    w_ext  = float(cfg_hpa["weight_extraction_strength"])
    w_hpa  = float(cfg_hpa["weight_hpa_strength"])
    w_spec = float(cfg_hpa["weight_hpa_specificity"])
    rank_src["final_joint_score"] = (
        w_ext  * rank_src["extraction_strength_scaled"]
        + w_hpa  * rank_src["hpa_strength_scaled"]
        + w_spec * rank_src["hpa_specificity"]
    )

    ranked = rank_src.sort_values(
        ["my_cell_type", "final_joint_score", "combined_rank_score",
         "hpa_specificity", "neg_log10_padj", "logfoldchange", "hpa_nCPM"],
        ascending=[True, False, True, False, False, False, False],
    ).reset_index(drop=True)
    ranked["final_rank_within_cell_type"] = ranked.groupby("my_cell_type").cumcount() + 1

    save_tsv(ranked, diagnostics_dir / "adult_hpa_ranked_candidates.tsv", "adult HPA ranked candidates")

    top_n_df = ranked[ranked["final_rank_within_cell_type"] <= top_n].copy()
    top_n_aux_path = aux_dir / f"adult_kidney_markers_top{top_n}_validated.tsv"
    top_n_main_path = paths["output_dir"] / f"adult_kidney_markers_top{top_n}_validated.tsv"
    save_tsv(top_n_df, top_n_aux_path, f"adult top{top_n} validated markers (aux)")
    save_tsv(top_n_df, top_n_main_path, f"adult top{top_n} validated markers (main output)")

    log.info("Phase 2 complete — adult HPA validation ranking finished.")
    return {
        "adult_top_validated_tsv": top_n_aux_path,
        "adult_top_validated_main_tsv": top_n_main_path,
        "adult_hpa_diagnostics_dir": diagnostics_dir,
    }


# =============================================================================
# SECTION 5 — Phase 3: Overlap + YAML Export
# =============================================================================

def _load_standardise_markers(
    path: Path, stage: str,
    cell_col_candidates: list, gene_col_candidates: list, rank_col_candidates: list,
    harmonization_lookup: dict,
) -> PandasDataFrame:
    df = pd.read_csv(path, sep="\t")

    def pick(candidates, label):
        for c in candidates:
            if c in df.columns:
                return c
        raise KeyError(f"Missing {label} column. Tried {candidates}; got {list(df.columns)}")

    cell_col = pick(cell_col_candidates, f"{stage} cell_type")
    gene_col = pick(gene_col_candidates, f"{stage} gene")
    rank_col = next((c for c in rank_col_candidates if c in df.columns), None)

    out = df.copy()
    out["stage"]         = stage
    out["cell_type"]     = out[cell_col].astype(str).str.strip()
    out["gene"]          = out[gene_col].astype(str).str.strip()
    out["cell_type_norm"] = out["cell_type"].map(lambda x: harmonize_cell_type(x, harmonization_lookup))
    out["rank"] = pd.to_numeric(out[rank_col], errors="coerce") if rank_col else out.groupby("cell_type").cumcount() + 1
    out = out.dropna(subset=["cell_type", "gene"])
    out = out[out["gene"] != ""]
    return out


def _ranked_genes_for_norm(df: PandasDataFrame, norm_name: str) -> list[str]:
    sub = df[df["cell_type_norm"] == norm_name].copy()
    if sub.empty:
        return []

    sub["gene_upper"] = sub["gene"].astype(str).str.upper()

    sort_specs: list[tuple[str, bool]] = []
    if "rank" in sub.columns:
        sort_specs.append(("rank", True))
    if "final_joint_score" in sub.columns:
        sort_specs.append(("final_joint_score", False))
    if "combined_rank_score" in sub.columns:
        sort_specs.append(("combined_rank_score", True))
    if "pval_adj" in sub.columns:
        sort_specs.append(("pval_adj", True))
    if "logfoldchange" in sub.columns:
        sort_specs.append(("logfoldchange", False))
    sort_specs.extend([("cell_type", True), ("gene_upper", True)])

    sort_cols = [c for c, _ in sort_specs if c in sub.columns]
    sort_ascending = [a for c, a in sort_specs if c in sub.columns]
    sub = sub.sort_values(sort_cols, ascending=sort_ascending, na_position="last")
    sub = sub.drop_duplicates("gene_upper", keep="first")
    return sub["gene_upper"].tolist()


def _grouped_ranked_genes(df: PandasDataFrame, top_n: int | None = None) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {}
    for norm in sorted(df["cell_type_norm"].astype(str).unique()):
        genes = _ranked_genes_for_norm(df, norm)
        grouped[norm] = genes if top_n is None else genes[:top_n]
    return grouped


def _ordered_intersection(primary_ranked: list[str], secondary_ranked: list[str]) -> list[str]:
    primary_set = set(primary_ranked)
    secondary_set = set(secondary_ranked)
    shared = [g for g in primary_ranked if g in secondary_set]
    seen = set(shared)
    shared.extend([g for g in secondary_ranked if g in primary_set and g not in seen])
    return shared


def _summarize_cell_types_from_h5ad(h5ad_path: Path, stage_label: str) -> PandasDataFrame:
    import anndata as ad

    adata = ad.read_h5ad(str(h5ad_path))
    if "cell_type" not in adata.obs.columns:
        raise KeyError(f"Missing 'cell_type' in {h5ad_path}")

    counts = adata.obs["cell_type"].astype(str).value_counts().rename("n_cells")
    summary = counts.to_frame()
    summary["fraction"] = summary["n_cells"] / summary["n_cells"].sum()
    summary["stage"] = stage_label
    return summary.reset_index().rename(columns={"index": "cell_type"})


_DEFAULT_MANUAL_MAPPING_SEED: list[dict[str, Any]] = [
    {
        "fetal_cell_type": "endothelial cell",
        "adult_cell_type_candidates": [
            "peritubular capillary endothelial cell",
            "glomerular capillary endothelial cell",
            "endothelial cell of lymphatic vessel",
        ],
        "mapping_type": "manual",
        "confidence": "medium",
        "is_broad_generic": False,
    },
    {
        "fetal_cell_type": "glomerular mesangial cell",
        "adult_cell_type_candidates": ["mural cell"],
        "mapping_type": "manual",
        "confidence": "low",
        "is_broad_generic": False,
    },
    {
        "fetal_cell_type": "kidney interstitial cell",
        "adult_cell_type_candidates": [
            "kidney interstitial fibroblast",
            "renal medullary fibroblast",
            "mural cell",
        ],
        "mapping_type": "manual",
        "confidence": "medium",
        "is_broad_generic": True,
    },
]


def _build_mapping_aware_comparison(
    adult_summary: PandasDataFrame,
    fetal_summary: PandasDataFrame,
    adult_marker_types: set[str],
    fetal_marker_types: set[str],
    manual_mapping_seed: list[dict[str, Any]] | None = None,
) -> tuple[PandasDataFrame, PandasDataFrame, PandasDataFrame, PandasDataFrame]:
    adult_types = set(adult_summary.index)
    fetal_types = set(fetal_summary.index)
    exact_types = sorted(adult_types & fetal_types)
    mapping_seed = manual_mapping_seed if manual_mapping_seed is not None else _DEFAULT_MANUAL_MAPPING_SEED

    mapping_rows = []
    for ct in exact_types:
        mapping_rows.append(
            {
                "fetal_cell_type": ct,
                "adult_cell_type": ct,
                "mapping_type": "exact",
                "confidence": "high",
                "is_broad_generic": False,
            }
        )

    for row in mapping_seed:
        for adult_ct in row["adult_cell_type_candidates"]:
            if row["fetal_cell_type"] in fetal_types and adult_ct in adult_types:
                mapping_rows.append(
                    {
                        "fetal_cell_type": row["fetal_cell_type"],
                        "adult_cell_type": adult_ct,
                        "mapping_type": row["mapping_type"],
                        "confidence": row["confidence"],
                        "is_broad_generic": bool(row["is_broad_generic"]),
                    }
                )

    mapping_df = pd.DataFrame(mapping_rows).drop_duplicates().reset_index(drop=True)
    if mapping_df.empty:
        mapping_df = pd.DataFrame(
            columns=[
                "fetal_cell_type",
                "adult_cell_type",
                "mapping_type",
                "confidence",
                "is_broad_generic",
                "pair_label",
            ]
        )
    mapping_df["pair_label"] = mapping_df["fetal_cell_type"].astype(str) + " <-> " + mapping_df["adult_cell_type"].astype(str)

    shared_adult_mapped = sorted(mapping_df["adult_cell_type"].dropna().astype(str).unique().tolist())
    shared_fetal_mapped = sorted(mapping_df["fetal_cell_type"].dropna().astype(str).unique().tolist())
    adult_only = sorted(adult_types - set(shared_adult_mapped))
    fetal_only = sorted(fetal_types - set(shared_fetal_mapped))

    sharing_summary_df = pd.DataFrame(
        [
            {
                "category": "shared_mapped",
                "n_cell_types": len(shared_fetal_mapped),
                "definition": "fetal cell types with >=1 mapped adult counterpart",
            },
            {
                "category": "adult_only",
                "n_cell_types": len(adult_only),
                "definition": "adult cell types without mapped fetal counterpart",
            },
            {
                "category": "fetal_only",
                "n_cell_types": len(fetal_only),
                "definition": "fetal cell types without mapped adult counterpart",
            },
            {
                "category": "shared_exact_labels",
                "n_cell_types": len(exact_types),
                "definition": "exact label intersection only (reference)",
            },
        ]
    )

    sharing_details_rows = []
    for ct in exact_types:
        sharing_details_rows.append({"category": "shared_exact_labels", "cell_type": ct})
    for ct in shared_fetal_mapped:
        sharing_details_rows.append({"category": "shared_mapped_fetal", "cell_type": ct})
    for ct in shared_adult_mapped:
        sharing_details_rows.append({"category": "shared_mapped_adult", "cell_type": ct})
    for ct in adult_only:
        sharing_details_rows.append({"category": "adult_only", "cell_type": ct})
    for ct in fetal_only:
        sharing_details_rows.append({"category": "fetal_only", "cell_type": ct})
    sharing_details_df = pd.DataFrame(sharing_details_rows)

    comparison_df = (
        mapping_df.merge(
            fetal_summary[["n_cells", "fraction"]].rename(columns={"n_cells": "n_cells_fetal", "fraction": "fraction_fetal"}),
            left_on="fetal_cell_type",
            right_index=True,
            how="left",
        )
        .merge(
            adult_summary[["n_cells", "fraction"]].rename(columns={"n_cells": "n_cells_adult", "fraction": "fraction_adult"}),
            left_on="adult_cell_type",
            right_index=True,
            how="left",
        )
    )
    comparison_df["fetal_has_markers"] = comparison_df["fetal_cell_type"].isin(fetal_marker_types)
    comparison_df["adult_has_markers"] = comparison_df["adult_cell_type"].isin(adult_marker_types)
    comparison_df["is_exact_label_match"] = comparison_df["fetal_cell_type"] == comparison_df["adult_cell_type"]

    return comparison_df, mapping_df, sharing_summary_df, sharing_details_df


def _build_cell_type_label_mapping_table(adult_df: PandasDataFrame, fetal_df: PandasDataFrame) -> PandasDataFrame:
    """Create an explicit adult/fetal raw-label to harmonized-vocabulary mapping table."""

    def _summarize(df: PandasDataFrame, stage: str) -> PandasDataFrame:
        tmp = df.copy()
        tmp["gene_upper"] = tmp["gene"].astype(str).str.upper()
        summary = (
            tmp.groupby(["cell_type", "cell_type_norm"], as_index=False)
            .agg(
                n_marker_rows=("gene_upper", "size"),
                n_unique_markers=("gene_upper", "nunique"),
                min_rank=("rank", "min"),
                median_rank=("rank", "median"),
            )
            .rename(
                columns={
                    "cell_type": "raw_cell_type",
                    "cell_type_norm": "harmonized_cell_type",
                }
            )
        )
        summary.insert(0, "stage", stage)
        return summary

    combined = pd.concat([
        _summarize(adult_df, "adult"),
        _summarize(fetal_df, "fetal"),
    ], ignore_index=True)
    combined = combined.sort_values(
        ["stage", "harmonized_cell_type", "n_marker_rows", "raw_cell_type"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)
    return combined


def _build_overlap_acceptance_report(
    overlap_df: PandasDataFrame,
    sharing_summary_df: PandasDataFrame,
    cfg_overlap: dict[str, Any],
) -> tuple[PandasDataFrame, list[str]]:
    """
    Build acceptance report and list of failing configured criteria.

    Supported criteria keys in overlap.acceptance_criteria:
      - min_best_mapped_jaccard
      - min_best_mapped_shared_genes
      - expected_shared_mapped_min / expected_shared_mapped_max
      - expected_adult_only_min / expected_adult_only_max
      - expected_fetal_only_min / expected_fetal_only_max
    """
    criteria = cfg_overlap.get("acceptance_criteria", {})
    if not isinstance(criteria, dict):
        criteria = {}

    summary_counts = {
        str(k): int(v)
        for k, v in sharing_summary_df.set_index("category")["n_cell_types"].to_dict().items()
    }
    shared_mapped = int(summary_counts.get("shared_mapped", 0))
    adult_only = int(summary_counts.get("adult_only", 0))
    fetal_only = int(summary_counts.get("fetal_only", 0))

    mapped_pairs = overlap_df[overlap_df["is_mapped_pair"]].copy()
    best_mapped_jaccard = float(mapped_pairs["jaccard"].max()) if not mapped_pairs.empty else 0.0
    best_mapped_shared = int(mapped_pairs["n_shared"].max()) if not mapped_pairs.empty else 0

    rows: list[dict[str, Any]] = []

    def _add_check(key: str, observed: float | int, operator: str, expected: Any) -> None:
        configured = expected is not None
        if not configured:
            passed = True
        elif operator == ">=":
            passed = float(observed) >= float(expected)
        elif operator == "<=":
            passed = float(observed) <= float(expected)
        else:
            raise ValueError(f"Unsupported operator in acceptance check: {operator}")

        rows.append(
            {
                "criterion": key,
                "observed": observed,
                "operator": operator,
                "expected": expected,
                "configured": configured,
                "passed": bool(passed),
            }
        )

    _add_check("min_best_mapped_jaccard", best_mapped_jaccard, ">=", criteria.get("min_best_mapped_jaccard"))
    _add_check("min_best_mapped_shared_genes", best_mapped_shared, ">=", criteria.get("min_best_mapped_shared_genes"))

    _add_check("expected_shared_mapped_min", shared_mapped, ">=", criteria.get("expected_shared_mapped_min"))
    _add_check("expected_shared_mapped_max", shared_mapped, "<=", criteria.get("expected_shared_mapped_max"))

    _add_check("expected_adult_only_min", adult_only, ">=", criteria.get("expected_adult_only_min"))
    _add_check("expected_adult_only_max", adult_only, "<=", criteria.get("expected_adult_only_max"))

    _add_check("expected_fetal_only_min", fetal_only, ">=", criteria.get("expected_fetal_only_min"))
    _add_check("expected_fetal_only_max", fetal_only, "<=", criteria.get("expected_fetal_only_max"))

    rows.extend(
        [
            {
                "criterion": "observed_shared_mapped",
                "observed": shared_mapped,
                "operator": "info",
                "expected": None,
                "configured": False,
                "passed": True,
            },
            {
                "criterion": "observed_adult_only",
                "observed": adult_only,
                "operator": "info",
                "expected": None,
                "configured": False,
                "passed": True,
            },
            {
                "criterion": "observed_fetal_only",
                "observed": fetal_only,
                "operator": "info",
                "expected": None,
                "configured": False,
                "passed": True,
            },
            {
                "criterion": "observed_best_mapped_jaccard",
                "observed": best_mapped_jaccard,
                "operator": "info",
                "expected": None,
                "configured": False,
                "passed": True,
            },
            {
                "criterion": "observed_best_mapped_shared_genes",
                "observed": best_mapped_shared,
                "operator": "info",
                "expected": None,
                "configured": False,
                "passed": True,
            },
        ]
    )

    report_df = pd.DataFrame(rows)
    failing = report_df[(report_df["configured"]) & (~report_df["passed"])].copy()
    failed_messages = [
        (
            f"{row['criterion']} failed: observed={row['observed']} "
            f"{row['operator']} expected={row['expected']}"
        )
        for _, row in failing.iterrows()
    ]
    return report_df, failed_messages


def _grouped_gene_sets_by_cell_type(df: PandasDataFrame) -> dict[str, set[str]]:
    grouped = {}
    for ct in sorted(df["cell_type"].astype(str).unique()):
        sub = df[df["cell_type"].astype(str) == ct].copy()
        sub = sub.sort_values(["rank", "gene"], na_position="last")
        genes = set(sub["gene"].astype(str).str.upper().drop_duplicates().tolist())
        grouped[ct] = genes
    return grouped


def _build_overlap_table_all_pairs(adult_df: PandasDataFrame, fetal_df: PandasDataFrame, mapping_df: PandasDataFrame) -> PandasDataFrame:
    adult_sets = _grouped_gene_sets_by_cell_type(adult_df)
    fetal_sets = _grouped_gene_sets_by_cell_type(fetal_df)

    overlap_rows = []
    for adult_ct, adult_genes in adult_sets.items():
        for fetal_ct, fetal_genes in fetal_sets.items():
            shared = adult_genes & fetal_genes
            union = adult_genes | fetal_genes
            jaccard = len(shared) / len(union) if len(union) > 0 else 0.0
            overlap_rows.append(
                {
                    "adult_cell_type": adult_ct,
                    "fetal_cell_type": fetal_ct,
                    "n_adult_markers": len(adult_genes),
                    "n_fetal_markers": len(fetal_genes),
                    "n_shared": len(shared),
                    "jaccard": round(jaccard, 6),
                    "shared_genes": ",".join(sorted(shared)),
                }
            )

    overlap_df = pd.DataFrame(overlap_rows)

    if not mapping_df.empty:
        map_cols = ["fetal_cell_type", "adult_cell_type", "mapping_type", "confidence", "is_broad_generic"]
        overlap_df = overlap_df.merge(mapping_df[map_cols].drop_duplicates(), on=["fetal_cell_type", "adult_cell_type"], how="left")
    else:
        overlap_df["mapping_type"] = np.nan
        overlap_df["confidence"] = np.nan
        overlap_df["is_broad_generic"] = np.nan

    overlap_df["is_mapped_pair"] = overlap_df["mapping_type"].notna()
    overlap_df["mapping_type"] = overlap_df["mapping_type"].fillna("unmapped")
    overlap_df["confidence"] = overlap_df["confidence"].fillna("none")
    overlap_df["is_broad_generic"] = overlap_df["is_broad_generic"].fillna(False)

    overlap_df = overlap_df.sort_values(["is_mapped_pair", "jaccard", "n_shared"], ascending=[False, False, False]).reset_index(drop=True)
    return overlap_df


def _plot_overlap_figure_mapping_aware(
    overlap_df: PandasDataFrame,
    mapping_df: PandasDataFrame,
    adult_summary: PandasDataFrame,
    fetal_summary: PandasDataFrame,
    out_path: Path,
    dpi: int = 300,
) -> None:
    adult_abundance = adult_summary["n_cells"].sort_values(ascending=False).index.tolist()
    fetal_abundance = fetal_summary["n_cells"].sort_values(ascending=False).index.tolist()

    mapped_adult_non_generic = []
    mapped_fetal_non_generic = []
    if not mapping_df.empty:
        mapped_adult_non_generic = mapping_df.loc[~mapping_df["is_broad_generic"], "adult_cell_type"].drop_duplicates().tolist()
        mapped_fetal_non_generic = mapping_df.loc[~mapping_df["is_broad_generic"], "fetal_cell_type"].drop_duplicates().tolist()

    nonzero_adult = set(overlap_df.loc[overlap_df["n_shared"] > 0, "adult_cell_type"])
    nonzero_fetal = set(overlap_df.loc[overlap_df["n_shared"] > 0, "fetal_cell_type"])

    adult_col_order = [ct for ct in mapped_adult_non_generic if ct in nonzero_adult]
    adult_col_order += [ct for ct in adult_abundance if ct in nonzero_adult and ct not in adult_col_order]
    adult_col_order += [ct for ct in sorted(nonzero_adult) if ct not in adult_col_order]

    fetal_row_order = [ct for ct in mapped_fetal_non_generic if ct in nonzero_fetal]
    fetal_row_order += [ct for ct in fetal_abundance if ct in nonzero_fetal and ct not in fetal_row_order]
    fetal_row_order += [ct for ct in sorted(nonzero_fetal) if ct not in fetal_row_order]

    if not adult_col_order:
        adult_col_order = [ct for ct in mapped_adult_non_generic if ct in set(overlap_df["adult_cell_type"])]
        adult_col_order += [ct for ct in adult_abundance if ct in set(overlap_df["adult_cell_type"]) and ct not in adult_col_order]
        adult_col_order += [ct for ct in sorted(set(overlap_df["adult_cell_type"])) if ct not in adult_col_order]

    if not fetal_row_order:
        fetal_row_order = [ct for ct in mapped_fetal_non_generic if ct in set(overlap_df["fetal_cell_type"])]
        fetal_row_order += [ct for ct in fetal_abundance if ct in set(overlap_df["fetal_cell_type"]) and ct not in fetal_row_order]
        fetal_row_order += [ct for ct in sorted(set(overlap_df["fetal_cell_type"])) if ct not in fetal_row_order]

    pivot_j = overlap_df.pivot(index="fetal_cell_type", columns="adult_cell_type", values="jaccard")
    pivot_n = overlap_df.pivot(index="fetal_cell_type", columns="adult_cell_type", values="n_shared")
    pivot_j = pivot_j.reindex(index=fetal_row_order, columns=adult_col_order)
    pivot_n = pivot_n.reindex(index=fetal_row_order, columns=adult_col_order)

    heat_values = pivot_j.fillna(0.0)
    shared_counts = pivot_n.fillna(0).astype(int)
    mask_zero = shared_counts == 0
    annot = shared_counts.astype(str).mask(mask_zero, "")

    if heat_values.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No overlap values available to plot.", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        save_figure(fig, out_path, "marker overlap figure", dpi=dpi)
        return

    max_jaccard = float(np.nanmax(heat_values.values)) if np.isfinite(heat_values.values).any() else 0.0

    wrap_width = 26
    x_labels = ["\n".join(textwrap.wrap(x, width=wrap_width)) for x in heat_values.columns]
    y_labels = ["\n".join(textwrap.wrap(y, width=wrap_width)) for y in heat_values.index]

    fig_w = max(10, len(adult_col_order) * 0.8 + 3)
    fig_h = max(6, len(fetal_row_order) * 0.8 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        heat_values,
        mask=mask_zero,
        cmap="YlGnBu",
        vmin=0,
        vmax=max(0.1, max_jaccard * 1.15),
        annot=annot,
        fmt="",
        linewidths=0.8,
        linecolor="#e5e5e5",
        cbar_kws={"label": "Jaccard similarity"},
        ax=ax,
        annot_kws={"fontsize": 8},
    )

    mapped_pair_lookup = set(zip(mapping_df["fetal_cell_type"], mapping_df["adult_cell_type"])) if not mapping_df.empty else set()
    for i, fetal_ct in enumerate(heat_values.index):
        for j, adult_ct in enumerate(heat_values.columns):
            if (fetal_ct, adult_ct) in mapped_pair_lookup:
                rect = plt.Rectangle((j, i), 1, 1, fill=False, edgecolor="#1a1a1a", linewidth=1.5)
                ax.add_patch(rect)

    ax.set_title("Adult vs fetal marker overlap (mapping-aware, tile text = # shared genes)", fontsize=13, pad=10)
    ax.set_xlabel("Adult cell type", fontsize=11)
    ax.set_ylabel("Fetal cell type", fontsize=11)
    ax.set_xticks(np.arange(len(x_labels)) + 0.5)
    ax.set_xticklabels(x_labels, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(y_labels)) + 0.5)
    ax.set_yticklabels(y_labels, rotation=0, fontsize=8)

    save_figure(fig, out_path, "marker overlap figure", dpi=dpi)


def _yaml_inline_list(items: list) -> str:
    return "[" + ", ".join(json.dumps(x) for x in items) + "]"


def _build_detailed_yaml(
    adult_df,
    fetal_df,
    cfg_overlap,
    adult_src,
    fetal_src,
    comparison_mode: str,
    adult_top_n: int = 25,
) -> str:
    sv = cfg_overlap.get("yaml_schema_version", 1)
    a_map = adult_df.groupby("cell_type_norm")["cell_type"].agg(lambda s: sorted(set(s))).to_dict()
    f_map = fetal_df.groupby("cell_type_norm")["cell_type"].agg(lambda s: sorted(set(s))).to_dict()
    a_ranked = _grouped_ranked_genes(adult_df, top_n=adult_top_n)
    f_ranked = _grouped_ranked_genes(fetal_df, top_n=None)
    lines = [
        f"schema_version: {sv}",
        'description: "Adult/Fetal kidney marker overlap export"',
        "comparison_mode: " + json.dumps(comparison_mode),
        "sources:",
        f"  adult_markers_tsv: {json.dumps(adult_src)}",
        f"  fetal_markers_tsv: {json.dumps(fetal_src)}",
        "cell_types:",
    ]
    for norm in sorted(set(a_ranked) | set(f_ranked)):
        ag = a_ranked.get(norm, [])
        fg = f_ranked.get(norm, [])
        fg_set = set(fg)
        ag_set = set(ag)
        shared = _ordered_intersection(ag, fg)
        lines += [
            "  - cell_type_norm: " + json.dumps(norm),
            "    in_adult: " + ("true" if norm in a_map else "false"),
            "    in_fetal: " + ("true" if norm in f_map else "false"),
            "    adult_labels: " + _yaml_inline_list(a_map.get(norm, [])),
            "    fetal_labels: " + _yaml_inline_list(f_map.get(norm, [])),
            "    adult_markers: " + _yaml_inline_list(ag),
            "    fetal_markers: " + _yaml_inline_list(fg),
            "    shared_markers: " + _yaml_inline_list(shared),
            "    adult_only_markers: " + _yaml_inline_list([g for g in ag if g not in fg_set]),
            "    fetal_only_markers: " + _yaml_inline_list([g for g in fg if g not in ag_set]),
        ]
    return "\n".join(lines) + "\n"


def _build_official_yaml(adult_df, fetal_df, cfg_overlap, comparison_mode: str, adult_top_n: int = 25) -> str:
    sv = cfg_overlap.get("yaml_schema_version", 1)
    a_ranked = _grouped_ranked_genes(adult_df, top_n=adult_top_n)
    f_ranked = _grouped_ranked_genes(fetal_df, top_n=None)
    entries = []
    for norm in sorted(set(a_ranked) | set(f_ranked)):
        overlap_genes = _ordered_intersection(a_ranked.get(norm, []), f_ranked.get(norm, []))
        if overlap_genes:
            entries.append((norm, overlap_genes))
    lines = [
        f"schema_version: {sv}",
        'description: "Official overlap marker YAML (adult vs fetal kidney)"',
        "comparison_mode: " + json.dumps(comparison_mode),
    ]
    if not entries:
        lines.append("cell_types: []")
    else:
        lines.append("cell_types:")
        for norm, overlap_genes in entries:
            lines += [
                "  - cell_type_norm: " + json.dumps(norm),
                "    overlap_genes: " + _yaml_inline_list(overlap_genes),
            ]
    return "\n".join(lines) + "\n"


def _validate_yaml_structure(yaml_text: str, label: str) -> None:
    try:
        data = yaml.safe_load(yaml_text)
    except yaml.YAMLError as e:
        raise ValueError(f"[{label}] YAML parse failure: {e}") from e

    if not isinstance(data, dict):
        raise ValueError(f"[{label}] YAML root must be a mapping.")
    if "schema_version" not in data:
        raise ValueError(f"[{label}] missing schema_version")
    if "cell_types" not in data:
        raise ValueError(f"[{label}] missing cell_types")
    if not isinstance(data["cell_types"], list):
        raise ValueError(f"[{label}] cell_types must be a list")
    if data["cell_types"] and "cell_type_norm" not in data["cell_types"][0]:
        raise ValueError(f"[{label}] first entry missing cell_type_norm")
    log.info("YAML validation passed [%s]: %d cell types", label, len(data["cell_types"]))


def run_phase3(cfg: dict, paths: dict[str, Path],
               phase1_artifacts: dict[str, Path], phase2_artifacts: dict[str, Path]) -> dict[str, Path]:
    require_runtime_packages(
        {
            "matplotlib": plt,
            "numpy": np,
            "pandas": pd,
            "seaborn": sns,
        },
        "Phase 3",
    )

    gate_phase({
        "fetal_markers_all_tsv":    phase1_artifacts["fetal_markers_all_tsv"],
        "adult_top_validated_tsv": phase2_artifacts["adult_top_validated_tsv"],
        "adult_h5ad": phase1_artifacts["adult_h5ad"],
        "fetal_h5ad": phase1_artifacts["fetal_h5ad"],
    }, "Phase 3")

    output_dir  = paths["output_dir"]
    aux_dir     = paths["aux_dir"]
    figures_dir = paths["figures_dir"]
    diagnostics_dir = paths["diagnostics_dir"] / "adult_fetal_overlap"
    ensure_dirs(output_dir, aux_dir, figures_dir, diagnostics_dir)

    cfg_overlap = cfg["overlap"]
    harmonization_lookup = build_harmonization_lookup(cfg_overlap["harmonization_map"])
    manual_mapping_seed = cfg_overlap.get("manual_mapping_seed", _DEFAULT_MANUAL_MAPPING_SEED)
    dpi = cfg.get("umap", {}).get("fig_dpi", 300)
    adult_yaml_top_n = 25
    comparison_mode = (
        f"adult top{adult_yaml_top_n} validated markers (rank-preserving after harmonization) "
        "vs fetal all markers (rank-preserving)"
    )

    adult_df = _load_standardise_markers(
        phase2_artifacts["adult_top_validated_tsv"], "adult",
        ["my_cell_type", "cell_type"], ["gene", "Gene"],
        ["final_rank_within_cell_type", "rank_filtered", "rank_raw", "rank"],
        harmonization_lookup,
    )
    fetal_df = _load_standardise_markers(
        phase1_artifacts["fetal_markers_all_tsv"], "fetal",
        ["cell_type", "my_cell_type"], ["gene", "Gene"],
        ["rank_filtered", "rank_raw", "final_rank_within_cell_type", "rank"],
        harmonization_lookup,
    )
    log.info("Phase 3 inputs — adult: %d rows / %d norm types; fetal: %d rows / %d norm types",
             len(adult_df), adult_df["cell_type_norm"].nunique(),
             len(fetal_df), fetal_df["cell_type_norm"].nunique())

    # Label mapping table (adult/fetal harmonized vocabulary) → aux output
    label_mapping_df = _build_cell_type_label_mapping_table(adult_df, fetal_df)
    label_mapping_tsv = aux_dir / "cell_type_label_mapping.tsv"
    save_tsv(label_mapping_df, label_mapping_tsv, "cell_type_label_mapping")

    adult_summary = _summarize_cell_types_from_h5ad(phase1_artifacts["adult_h5ad"], "adult_preprocessed").set_index("cell_type")
    fetal_summary = _summarize_cell_types_from_h5ad(phase1_artifacts["fetal_h5ad"], "fetal_preprocessed").set_index("cell_type")

    comparison_df, mapping_df, sharing_summary_df, sharing_details_df = _build_mapping_aware_comparison(
        adult_summary,
        fetal_summary,
        set(adult_df["cell_type"].astype(str).unique()),
        set(fetal_df["cell_type"].astype(str).unique()),
        manual_mapping_seed=manual_mapping_seed,
    )

    # Cell-type comparison → main output
    comparison_tsv = output_dir / "cell_type_comparison.tsv"
    save_tsv(comparison_df, comparison_tsv, "cell_type_comparison")

    save_tsv(mapping_df, diagnostics_dir / "cell_type_shared_pairs.tsv", "cell_type shared pairs")
    save_tsv(sharing_summary_df, diagnostics_dir / "cell_type_presence_summary.tsv", "cell_type presence summary")
    save_tsv(sharing_details_df, diagnostics_dir / "cell_type_presence_details.tsv", "cell_type presence details")

    # Marker overlap table → main output
    overlap_df = _build_overlap_table_all_pairs(adult_df, fetal_df, mapping_df)
    overlap_tsv = output_dir / "marker_overlap_adult_vs_fetal.tsv"
    save_tsv(overlap_df, overlap_tsv, "marker_overlap_adult_vs_fetal")

    acceptance_report_df, acceptance_failures = _build_overlap_acceptance_report(
        overlap_df,
        sharing_summary_df,
        cfg_overlap,
    )
    acceptance_report_tsv = diagnostics_dir / "overlap_acceptance_report.tsv"
    save_tsv(acceptance_report_df, acceptance_report_tsv, "overlap acceptance report")
    if acceptance_failures:
        details = "\n  - " + "\n  - ".join(acceptance_failures)
        raise ValueError(
            "Overlap acceptance criteria failed:" + details
            + "\nUpdate overlap.acceptance_criteria or revisit harmonization/marker settings."
        )

    # Overlap figure → main output root
    overlap_fig = output_dir / "marker_overlap_figure.png"
    _plot_overlap_figure_mapping_aware(overlap_df, mapping_df, adult_summary, fetal_summary, overlap_fig, dpi=dpi)

    # SnapSeed YAMLs → aux outputs
    adult_src = str(phase2_artifacts["adult_top_validated_tsv"])
    fetal_src = str(phase1_artifacts["fetal_markers_all_tsv"])

    detailed_text = _build_detailed_yaml(
        adult_df,
        fetal_df,
        cfg_overlap,
        adult_src,
        fetal_src,
        comparison_mode,
        adult_top_n=adult_yaml_top_n,
    )
    official_text = _build_official_yaml(
        adult_df,
        fetal_df,
        cfg_overlap,
        comparison_mode,
        adult_top_n=adult_yaml_top_n,
    )
    _validate_yaml_structure(detailed_text, "detailed")
    _validate_yaml_structure(official_text, "official")

    detailed_yaml = aux_dir / "snapseed_markers_adult_fetal_overlap_detailed.yaml"
    official_yaml = aux_dir / "snapseed_markers_adult_fetal_overlap_official.yaml"
    save_text(detailed_text, detailed_yaml, "snapseed detailed YAML")
    save_text(official_text, official_yaml, "snapseed official YAML")

    log.info("Phase 3 complete.")
    return {
        "cell_type_comparison_tsv": comparison_tsv,
        "cell_type_label_mapping_tsv": label_mapping_tsv,
        "marker_overlap_tsv":       overlap_tsv,
        "overlap_acceptance_report_tsv": acceptance_report_tsv,
        "overlap_figure":           overlap_fig,
        "detailed_yaml":            detailed_yaml,
        "official_yaml":            official_yaml,
        "overlap_diagnostics_dir":  diagnostics_dir,
    }


# =============================================================================
# SECTION 6 — Dry-run & Verify
# =============================================================================

def dry_run(cfg: dict, paths: dict[str, Path]) -> None:
    log.info("=== DRY-RUN: validating config and inputs ===")
    require_file(paths["adult_h5ad"], "adult h5ad")
    require_file(paths["fetal_h5ad"], "fetal h5ad")
    require_file(paths["hpa_tsv"],    "HPA reference TSV")
    log.info("All required input files exist.")
    for name in ("output_dir", "aux_dir", "figures_dir", "diagnostics_dir"):
        p = paths[name]
        try:
            p.mkdir(parents=True, exist_ok=True)
            log.info("Output directory writable: %s", p)
        except OSError as e:
            raise RuntimeError(f"Cannot create output directory [{name}]: {p}\n{e}") from e
    log.info("Seed: %d | Top-N markers: %d | Top-N validated: %d",
             cfg.get("seed", 42), cfg["marker_extraction"]["top_n_markers"],
             cfg["hpa_validation"]["top_n_validated"])
    log.info("Dry-run passed — ready to execute.")


def verify_outputs(cfg: dict, paths: dict[str, Path]) -> None:
    output_dir  = paths["output_dir"]
    aux_dir     = paths["aux_dir"]
    figures_dir = paths["figures_dir"]
    diagnostics_dir = paths["diagnostics_dir"]
    top_n_val = cfg["hpa_validation"]["top_n_validated"]

    checks: list[tuple[str, Path]] = [
        # Run metadata outputs
        ("aux/dataset_manifest.tsv",                     aux_dir / "dataset_manifest.tsv"),
        ("aux/marker_method_params.yaml",                aux_dir / "marker_method_params.yaml"),
        # Phase 1
        (f"fetal_kidney_markers_top{top_n_val}.tsv",      output_dir / f"fetal_kidney_markers_top{top_n_val}.tsv"),
        ("aux/adult_kidney_markers_all.tsv",               aux_dir / "adult_kidney_markers_all.tsv"),
        ("aux/fetal_kidney_markers_all.tsv",               aux_dir / "fetal_kidney_markers_all.tsv"),
        ("aux/fetal_generic_kidney_rescue_assignments.tsv",aux_dir / "fetal_generic_kidney_rescue_assignments.tsv"),
        ("aux/adult_preprocessed.h5ad",                    aux_dir / "adult_preprocessed.h5ad"),
        ("aux/fetal_preprocessed.h5ad",                    aux_dir / "fetal_preprocessed.h5ad"),
        ("aux/figures/umap_adult_fetal_two_panel.png",     figures_dir / "umap_adult_fetal_two_panel.png"),
        # Phase 2
        (f"adult_kidney_markers_top{top_n_val}_validated.tsv",
         output_dir / f"adult_kidney_markers_top{top_n_val}_validated.tsv"),
        (f"aux/adult_kidney_markers_top{top_n_val}_validated.tsv",
         aux_dir / f"adult_kidney_markers_top{top_n_val}_validated.tsv"),
        ("aux/diagnostics/adult_hpa_validation/adult_hpa_validation_audit.tsv",
         diagnostics_dir / "adult_hpa_validation" / "adult_hpa_validation_audit.tsv"),
        ("aux/diagnostics/adult_hpa_validation/adult_hpa_ranked_candidates.tsv",
         diagnostics_dir / "adult_hpa_validation" / "adult_hpa_ranked_candidates.tsv"),
        # Phase 3
        ("cell_type_comparison.tsv",               output_dir / "cell_type_comparison.tsv"),
        ("aux/cell_type_label_mapping.tsv",        aux_dir / "cell_type_label_mapping.tsv"),
        ("marker_overlap_adult_vs_fetal.tsv",      output_dir / "marker_overlap_adult_vs_fetal.tsv"),
        ("marker_overlap_figure.png",              output_dir / "marker_overlap_figure.png"),
        ("aux/diagnostics/adult_fetal_overlap/cell_type_shared_pairs.tsv",
         diagnostics_dir / "adult_fetal_overlap" / "cell_type_shared_pairs.tsv"),
        ("aux/diagnostics/adult_fetal_overlap/cell_type_presence_summary.tsv",
         diagnostics_dir / "adult_fetal_overlap" / "cell_type_presence_summary.tsv"),
        ("aux/diagnostics/adult_fetal_overlap/cell_type_presence_details.tsv",
         diagnostics_dir / "adult_fetal_overlap" / "cell_type_presence_details.tsv"),
        ("aux/diagnostics/adult_fetal_overlap/overlap_acceptance_report.tsv",
         diagnostics_dir / "adult_fetal_overlap" / "overlap_acceptance_report.tsv"),
        ("aux/snapseed_...detailed.yaml",          aux_dir / "snapseed_markers_adult_fetal_overlap_detailed.yaml"),
        ("aux/snapseed_...official.yaml",          aux_dir / "snapseed_markers_adult_fetal_overlap_official.yaml"),
    ]

    passed, failed = 0, 0
    log.info("=== OUTPUT CONTRACT VERIFICATION ===")
    for label, path in checks:
        if path.exists() and path.stat().st_size > 0:
            log.info("  PASS  %s", label); passed += 1
        else:
            log.warning("  FAIL  %s  (%s)", label, "not found" if not path.exists() else "empty"); failed += 1

    # YAML structural checks
    for yaml_label, yaml_path in [
        ("detailed YAML", aux_dir / "snapseed_markers_adult_fetal_overlap_detailed.yaml"),
        ("official YAML", aux_dir / "snapseed_markers_adult_fetal_overlap_official.yaml"),
    ]:
        if yaml_path.exists():
            try:
                data = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
                if not isinstance(data, dict):
                    raise ValueError("YAML root must be a mapping.")
                if "schema_version" not in data or "cell_types" not in data:
                    raise ValueError("Missing required YAML keys: schema_version/cell_types.")
                if not isinstance(data["cell_types"], list):
                    raise ValueError("YAML key 'cell_types' must be a list.")
                log.info("  PASS  %s parses correctly (%d cell types)", yaml_label, len(data.get("cell_types", [])))
                passed += 1
            except Exception as e:
                log.warning("  FAIL  %s parse error: %s", yaml_label, e); failed += 1
        else:
            log.warning("  FAIL  %s not found", yaml_label); failed += 1

    # Legacy notebook artifacts must NOT exist
    for p in [aux_dir / "sectionA_top20_validated_markers.tsv",
              aux_dir / "sectionA_top20_validation_audit.tsv",
              aux_dir / "sectionA_top20_validation_stats.tsv"]:
        if p.exists():
            log.warning("  FAIL  Legacy artifact must not exist: %s", p); failed += 1
        else:
            log.info("  PASS  Legacy artifact absent (correct): %s", p.name); passed += 1

    log.info("Verification: %d passed, %d failed.", passed, failed)
    sys.exit(0 if failed == 0 else 1)


# =============================================================================
# SECTION 7 — CLI
# =============================================================================

def _reconstruct_phase_artifacts(cfg: dict, paths: dict[str, Path]) -> tuple[dict, dict]:
    """Build artifact path dicts from expected disk locations for partial runs."""
    top_n_val = cfg["hpa_validation"]["top_n_validated"]
    output_dir, aux_dir, figures_dir = paths["output_dir"], paths["aux_dir"], paths["figures_dir"]

    phase1 = {
        "adult_markers_all_tsv": aux_dir / "adult_kidney_markers_all.tsv",
        "fetal_markers_all_tsv": aux_dir / "fetal_kidney_markers_all.tsv",
        "fetal_markers_top_tsv": output_dir / f"fetal_kidney_markers_top{top_n_val}.tsv",
        "fetal_rescue_tsv":      aux_dir / "fetal_generic_kidney_rescue_assignments.tsv",
        "adult_h5ad":            aux_dir / "adult_preprocessed.h5ad",
        "fetal_h5ad":            aux_dir / "fetal_preprocessed.h5ad",
        "umap_figure":           figures_dir / "umap_adult_fetal_two_panel.png",
    }
    phase2 = {
        "adult_top_validated_tsv": aux_dir / f"adult_kidney_markers_top{top_n_val}_validated.tsv",
        "adult_top_validated_main_tsv": output_dir / f"adult_kidney_markers_top{top_n_val}_validated.tsv",
    }
    return phase1, phase2


def _reconstruct_phase3_artifacts(paths: dict[str, Path]) -> dict[str, Path]:
    output_dir, aux_dir = paths["output_dir"], paths["aux_dir"]
    diagnostics_dir = paths["diagnostics_dir"] / "adult_fetal_overlap"
    return {
        "cell_type_comparison_tsv": output_dir / "cell_type_comparison.tsv",
        "cell_type_label_mapping_tsv": aux_dir / "cell_type_label_mapping.tsv",
        "marker_overlap_tsv":       output_dir / "marker_overlap_adult_vs_fetal.tsv",
        "overlap_acceptance_report_tsv": diagnostics_dir / "overlap_acceptance_report.tsv",
        "overlap_figure":           output_dir / "marker_overlap_figure.png",
        "detailed_yaml":            aux_dir / "snapseed_markers_adult_fetal_overlap_detailed.yaml",
        "official_yaml":            aux_dir / "snapseed_markers_adult_fetal_overlap_official.yaml",
        "overlap_diagnostics_dir":  diagnostics_dir,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Unified Adult-Fetal Kidney Reference Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--config", "-c", default="config.yaml",
                        help="Path to config.yaml (default: ./config.yaml)")
    parser.add_argument("--phases", "--only-phases", "-p", nargs="+", type=int, choices=[1, 2, 3],
                        default=[1, 2, 3], dest="phases",
                        help="Only run selected phases (default: 1 2 3)")
    parser.add_argument("--skip-phases", nargs="+", type=int, choices=[1, 2, 3], default=[],
                        help="Skip specific phases from the selected set")
    parser.add_argument("--log-file", default=None,
                        help="Optional log file path (default: <aux_dir>/logs/pipeline_<run_id>.log)")
    parser.add_argument("--no-resume", action="store_true",
                        help="Disable .meta.json resume checks and force selected phases to rerun")
    parser.add_argument("--dry-run",  action="store_true", help="Validate inputs without processing")
    parser.add_argument("--verify",   action="store_true", help="Verify output contract")
    return parser.parse_args()


def emit_run_metadata_outputs(cfg: dict[str, Any], paths: dict[str, Path]) -> dict[str, Path]:
    """Emit run-level metadata outputs that are independent of individual phases."""
    ensure_dirs(paths["output_dir"])
    out: dict[str, Path] = {}

    dataset_manifest = write_dataset_manifest(cfg, paths)
    if dataset_manifest is not None:
        out["dataset_manifest_tsv"] = dataset_manifest

    marker_params = write_marker_method_params(cfg, paths)
    out["marker_method_params_yaml"] = marker_params
    return out


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    cfg = load_config(config_path)
    paths = build_paths(cfg, config_path.parent)

    run_id = _new_run_id()
    default_log_file = paths["aux_dir"] / "logs" / f"pipeline_{run_id}.log"
    log_file = Path(args.log_file).resolve() if args.log_file else default_log_file
    add_file_logging(log_file)

    selected_phases = set(args.phases)
    skipped_phases = set(args.skip_phases)
    phases = sorted(selected_phases - skipped_phases)
    if not phases:
        raise ValueError("No phases selected after applying --skip-phases. Please select at least one phase.")

    write_run_manifest(
        manifest_path=paths["aux_dir"] / "run_manifests" / f"run_manifest_{run_id}.json",
        run_id=run_id,
        config_path=config_path,
        cfg=cfg,
        paths=paths,
        phases=phases,
        args=args,
    )

    if args.dry_run:
        dry_run(cfg, paths); return
    if args.verify:
        verify_outputs(cfg, paths); return

    run_metadata_artifacts = emit_run_metadata_outputs(cfg, paths)

    log.info("=== Kidney Reference Pipeline — phases %s ===", phases)
    t_start = time.perf_counter()

    resume_meta_path = paths["aux_dir"] / ".meta.json"
    resume_meta = _load_resume_meta(resume_meta_path)
    if args.no_resume:
        log.info("Resume disabled (--no-resume): selected phases will run even if prior artifacts exist.")
    else:
        log.info("Resume metadata path: %s", resume_meta_path)

    script_sha256 = _sha256_file(Path(__file__).resolve())
    phase1_fingerprint = _build_phase1_fingerprint(cfg, paths, script_sha256)
    phase2_fingerprint = _build_phase2_fingerprint(cfg, paths, phase1_fingerprint, script_sha256)
    phase3_fingerprint = _build_phase3_fingerprint(cfg, phase1_fingerprint, phase2_fingerprint, script_sha256)

    phase1_artifacts, phase2_artifacts = _reconstruct_phase_artifacts(cfg, paths)
    phase3_expected_artifacts = _reconstruct_phase3_artifacts(paths)
    phase3_artifacts: dict[str, Path] = {}

    if 2 in phases and 1 not in phases:
        _ensure_upstream_compatibility(resume_meta, 1, phase1_fingerprint, "Phase 2")
    if 3 in phases and 1 not in phases:
        _ensure_upstream_compatibility(resume_meta, 1, phase1_fingerprint, "Phase 3")
    if 3 in phases and 2 not in phases:
        _ensure_upstream_compatibility(resume_meta, 2, phase2_fingerprint, "Phase 3")

    if 1 in phases:
        run_phase = True
        if args.no_resume:
            log.info("--- Phase 1: Preprocessing + Marker Extraction (forced run) ---")
        else:
            resume_ok, reason = _can_resume_phase(1, phase1_fingerprint, phase1_artifacts, resume_meta)
            if resume_ok:
                run_phase = False
                log.info("--- Phase 1: skipped via resume (%s) ---", reason)
            else:
                log.info("--- Phase 1: Preprocessing + Marker Extraction ---")
                if "fingerprint changed" in reason and _has_any_existing_outputs(phase1_artifacts):
                    log.info("Resume invalidated for Phase 1 (%s). Existing artifacts will be overwritten.", reason)
                else:
                    log.info("Resume miss for Phase 1 (%s). Running phase.", reason)

        if run_phase:
            t1 = time.perf_counter()
            phase1_artifacts = run_phase1(cfg, paths)
            _record_phase_resume_state(resume_meta, 1, phase1_fingerprint, phase1_artifacts, run_id)
            _save_resume_meta(resume_meta_path, resume_meta)
            log.info("Phase 1 done in %.1f s", time.perf_counter() - t1)

    if 2 in phases:
        run_phase = True
        if args.no_resume:
            log.info("--- Phase 2: Adult HPA Validation + Ranking (forced run) ---")
        else:
            resume_ok, reason = _can_resume_phase(2, phase2_fingerprint, phase2_artifacts, resume_meta)
            if resume_ok:
                run_phase = False
                log.info("--- Phase 2: skipped via resume (%s) ---", reason)
            else:
                log.info("--- Phase 2: Adult HPA Validation + Ranking ---")
                if "fingerprint changed" in reason and _has_any_existing_outputs(phase2_artifacts):
                    log.info("Resume invalidated for Phase 2 (%s). Existing artifacts will be overwritten.", reason)
                else:
                    log.info("Resume miss for Phase 2 (%s). Running phase.", reason)

        if run_phase:
            t2 = time.perf_counter()
            phase2_artifacts = run_phase2(cfg, paths, phase1_artifacts)
            _record_phase_resume_state(resume_meta, 2, phase2_fingerprint, phase2_artifacts, run_id)
            _save_resume_meta(resume_meta_path, resume_meta)
            log.info("Phase 2 done in %.1f s", time.perf_counter() - t2)

    if 3 in phases:
        run_phase = True
        if args.no_resume:
            log.info("--- Phase 3: Overlap Analysis + YAML Export (forced run) ---")
        else:
            resume_ok, reason = _can_resume_phase(3, phase3_fingerprint, phase3_expected_artifacts, resume_meta)
            if resume_ok:
                run_phase = False
                phase3_artifacts = phase3_expected_artifacts
                log.info("--- Phase 3: skipped via resume (%s) ---", reason)
            else:
                log.info("--- Phase 3: Overlap Analysis + YAML Export ---")
                if "fingerprint changed" in reason and _has_any_existing_outputs(phase3_expected_artifacts):
                    log.info("Resume invalidated for Phase 3 (%s). Existing artifacts will be overwritten.", reason)
                else:
                    log.info("Resume miss for Phase 3 (%s). Running phase.", reason)

        if run_phase:
            t3 = time.perf_counter()
            phase3_artifacts = run_phase3(cfg, paths, phase1_artifacts, phase2_artifacts)
            _record_phase_resume_state(resume_meta, 3, phase3_fingerprint, phase3_artifacts, run_id)
            _save_resume_meta(resume_meta_path, resume_meta)
            log.info("Phase 3 done in %.1f s", time.perf_counter() - t3)

    log.info("=== Pipeline complete in %.1f s ===", time.perf_counter() - t_start)
    log.info("--- Output summary ---")
    for name, path in sorted({**run_metadata_artifacts, **phase1_artifacts, **phase2_artifacts, **phase3_artifacts}.items()):
        log.info("  [%s]  %-45s  %s", "OK" if path.exists() else "MISSING", name, path)


if __name__ == "__main__":
    main()
