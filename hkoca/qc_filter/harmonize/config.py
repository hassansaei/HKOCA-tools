"""Harmonize config loading (stdlib only; safe without scanpy)."""

from __future__ import annotations

import configparser
import os

CONFIG_FILENAME = "harmonize.config"


def config_path(cli_config: str | None = None) -> str:
    if cli_config and cli_config.strip():
        return cli_config.strip()
    cwd_cfg = os.path.join(os.getcwd(), CONFIG_FILENAME)
    if os.path.isfile(cwd_cfg):
        return cwd_cfg
    pkg_cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), CONFIG_FILENAME)
    if os.path.isfile(pkg_cfg):
        return pkg_cfg
    return cwd_cfg


def load_config(path: str) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    if os.path.isfile(path):
        cfg.read(path, encoding="utf-8")
    return cfg


def resolve_paths(args, cfg: configparser.ConfigParser) -> dict:
    """CLI > config [paths] > environment variables."""

    def get_path(cli_val, cfg_key: str, env_key: str, default: str = "") -> str:
        if cli_val is not None and str(cli_val).strip():
            return str(cli_val).strip()
        if cfg.has_section("paths") and cfg.has_option("paths", cfg_key):
            v = cfg.get("paths", cfg_key).strip()
            if v:
                return v
        v = os.environ.get(env_key, "").strip()
        if v:
            return v
        return default

    working_dir = get_path(
        getattr(args, "working_dir", None), "working_dir", "WORKING_DIR", os.getcwd()
    )
    if not os.path.isabs(working_dir):
        working_dir = os.path.abspath(working_dir)

    return {
        "gtf_file": get_path(getattr(args, "gtf", None), "gtf_file", "GTF_FILE"),
        "metadata_csv": get_path(getattr(args, "csv", None), "metadata_csv", "METADATA_CSV"),
        "output_root": get_path(getattr(args, "output", None), "output_root", "OUTPUT_ROOT"),
        "working_dir": working_dir,
    }


def get_summary_options(cfg: configparser.ConfigParser) -> dict:
    defaults = {
        "figure_dpi": 300,
        "figure_extensions": ["png", "pdf"],
        "report_subdir": "reports/atlas_summary",
        "age_plot_top_n": 15,
    }
    if not cfg.has_section("summary"):
        return defaults
    out = dict(defaults)
    if cfg.has_option("summary", "figure_dpi"):
        try:
            _v = cfg.getint("summary", "figure_dpi")
            out["figure_dpi"] = _v if _v > 0 else defaults["figure_dpi"]
        except ValueError:
            pass
    if cfg.has_option("summary", "figure_extensions"):
        raw = cfg.get("summary", "figure_extensions").strip()
        if raw:
            parsed = [x.strip() for x in raw.split(",") if x.strip()]
            out["figure_extensions"] = parsed or defaults["figure_extensions"]
    if cfg.has_option("summary", "report_subdir"):
        v = cfg.get("summary", "report_subdir").strip()
        if v:
            out["report_subdir"] = v
    if cfg.has_option("summary", "age_plot_top_n"):
        try:
            _v = cfg.getint("summary", "age_plot_top_n")
            out["age_plot_top_n"] = _v if _v > 0 else defaults["age_plot_top_n"]
        except ValueError:
            pass
    return out


def resolve_transgenes(args, cfg: configparser.ConfigParser) -> set[str]:
    if getattr(args, "transgenes", None):
        return {t.strip() for t in args.transgenes.split(",") if t.strip()}
    if cfg.has_section("transgenes") and cfg.has_option("transgenes", "names"):
        raw = cfg.get("transgenes", "names").strip()
        return {t.strip() for t in raw.split(",") if t.strip()}
    return set()
