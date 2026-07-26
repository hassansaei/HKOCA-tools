"""Load and resolve CellBender configuration."""

from __future__ import annotations

import configparser
import os
from dataclasses import dataclass
from typing import Sequence

CONFIG_FILENAME = "cellbender.config"


@dataclass(frozen=True)
class CellBenderParams:
    mode: str = "full"
    cuda: bool = True
    epochs: int = 100
    posterior_batch_size: int = 64
    total_droplets_included: int = 35000
    learning_rate: float = 0.00005
    cpu_threads: int = 24


@dataclass(frozen=True)
class CellBenderConfig:
    working_dir: str
    samples_dir: str
    samples: tuple[str, ...]
    h5_filename: str
    mtx_dirname: str
    output_suffix: str
    params: CellBenderParams
    config_path: str


def default_config_path(cli_config: str | None = None) -> str:
    if cli_config and cli_config.strip():
        return cli_config.strip()
    cwd = os.path.join(os.getcwd(), CONFIG_FILENAME)
    if os.path.isfile(cwd):
        return cwd
    pkg_cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), CONFIG_FILENAME)
    if os.path.isfile(pkg_cfg):
        return pkg_cfg
    return cwd


def _as_bool(value: str, default: bool = False) -> bool:
    if value is None or str(value).strip() == "":
        return default
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _parse_samples(raw: str) -> tuple[str, ...]:
    if not raw or not raw.strip():
        return ()
    return tuple(s.strip() for s in raw.split(",") if s.strip())


def load_config(config_path: str) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    if os.path.isfile(config_path):
        cfg.read(config_path, encoding="utf-8")
    return cfg


def resolve_config(
    config_path: str,
    *,
    working_dir: str | None = None,
    samples_dir: str | None = None,
    samples: Sequence[str] | None = None,
    epochs: int | None = None,
    total_droplets_included: int | None = None,
    learning_rate: float | None = None,
    cpu_threads: int | None = None,
    cuda: bool | None = None,
) -> CellBenderConfig:
    cfg = load_config(config_path)

    def get(section: str, key: str, default: str = "") -> str:
        if cfg.has_option(section, key):
            return cfg.get(section, key).strip()
        return default

    wd = working_dir or get("paths", "working_dir") or os.getcwd()
    if not os.path.isabs(wd):
        wd = os.path.abspath(wd)

    sd = samples_dir if samples_dir is not None else get("paths", "samples_dir")
    if sd and not os.path.isabs(sd):
        sd = os.path.normpath(os.path.join(wd, sd))

    sample_list: tuple[str, ...]
    if samples:
        sample_list = tuple(samples)
    else:
        sample_list = _parse_samples(get("paths", "samples"))

    params = CellBenderParams(
        mode=get("cellbender", "mode", "full") or "full",
        cuda=_as_bool(get("cellbender", "cuda", "true"), True) if cuda is None else cuda,
        epochs=int(epochs if epochs is not None else get("cellbender", "epochs", "100") or 100),
        posterior_batch_size=int(get("cellbender", "posterior_batch_size", "64") or 64),
        total_droplets_included=int(
            total_droplets_included
            if total_droplets_included is not None
            else get("cellbender", "total_droplets_included", "35000") or 35000
        ),
        learning_rate=float(
            learning_rate
            if learning_rate is not None
            else get("cellbender", "learning_rate", "0.00005") or 0.00005
        ),
        cpu_threads=int(
            cpu_threads
            if cpu_threads is not None
            else get("cellbender", "cpu_threads", "24") or 24
        ),
    )

    return CellBenderConfig(
        working_dir=wd,
        samples_dir=sd or "",
        samples=sample_list,
        h5_filename=get("input", "h5_filename", "raw_feature_bc_matrix.h5")
        or "raw_feature_bc_matrix.h5",
        mtx_dirname=get("input", "mtx_dirname", "raw_feature_bc_matrix")
        or "raw_feature_bc_matrix",
        output_suffix=get("input", "output_suffix", "_filtered.h5") or "_filtered.h5",
        params=params,
        config_path=config_path,
    )


def resolve_sample_dir(cfg: CellBenderConfig, sample: str) -> str:
    """Return absolute path to a sample directory."""
    if os.path.isabs(sample):
        return sample
    if cfg.samples_dir:
        return os.path.normpath(os.path.join(cfg.samples_dir, sample))
    return os.path.normpath(os.path.join(cfg.working_dir, sample))


def sample_id_from_path(sample_dir: str) -> str:
    return os.path.basename(os.path.normpath(sample_dir))
