# syntax=docker/dockerfile:1.7
# =============================================================================
# HKOCA / NephAtlas toolkit image
# =============================================================================
# Builds stage-specific conda environments from conda/*.yaml, then installs the
# hkoca package so the CLI is available inside the container.
#
# Build (from repository root). Pin linux/amd64 for conda/bioconda + CellBender
# CUDA wheels; required on Apple Silicon hosts.
#   docker build --platform linux/amd64 -f docker/Dockerfile -t hkoca:local .
#
# Run (match the build platform):
#   docker run --rm --platform linux/amd64 hkoca:local --version
#   docker run --rm --platform linux/amd64 \
#       -v "$PWD/data:/data" -v "$PWD/out:/out" hkoca:local \
#       qc-filter --csv /data/meta.csv --gtf /data/genes.gtf --output /out
#
# Stage environments (set HKOCA_ENV to switch; default is hkoca_harmonize):
#   hkoca_harmonize   - gene-space harmonization / default CLI (incl. qc-filter)
#   hkoca_cellbender  - CellBender ambient RNA removal (GPU optional)
#   sc_qc_pipeline    - QC R stack; qc-filter hands off to this env automatically
#
#   docker run --rm --platform linux/amd64 -e HKOCA_ENV=hkoca_cellbender \
#       --gpus all hkoca:local \
#       cellbender h5 --samples-dir /data/samples --samples sampleA
#
# NOTE ON CACHE MOUNTS: deliberately unused. BuildKit does not export
# `RUN --mount=type=cache` into remote cache backends; mounting a cache at
# /opt/conda/pkgs also breaks conda hardlinking across filesystems and can
# inflate the image.
# =============================================================================

ARG TARGETPLATFORM=linux/amd64
ARG MINIFORGE=condaforge/miniforge3:24.11.3-0@sha256:bcfc04b3562faa8dfd48f3767e2bd7b73a8529f9611b508d42e5bd227db1342d

# =============================================================================
# Stage: envs - conda environments from conda/
# =============================================================================
FROM --platform=$TARGETPLATFORM ${MINIFORGE} AS envs

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# miniforge 24.11.3-0 already ships mamba with libmamba. Do not reinstall mamba
# from conda-forge here; that floats the toolchain and breaks the digest pin.

COPY conda/ /tmp/conda/
RUN mamba env create -f /tmp/conda/environment_harmonize.yaml && \
    mamba env create -f /tmp/conda/environment_cellbender.yaml && \
    mamba env create -f /tmp/conda/environment_qc.yaml && \
    conda clean -afy && \
    find /opt/conda/envs -type d -name '__pycache__' -prune -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda/envs -type f \( -name '*.pyc' -o -name '*.pyo' -o -name '*.a' \) -delete 2>/dev/null || true && \
    rm -rf /tmp/conda

# =============================================================================
# Stage: runtime - slim image with envs + hkoca package
# =============================================================================
FROM --platform=$TARGETPLATFORM ${MINIFORGE} AS runtime

LABEL org.opencontainers.image.title="HKOCA" \
      org.opencontainers.image.description="NephAtlas / HKOCA modular single-cell kidney organoid analysis toolkit" \
      org.opencontainers.image.vendor="HKOCA / NephAtlas" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/hassansaei/HKOCA-tools"

ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    HKOCA_ENV=hkoca_harmonize \
    PATH=/opt/conda/envs/hkoca_harmonize/bin:$PATH \
    DEFAULT_INPUT_DIR=/opt/hkoca/input \
    DEFAULT_OUTPUT_DIR=/opt/hkoca/output

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ARG USERNAME=appuser
ARG USER_UID=1001
ARG USER_GID=1001

RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID --shell /bin/bash --create-home $USERNAME && \
    install -d -o $USERNAME -g $USERNAME \
        $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR /data /out

COPY --from=envs /opt/conda/envs /opt/conda/envs

WORKDIR /opt/hkoca

# Application source from the build context (not a git clone), so the image
# matches the commit being built. --no-deps: stage deps live in conda envs.
COPY pyproject.toml README.md ./
COPY hkoca/ ./hkoca/
RUN for env in hkoca_harmonize hkoca_cellbender sc_qc_pipeline; do \
        conda run -n "$env" pip install --no-cache-dir --no-deps . ; \
    done && \
    conda run -n hkoca_harmonize hkoca --version && \
    rm -rf /opt/hkoca/hkoca /opt/hkoca/pyproject.toml /opt/hkoca/README.md \
           /opt/hkoca/*.egg-info /opt/hkoca/build 2>/dev/null || true

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

USER appuser

WORKDIR /data

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]
