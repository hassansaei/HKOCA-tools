# syntax=docker/dockerfile:1.7
# HKOCA tool image.
#
# Installs the hkoca package onto hkoca.base. Build the base first, and again
# only when the conda environments change:
#   docker build --platform linux/amd64 -f docker/Dockerfile.base -t hkoca.base:latest .
#   docker build --platform linux/amd64 -f docker/Dockerfile -t hkoca:local .
#
# A published base can be named explicitly:
#   docker build --platform linux/amd64 -f docker/Dockerfile \
#       --build-arg BASE_IMAGE=ghcr.io/hassansaei/hkoca.base:latest \
#       -t hkoca:local .
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
#   hkoca_integration - Seurat integration (SCT, Harmony, clustree, ...)
#   hkoca_projection  - scPoli atlas projection (PyTorch + scArches; GPU)
#
#   docker run --rm --platform linux/amd64 -e HKOCA_ENV=hkoca_cellbender \
#       --gpus all hkoca:local \
#       cellbender h5 --samples-dir /data/samples --samples sampleA

ARG BASE_IMAGE=hkoca.base:latest
ARG TARGETPLATFORM=linux/amd64

FROM --platform=$TARGETPLATFORM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="HKOCA" \
      org.opencontainers.image.description="NephAtlas / HKOCA modular single-cell kidney organoid analysis toolkit" \
      org.opencontainers.image.vendor="HKOCA / NephAtlas" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.source="https://github.com/hassansaei/HKOCA-tools"

ARG HKOCA_VERSION=dev
LABEL org.opencontainers.image.version="${HKOCA_VERSION}"

USER root
WORKDIR /opt/hkoca

# Keep the source tree at /opt/hkoca so sibling conda envs can import hkoca via
# HKOCA_ROOT/PYTHONPATH without installing the package into each env.
# --no-deps: stage deps live in the base conda envs; the CLI is installed
# only in hkoca_harmonize.
COPY pyproject.toml README.md ./
COPY hkoca/ ./hkoca/
ENV HKOCA_ROOT=/opt/hkoca
RUN conda run -n hkoca_harmonize pip install --no-cache-dir --no-deps . && \
    conda run -n hkoca_harmonize hkoca --version && \
    conda clean -afy && \
    rm -rf /opt/hkoca/*.egg-info /opt/hkoca/build 2>/dev/null || true

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

USER appuser
WORKDIR /data

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]
