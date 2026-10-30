process QC {
    tag { "qc" }
    label 'process_high'

    conda "${params.conda_qc}"

    publishDir "${params.outdir}/qc_filter", mode: params.publish_mode

    input:
    path harmonize_root

    output:
    path "qc_out", emit: output_dir
    path "qc_out/qc_filtered_rds/*_filtered.rds", emit: filtered_rds
    path "qc_out/h5ad_converted/*_filtered.h5ad", emit: filtered_h5ad, optional: true
    path "versions.yml", emit: versions

    script:
    def qc_config = params.qc_config ?: ''
    """
    set -euo pipefail
    export HKOCA_ROOT="${projectDir}"
    export PYTHONPATH="${projectDir}\${PYTHONPATH:+:\$PYTHONPATH}"
    export HKOCA_QC_ENV="${params.qc_env_name}"
    export HKOCA_QC_RDS_DIR="${harmonize_root}"
    export HKOCA_QC_OUT="qc_out"
    export HKOCA_QC_RUN="${params.qc_run}"
    export HKOCA_QC_RDS_PATTERN='${params.qc_rds_pattern}'
    export HKOCA_QC_FORCE="${params.force ? '1' : '0'}"
    export HKOCA_QC_CONFIG="${qc_config}"

    mkdir -p qc_out

    python << 'PY'
    import os
    import sys
    from hkoca.qc_filter.qc_runner import run_qc

    cfg = os.environ.get("HKOCA_QC_CONFIG", "").strip() or None
    rc = run_qc(
        rds_dir=os.environ["HKOCA_QC_RDS_DIR"],
        output_dir=os.environ["HKOCA_QC_OUT"],
        config=cfg,
        run=os.environ.get("HKOCA_QC_RUN", "all"),
        recursive=True,
        rds_pattern=os.environ["HKOCA_QC_RDS_PATTERN"],
        force_overwrite=os.environ.get("HKOCA_QC_FORCE", "0") == "1",
        dry_run=False,
    )
    sys.exit(rc)
    PY

    test -d qc_out/qc_filtered_rds
    ls qc_out/qc_filtered_rds/*_filtered.rds >/dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hkoca: \$(python -c "from hkoca._version import __version__; print(__version__)" 2>/dev/null || echo "unknown")
        rscript: \$(Rscript --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
