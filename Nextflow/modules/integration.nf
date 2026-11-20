process INTEGRATION {
    tag { filtered_rds.simpleName }
    label 'process_high'
    label 'process_long'

    conda "${params.conda_integration}"

    publishDir path: { "${params.outdir}/integration/${filtered_rds.simpleName}" }, mode: params.publish_mode

    input:
    tuple path(filtered_rds), path(annotated_h5ad)

    output:
    tuple path(filtered_rds), path("integration_out/prep/sct_prepared.rds"), emit: prepared
    path "integration_out", emit: output_dir
    path "integration_out/objects/integrated_*.rds", emit: integrated_rds, optional: true
    path "versions.yml", emit: versions

    script:
    def config_arg   = params.integration_config ? "--config ${params.integration_config}" : ""
    def methods_arg  = params.integration_methods ? "--methods ${params.integration_methods}" : ""
    def force_flag   = params.force ? "--force" : ""
    def ann_arg      = annotated_h5ad.name != 'NO_FILE' ? "--annotated-h5ad ${annotated_h5ad}" : ""
    """
    set -euo pipefail
    export HKOCA_ROOT="${projectDir}"
    export PYTHONPATH="${projectDir}\${PYTHONPATH:+:\$PYTHONPATH}"
    export R_MAX_VSIZE="${params.r_max_vsize}"
    export RETICULATE_PYTHON="\$(command -v python)"
    export HKOCA_ANNOTATION_PYTHON="\${RETICULATE_PYTHON}"

    mkdir -p integration_out

    python -m hkoca.integration prep \\
        --input-rds "${filtered_rds}" \\
        --output-dir "\$(pwd)/integration_out" \\
        ${ann_arg} \\
        ${config_arg} \\
        ${force_flag}

    python -m hkoca.integration run \\
        --prepared-rds "\$(pwd)/integration_out/prep/sct_prepared.rds" \\
        --output-dir "\$(pwd)/integration_out" \\
        ${methods_arg} \\
        ${config_arg} \\
        ${force_flag}

    test -f integration_out/prep/sct_prepared.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hkoca: \$(python -c "from hkoca._version import __version__; print(__version__)" 2>/dev/null || echo "unknown")
        rscript: \$(Rscript --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
