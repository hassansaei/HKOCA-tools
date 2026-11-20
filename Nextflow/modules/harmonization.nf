process HARMONIZATION {
    tag { "harmonize" }
    label 'process_high'

    conda "${params.conda_harmonize}"

    publishDir "${params.outdir}/harmonized", mode: params.publish_mode

    input:
    path harmonize_csv
    path gtf

    output:
    path "harmonize_out", emit: output_root
    path "harmonize_out/rds/*_harmonized.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    def config_arg     = params.harmonize_config ? "--config ${params.harmonize_config}" : ""
    def transgene_arg  = params.transgenes ? "--transgenes ${params.transgenes}" : ""
    def force_flag     = params.force ? "--force" : "--skip-existing"
    """
    set -euo pipefail
    export HKOCA_ROOT="${projectDir}"
    export PYTHONPATH="${projectDir}\${PYTHONPATH:+:\$PYTHONPATH}"

    mkdir -p harmonize_out

    python -m hkoca.qc_filter.harmonize \\
        --csv "${harmonize_csv}" \\
        --gtf "${gtf}" \\
        --output "\$(pwd)/harmonize_out" \\
        --working-dir "\$(pwd)" \\
        --to-rds \\
        ${config_arg} \\
        ${transgene_arg} \\
        ${force_flag}

    test -d harmonize_out/rds
    ls harmonize_out/rds/*_harmonized.rds >/dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hkoca: \$(python -c "from hkoca._version import __version__; print(__version__)" 2>/dev/null || echo "unknown")
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
