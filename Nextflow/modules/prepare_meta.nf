process PREPARE_HARMONIZE_META {
    tag { "meta" }
    label 'process_single'

    publishDir "${params.outdir}/.hkoca", mode: params.publish_mode

    input:
    path sample_info
    path cellbender_map

    output:
    path "harmonize_metadata.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    def skip_flag = params.skip_cellbender ? "--skip-cellbender" : ""
    """
    set -euo pipefail

    python "${projectDir}/modules/prepare_harmonize_csv.py" \\
        --sample-info "${sample_info}" \\
        --output harmonize_metadata.csv \\
        --cellbender-map "${cellbender_map}" \\
        --working-dir "${params.working_dir}" \\
        --output-suffix "${params.cellbender_output_suffix}" \\
        ${skip_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_harmonize_csv: "1.0"
    END_VERSIONS
    """
}
