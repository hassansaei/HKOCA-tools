process ANNOTATION {
    tag { h5ad.simpleName }
    label 'process_medium'

    conda "${params.conda_harmonize}"

    publishDir "${params.outdir}/annotation", mode: params.publish_mode

    input:
    path h5ad

    output:
    tuple path(h5ad), path("annotation_out/annotated_obj/*_annotated.h5ad"), emit: annotated
    path "annotation_out", emit: output_dir
    path "versions.yml", emit: versions

    script:
    def markers = params.annotation_markers ?: "${projectDir}/hkoca/config/snapseed_markers_v4.yaml"
    def config_arg = params.annotation_config ? "--config ${params.annotation_config}" : ""
    def res_arg    = params.annotation_resolutions ? "--resolutions ${params.annotation_resolutions}" : ""
    def force_flag = params.force ? "--force" : ""
    def plots_flag = params.annotation_save_plots ? "--save-plots" : "--no-save-plots"
    """
    set -euo pipefail
    export HKOCA_ROOT="${projectDir}"
    export PYTHONPATH="${projectDir}\${PYTHONPATH:+:\$PYTHONPATH}"

    mkdir -p annotation_out

    python -m hkoca.annotation run \\
        --input "${h5ad}" \\
        --output-dir "\$(pwd)/annotation_out" \\
        --markers "${markers}" \\
        ${config_arg} \\
        ${res_arg} \\
        ${plots_flag} \\
        ${force_flag}

    ls annotation_out/annotated_obj/*_annotated.h5ad >/dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hkoca: \$(python -c "from hkoca._version import __version__; print(__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
