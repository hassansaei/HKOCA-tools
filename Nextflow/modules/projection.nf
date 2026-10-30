process PROJECTION {
    tag { prepared_rds.simpleName }
    label 'gpu'
    label 'process_high'

    conda "${params.conda_projection}"

    publishDir path: { "${params.outdir}/projection/${prepared_rds.simpleName}" }, mode: params.publish_mode

    input:
    path prepared_rds
    path atlas_h5ad
    path model_dir

    output:
    path "projection_out", emit: output_dir
    path "projection_out/projected_obj/*_projected.h5ad", emit: projected_h5ad
    path "versions.yml", emit: versions

    script:
    def config_arg = params.projection_config ? "--config ${params.projection_config}" : ""
    def force_flag = params.force ? "--force" : ""
    def joint_flag = params.projection_joint_umap ? "--joint-umap" : ""
    """
    set -euo pipefail
    export HKOCA_ROOT="${projectDir}"
    export PYTHONPATH="${projectDir}\${PYTHONPATH:+:\$PYTHONPATH}"
    export HKOCA_PROJECTION_ENV="${params.projection_env_name}"
    export HKOCA_PROJECTION_INNER=1

    mkdir -p projection_out

    python -m hkoca.projection map \\
        --query "${prepared_rds}" \\
        --atlas "${atlas_h5ad}" \\
        --model-dir "${model_dir}" \\
        --output-dir "\$(pwd)/projection_out" \\
        ${config_arg} \\
        ${joint_flag} \\
        ${force_flag}

    ls projection_out/projected_obj/*_projected.h5ad >/dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hkoca: \$(python -c "from hkoca._version import __version__; print(__version__)" 2>/dev/null || echo "unknown")
        torch: \$(python -c "import torch; print(torch.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
