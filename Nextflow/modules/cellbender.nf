process CELLBENDER {
    tag { sample_id }
    label 'gpu_optional'
    label 'process_high'

    conda "${params.conda_cellbender}"

    publishDir path: { "${params.outdir}/cellbender/${sample_id}" }, mode: params.publish_mode, enabled: params.publish_cellbender

    input:
    tuple val(sample_id), path(raw_input), val(input_mode)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_filtered.h5"), emit: filtered_h5
    path "${sample_id}_filtered.h5", emit: full_h5, optional: true
    path "versions.yml", emit: versions

    script:
    def config_arg = params.cellbender_config ? "--config ${params.cellbender_config}" : ""
    def cuda_flag  = params.cellbender_cuda ? "--cuda" : "--no-cuda"
    def h5_name    = params.cellbender_h5_filename
    def mtx_name   = params.cellbender_mtx_dirname
    """
    set -euo pipefail
    export HKOCA_ROOT="${projectDir}"
    export PYTHONPATH="${projectDir}\${PYTHONPATH:+:\$PYTHONPATH}"

    mkdir -p "${sample_id}"
    if [ "${input_mode}" = "h5" ]; then
        ln -sf "\$(readlink -f ${raw_input})" "${sample_id}/${h5_name}"
    else
        mkdir -p "${sample_id}/${mtx_name}"
        # Stage MTX directory contents (or the directory itself)
        if [ -d "${raw_input}" ]; then
            cp -a "${raw_input}/." "${sample_id}/${mtx_name}/"
        else
            echo "ERROR: MTX input must be a directory: ${raw_input}" >&2
            exit 1
        fi
    fi

    python -m hkoca.cellbender ${input_mode} \\
        --samples "\$(pwd)/${sample_id}" \\
        --working-dir "\$(pwd)" \\
        ${config_arg} \\
        ${cuda_flag} \\
        --epochs ${params.cellbender_epochs} \\
        --total-droplets-included ${params.cellbender_total_droplets} \\
        --learning-rate ${params.cellbender_learning_rate} \\
        --cpu-threads ${task.cpus} \\
        --force

    # Emit the cell-called matrix expected by harmonize_data_path
    if [ ! -f "${sample_id}/${sample_id}_filtered_filtered.h5" ]; then
        # Fallback: some CellBender builds only write the --output path
        if [ -f "${sample_id}/${sample_id}_filtered.h5" ]; then
            cp "${sample_id}/${sample_id}_filtered.h5" "${sample_id}_filtered_filtered.h5"
        else
            echo "ERROR: CellBender output missing for ${sample_id}" >&2
            exit 1
        fi
    else
        cp "${sample_id}/${sample_id}_filtered_filtered.h5" .
    fi
    if [ -f "${sample_id}/${sample_id}_filtered.h5" ]; then
        cp "${sample_id}/${sample_id}_filtered.h5" .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hkoca: \$(python -c "from hkoca._version import __version__; print(__version__)" 2>/dev/null || echo "unknown")
        cellbender: \$(cellbender --version 2>/dev/null | head -1 || echo "unknown")
    END_VERSIONS
    """
}
