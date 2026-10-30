#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HKOCA } from './workflows/hkoca.nf'

def helpMessage() {
    log.info """
    HKOCA-tools Nextflow

    Required:
      --sample_info   Sample info CSV
      --gtf           GRCh38 GTF
      --outdir        Results directory

    Optional:
      --skip_cellbender
      --run_projection --atlas <h5ad> --model_dir <scPoli model dir>

    See docs/nextflow.md
    """.stripIndent()
}

workflow {
    if (params.help) {
        helpMessage()
        return
    }

    if (!params.sample_info) error "Missing --sample_info"
    if (!params.gtf)         error "Missing --gtf"
    if (!params.outdir)      error "Missing --outdir"

    if (params.run_projection) {
        if (!params.atlas)     error "--run_projection requires --atlas"
        if (!params.model_dir) error "--run_projection requires --model_dir"
    }

    sample_info_ch = Channel.fromPath(params.sample_info, checkIfExists: true)
    gtf_ch         = Channel.fromPath(params.gtf, checkIfExists: true)

    if (params.run_projection) {
        atlas_ch     = Channel.fromPath(params.atlas, checkIfExists: true)
        model_dir_ch = Channel.fromPath(params.model_dir, checkIfExists: true, type: 'dir')
    } else {
        atlas_ch     = Channel.empty()
        model_dir_ch = Channel.empty()
    }

    HKOCA(sample_info_ch, gtf_ch, atlas_ch, model_dir_ch)
}
