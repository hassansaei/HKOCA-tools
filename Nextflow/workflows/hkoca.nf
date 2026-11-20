include { CELLBENDER              } from '../modules/cellbender.nf'
include { PREPARE_HARMONIZE_META  } from '../modules/prepare_meta.nf'
include { HARMONIZATION           } from '../modules/harmonization.nf'
include { QC                      } from '../modules/qc.nf'
include { ANNOTATION              } from '../modules/annotation.nf'
include { INTEGRATION             } from '../modules/integration.nf'
include { PROJECTION              } from '../modules/projection.nf'

def parseBool(val, defaultVal) {
    if (val == null) return defaultVal
    def s = val.toString().trim()
    if (s == '') return defaultVal
    return s.toLowerCase() in ['1', 'true', 't', 'yes', 'y']
}

workflow HKOCA {

    take:
    sample_info_ch
    gtf_ch
    atlas_ch
    model_dir_ch

    main:

    samples_ch = sample_info_ch
        .splitCsv(header: true)
        .map { row ->
            def sid = row.sample_id.toString().trim()
            def run_cb = parseBool(row.run_cellbender, true)
            def data_path = row.data_path ? row.data_path.toString().trim() : ''
            def sample_dir = row.sample_dir ? row.sample_dir.toString().trim() : ''
            if (!sample_dir) {
                if (data_path.toLowerCase().endsWith('.h5') ||
                    data_path.toLowerCase().endsWith('.h5ad') ||
                    data_path.toLowerCase().endsWith('.hdf5')) {
                    sample_dir = new File(data_path).parent
                } else {
                    sample_dir = data_path
                }
            }
            def mode = params.cellbender_mode ?: 'h5'
            def raw
            if (mode == 'h5') {
                def candidate = (data_path && data_path.endsWith('.h5') && !data_path.contains('_filtered'))
                    ? data_path
                    : "${sample_dir}/${params.cellbender_h5_filename}"
                raw = file(candidate, checkIfExists: false)
            } else {
                raw = file("${sample_dir}/${params.cellbender_mtx_dirname}", checkIfExists: false)
            }
            tuple(sid, raw, mode, run_cb)
        }

    no_file_ch = Channel.of('').collectFile(name: 'NO_FILE')

    if (!params.skip_cellbender) {
        CELLBENDER(
            samples_ch
                .filter { _sid, _raw, _mode, run_cb -> run_cb }
                .map { sid, raw, mode, _run -> tuple(sid, raw, mode) }
        )
        cb_map_ch = CELLBENDER.out.filtered_h5
            .map { sid, h5 -> "${sid}\t${h5.toString()}" }
            .collectFile(name: 'cellbender_map.tsv', newLine: true, sort: true)
    } else {
        cb_map_ch = Channel.of('').collectFile(name: 'empty_cellbender_map.tsv')
    }

    PREPARE_HARMONIZE_META(sample_info_ch, cb_map_ch)
    HARMONIZATION(PREPARE_HARMONIZE_META.out.csv, gtf_ch)
    QC(HARMONIZATION.out.output_root)
    ANNOTATION(QC.out.filtered_h5ad)

    rds_keyed = QC.out.filtered_rds.map { rds -> tuple(rds.simpleName, rds) }
    ann_keyed = ANNOTATION.out.annotated.map { _in_h5ad, ann ->
        tuple(ann.simpleName.replaceFirst(/_annotated$/, ''), ann)
    }

    integration_in = rds_keyed
        .join(ann_keyed, by: 0, remainder: true)
        .combine(no_file_ch)
        .map { _key, rds, ann, stub ->
            tuple(rds, ann != null ? ann : stub)
        }

    INTEGRATION(integration_in)

    if (params.run_projection) {
        prep_ch = INTEGRATION.out.prepared.map { _filtered, prep -> prep }
        PROJECTION(prep_ch, atlas_ch.first(), model_dir_ch.first())
        projection_out = PROJECTION.out.projected_h5ad
    } else {
        projection_out = Channel.empty()
    }

    emit:
    harmonize_csv    = PREPARE_HARMONIZE_META.out.csv
    harmonize_rds    = HARMONIZATION.out.rds
    qc_filtered_rds  = QC.out.filtered_rds
    qc_filtered_h5ad = QC.out.filtered_h5ad
    annotated        = ANNOTATION.out.annotated
    integrated_prep  = INTEGRATION.out.prepared
    projected        = projection_out
}
