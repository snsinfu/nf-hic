/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Call CALDER2 compartments and nested sub-domains from a balanced multiresolution cooler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALDER2 } from '../../../modules/local/calder2'

workflow CALDER2_HIC {
    take:
    ch_cool            // channel: [ val(meta), path(mcool) ] ; meta.cool_id is the output basename
    resolution         // integer: resolution group inside the .mcool
    ch_feature_track   // channel: path(bed) or [] (empty = use CALDER's built-in reference)

    main:
    CALDER2(ch_cool, resolution, ch_feature_track)

    emit:
    CALDER2.out.output_folder
}
