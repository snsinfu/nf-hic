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
    save_intermediates // boolean: keep and publish the full intermediate_data/ folder

    main:
    CALDER2(ch_cool, resolution, ch_feature_track, save_intermediates)

    emit:
    CALDER2.out.output_folder
}
