/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Call CALDER2 compartments and nested sub-domains from a balanced multiresolution cooler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALDER2 } from '../../../modules/nf-core/calder2'

workflow CALDER2_HIC {
    take:
    ch_cool         // channel: [ val(meta), path(mcool) ] ; meta.cool_id is the output basename
    resolution      // integer: resolution group inside the .mcool

    main:
    CALDER2(ch_cool, resolution)

    emit:
    CALDER2.out.output_folder
}
