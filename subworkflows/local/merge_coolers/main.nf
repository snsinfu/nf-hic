/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Merge a group of contact maps and produce a multiresolution cooler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COOLER_MERGE   } from '../../../modules/nf-core/cooler/merge'
include { COOLER_ZOOMIFY } from '../../../modules/nf-core/cooler/zoomify'

workflow MERGE_COOLERS {
    take:
    ch_cool        // channel: [ val(meta), [ cool, ... ] ] ; meta.cool_id is the output basename

    main:
    COOLER_MERGE(ch_cool)
    COOLER_ZOOMIFY(COOLER_MERGE.out.cool)

    emit:
    cool  = COOLER_MERGE.out.cool
    mcool = COOLER_ZOOMIFY.out.mcool
}
