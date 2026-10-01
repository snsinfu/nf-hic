/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Merge a group of contact maps and produce a multiresolution cooler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Merging is done on unbalanced raw counts; balancing is applied afterwards only.
*/

include { COOLER_MERGE   } from '../../../modules/nf-core/cooler/merge'
include { COOLER_BALANCE } from '../../../modules/local/cooler_balance'
include { COOLER_ZOOMIFY } from '../../../modules/nf-core/cooler/zoomify'

workflow MERGE_COOLERS {
    take:
    ch_cool        // channel: [ val(meta), [ cool, ... ] ] ; meta.cool_id is the output basename

    main:
    //
    // Pool raw counts
    //
    COOLER_MERGE(ch_cool)

    //
    // Balance a derived copy of the merged raw cooler (published) and zoomify the raw cooler
    //
    ch_cool_balanced = channel.empty()
    if (params.balance) {
        COOLER_BALANCE(COOLER_MERGE.out.cool.map { meta, cool -> [ meta, cool, '' ] })
        ch_cool_balanced = COOLER_BALANCE.out.cool
    }
    else {
        ch_cool_balanced = COOLER_MERGE.out.cool
    }
    COOLER_ZOOMIFY(COOLER_MERGE.out.cool)

    emit:
    cool     = ch_cool_balanced
    cool_raw = COOLER_MERGE.out.cool
    mcool    = COOLER_ZOOMIFY.out.mcool
}
