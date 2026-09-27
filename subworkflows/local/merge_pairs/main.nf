/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Merge pairs from several libraries into one file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PAIRTOOLS_MERGE } from '../../../modules/nf-core/pairtools/merge'

workflow MERGE_PAIRS {
    take:
    ch_pairs    // channel: [ val(meta), [ pairs, ... ] ] ; meta.cool_id is the output basename

    main:
    PAIRTOOLS_MERGE(ch_pairs)

    emit:
    PAIRTOOLS_MERGE.out.pairs
}
