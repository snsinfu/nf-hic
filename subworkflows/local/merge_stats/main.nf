/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Merge pairtools stats from several libraries into one file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PAIRTOOLS_STATS } from '../../../modules/nf-core/pairtools/stats'

workflow MERGE_STATS {
    take:
    ch_stats    // channel: [ val(meta), [ stat, ... ] ] ; meta.cool_id is the output basename

    main:
    PAIRTOOLS_STATS(ch_stats)

    emit:
    PAIRTOOLS_STATS.out.stats
}
