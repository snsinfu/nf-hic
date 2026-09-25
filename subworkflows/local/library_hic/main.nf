/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-library Hi-C processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align each sequencing run, parse/sort pairs, merge all runs of a library, remove
    duplicates, emit a deduped BAM + pairs, then build contact maps (one per MAPQ filter).
----------------------------------------------------------------------------------------
*/

include { BWA_MEM      } from '../../../modules/nf-core/bwa/mem'
include { BWAMEM2_MEM  } from '../../../modules/nf-core/bwamem2/mem'

include { PAIRTOOLS_PARSE  } from '../../../modules/nf-core/pairtools/parse'
include { PAIRTOOLS_SORT   } from '../../../modules/nf-core/pairtools/sort'
include { PAIRTOOLS_MERGE  } from '../../../modules/nf-core/pairtools/merge'
include { PAIRTOOLS_DEDUP  } from '../../../modules/nf-core/pairtools/dedup'
include { PAIRTOOLS_SPLIT  } from '../../../modules/nf-core/pairtools/split'
include { PAIRTOOLS_SELECT } from '../../../modules/nf-core/pairtools/select'

include { COOLER_CLOAD  } from '../../../modules/nf-core/cooler/cload'
include { COOLER_BALANCE } from '../../../modules/nf-core/cooler/balance'
include { COOLER_ZOOMIFY as COOLER_ZOOMIFY_LIBRARY } from '../../../modules/nf-core/cooler/zoomify'

workflow LIBRARY_HIC {
    take:
    ch_reads         // channel: [ val(meta), fastq_1, fastq_2 ]
    ch_bwa_index     // channel: [ val(meta), path(index) ]
    ch_fasta         // channel: path(fasta)
    ch_chrom_sizes   // channel: path(chrom_sizes)
    bin_size         // integer
    mapq_filters     // list of integers, e.g. [0] or [0, 30]

    main:
    //
    // Align each sequencing run (split runs share the same library id)
    //
    ch_bam = channel.empty()
    if (params.aligner == 'bwa') {
        BWA_MEM(ch_reads, ch_bwa_index, ch_fasta.map { fasta -> [ [:], fasta ] }, false)
        ch_bam = BWA_MEM.out.bam
    }
    else {
        BWAMEM2_MEM(ch_reads, ch_bwa_index, ch_fasta.map { fasta -> [ [:], fasta ] }, false)
        ch_bam = BWAMEM2_MEM.out.bam
    }

    //
    // Parse and sort each run, then merge all runs belonging to the same library
    // (drop the per-run index so runs group by library)
    //
    PAIRTOOLS_PARSE(ch_bam, ch_chrom_sizes)
    PAIRTOOLS_SORT(PAIRTOOLS_PARSE.out.pairsam)
    ch_sorted = PAIRTOOLS_SORT.out.sorted
        .map { meta, pairs -> [ meta.findAll { key, _value -> key != 'run' }, pairs ] }
    PAIRTOOLS_MERGE(ch_sorted.groupTuple(by: [0]))

    //
    // Remove duplicates once per library (split runs are deduplicated together)
    // and split the deduped pairsam into a BAM + pairs
    //
    PAIRTOOLS_DEDUP(PAIRTOOLS_MERGE.out.pairs)
    PAIRTOOLS_SPLIT(PAIRTOOLS_DEDUP.out.pairs)

    //
    // One contact-map set per MAPQ filter; Q0 keeps the unsuffixed base name
    //
    ch_pairs_filters = PAIRTOOLS_SPLIT.out.pairs
        .combine(channel.fromList(mapq_filters))
        .map { meta, pairs, mapq ->
            def suffix = mapq > 0 ? ".Q${mapq}" : ''
            [ meta + [mapq: mapq, cool_id: "${meta.library}${suffix}"], pairs ]
        }
    PAIRTOOLS_SELECT(ch_pairs_filters)
    COOLER_CLOAD(
        PAIRTOOLS_SELECT.out.selected.map { meta, pairs -> [ meta, pairs, [] ] },
        ch_chrom_sizes.map { sizes -> [ [:], sizes ] },
        'pairs',
        bin_size
    )

    //
    // Balance a derived copy of the raw cooler (published); zoomify the raw cooler
    //
    ch_cool = channel.empty()
    if (params.balance) {
        COOLER_BALANCE(COOLER_CLOAD.out.cool.map { meta, cool -> [ meta, cool, '' ] })
        ch_cool = COOLER_BALANCE.out.cool
    }
    else {
        ch_cool = COOLER_CLOAD.out.cool
    }
    COOLER_ZOOMIFY_LIBRARY(COOLER_CLOAD.out.cool)

    emit:
    bam      = PAIRTOOLS_SPLIT.out.bam
    pairs    = PAIRTOOLS_SELECT.out.selected
    cool     = ch_cool
    cool_raw = COOLER_CLOAD.out.cool
    mcool    = COOLER_ZOOMIFY_LIBRARY.out.mcool
}
