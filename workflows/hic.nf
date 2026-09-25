/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-hic main analysis workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'

include { PREPARE_GENOME                    } from '../subworkflows/local/prepare_genome'
include { LIBRARY_HIC                       } from '../subworkflows/local/library_hic'
include { MERGE_COOLERS as MERGE_LIBRARY    } from '../subworkflows/local/merge_coolers'
include { MERGE_COOLERS as MERGE_REPLICATE  } from '../subworkflows/local/merge_coolers'

workflow HIC {
    main:
    if (!(params.aligner in ['bwa', 'bwa-mem2'])) {
        error("Invalid --aligner '${params.aligner}'. Use 'bwa' or 'bwa-mem2'.")
    }
    def mapq_filters = (params.min_mapq && params.min_mapq > 0) ? [0, params.min_mapq as int] : [0]

    //
    // Build one input item per samplesheet row and normalize the (optional) tech replicate
    //
    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    ch_rows = ch_samplesheet.map { meta, fastq_1, fastq_2 ->
        def tech    = meta.tech_replicate ?: 1
        def library = "${meta.id}_REP${meta.bio_replicate}_T${tech}"
        [ meta + [tech_replicate: tech, library: library, cool_id: library], [ fastq_1, fastq_2 ] ]
    }

    //
    // Assign a run index to every sequencing run, so split runs of the same library
    // get unique intermediate file names before they are merged
    //
    ch_reads = ch_rows
        .groupTuple(by: [0])
        .flatMap { meta, fastqs ->
            def runs = []
            fastqs.eachWithIndex { fastq, idx ->
                runs << [ meta + [run: idx + 1], fastq ]
            }
            runs
        }

    //
    // Reference preparation
    //
    PREPARE_GENOME(params.fasta, params.bwa_index, params.chrom_sizes)

    //
    // Per-library alignment, deduplication and contact maps
    //
    LIBRARY_HIC(
        ch_reads,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.chrom_sizes,
        params.bin_size,
        mapq_filters
    )

    //
    // Merged library: merge technical replicates within a biological replicate
    //
    ch_mlb = LIBRARY_HIC.out.cool_raw
        .map { meta, cool ->
            def suffix = meta.mapq > 0 ? ".Q${meta.mapq}" : ''
            [
                [
                    id             : meta.id,
                    bio_replicate  : meta.bio_replicate,
                    mapq           : meta.mapq,
                    cool_id        : "${meta.id}_REP${meta.bio_replicate}.mLb${suffix}"
                ],
                cool
            ]
        }
        .groupTuple(by: [0])
    MERGE_LIBRARY(ch_mlb)

    //
    // Merged replicate: merge biological replicates of a sample
    //
    ch_mrp = MERGE_LIBRARY.out.cool_raw
        .map { meta, cool ->
            def suffix = meta.mapq > 0 ? ".Q${meta.mapq}" : ''
            [
                [
                    id      : meta.id,
                    mapq    : meta.mapq,
                    cool_id : "${meta.id}.mRp${suffix}"
                ],
                cool
            ]
        }
        .groupTuple(by: [0])
    MERGE_REPLICATE(ch_mrp)
}
