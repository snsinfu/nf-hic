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
    take:
    fasta          // path: resolved reference fasta (from --fasta or --genome)
    index          // path: resolved aligner index directory (from --bwa_index/--bwamem2_index or --genome); may be null
    chrom_sizes    // path: resolved chrom.sizes (from --chrom_sizes or --genome); may be null
    aligner        // string: selected aligner, e.g. 'bwa' or 'bwa-mem2'

    main:
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def keys = params.genomes.keySet().join(', ')
        error("Genome '${params.genome}' not found in any config file provided to the pipeline. Available genome keys: ${keys}")
    }
    if (!fasta) {
        error("Genome fasta file not specified: use --genome <key>, --fasta <file.fa>, or a custom config.")
    }
    if (!(aligner in ['bwa', 'bwa-mem2'])) {
        error("Invalid --aligner '${aligner}'. Use 'bwa' or 'bwa-mem2'.")
    }
    //
    // MAPQ filters: always keep the unfiltered (Q0) set, then add one per requested
    // threshold. --min_mapq accepts a comma-separated list, e.g. --min_mapq 0,30,60.
    //
    def mapq_thresholds = []
    if (params.min_mapq != null) {
        def raw   = params.min_mapq
        def items = (raw instanceof List) ? raw : raw.toString().split(',')
        mapq_thresholds = items
            .collect { item -> item.toString().trim() }
            .findAll { item -> item }
            .collect { token ->
                if (!(token ==~ /\d+/)) {
                    error("Invalid --min_mapq value '${token}'. Use a comma-separated list of non-negative integers, e.g. --min_mapq 0,30,60")
                }
                token as int
            }
    }
    def mapq_filters = ([0] + mapq_thresholds).unique().sort()

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
    PREPARE_GENOME(fasta, index, chrom_sizes, aligner)

    //
    // Per-library alignment, deduplication and contact maps
    //
    LIBRARY_HIC(
        ch_reads,
        PREPARE_GENOME.out.index,
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
        // One bio replicate has nothing to pool: skip it, its .mLb is the sample map
        .filter { _meta, cools -> cools.size() > 1 }
    MERGE_REPLICATE(ch_mrp)
}
