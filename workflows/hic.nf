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
include { MERGE_STATS as MERGE_LIBRARY_STATS   } from '../subworkflows/local/merge_stats'
include { MERGE_STATS as MERGE_REPLICATE_STATS } from '../subworkflows/local/merge_stats'
include { MERGE_PAIRS as MERGE_LIBRARY_PAIRS   } from '../subworkflows/local/merge_pairs'
include { MERGE_PAIRS as MERGE_REPLICATE_PAIRS } from '../subworkflows/local/merge_pairs'
include { CALDER2_HIC as MERGE_LIBRARY_CALDER2   } from '../subworkflows/local/calder2'
include { CALDER2_HIC as MERGE_REPLICATE_CALDER2 } from '../subworkflows/local/calder2'
include { CALDER2_GENE_DENSITY } from '../modules/local/calder2_gene_density'
include { ONTAD_HIC as MERGE_LIBRARY_ONTAD   } from '../subworkflows/local/ontad'
include { ONTAD_HIC as MERGE_REPLICATE_ONTAD } from '../subworkflows/local/ontad'

workflow HIC {
    take:
    fasta          // path: resolved reference fasta (from --fasta or --genome)
    index          // path: resolved aligner index directory (from --bwa_index/--bwamem2_index/--bwamem3_index or --genome); may be null
    chrom_sizes    // path: resolved chrom.sizes (from --chrom_sizes or --genome); may be null
    gtf            // path: resolved GTF annotation (from --gtf or --genome); may be null
    aligner        // string: selected aligner, e.g. 'bwa', 'bwa-mem2' or 'bwa-mem3'

    main:
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def keys = params.genomes.keySet().join(', ')
        error("Genome '${params.genome}' not found in any config file provided to the pipeline. Available genome keys: ${keys}")
    }
    if (!fasta) {
        error("Genome fasta file not specified: use --genome <key>, --fasta <file.fa>, or a custom config.")
    }
    if (!(aligner in ['bwa', 'bwa-mem2', 'bwa-mem3'])) {
        error("Invalid --aligner '${aligner}'. Use 'bwa', 'bwa-mem2' or 'bwa-mem3'.")
    }

    //
    // CALDER2 A/B phasing reference. Precedence:
    //   1. explicit --calder2_feature_track
    //   2. auto-generated gene density (non-built-in genome; requires a GTF)
    //   3. none -> CALDER's built-in reference compartments (hg19|hg38|mm9|mm10)
    //
    def calder2_genome  = params.calder2_genome ?: params.genome
    def calder2_builtin = ['hg19', 'hg38', 'mm9', 'mm10'].contains(calder2_genome)
    def calder2_enabled = !params.skip_calder2 && params.balance
    def autogen_track   = calder2_enabled && !calder2_builtin && !params.calder2_feature_track
    if (autogen_track && !gtf) {
        error("CALDER2 A/B phasing for genome '${calder2_genome}' needs a gene-density track, but no GTF annotation was found. " +
              "Use --gtf <genes.gtf> or a --genome with a 'gtf' catalog entry.")
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
        def library = "${meta.id}_REP${meta.replicate}_T${tech}"
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
                    replicate      : meta.replicate,
                    mapq           : meta.mapq,
                    cool_id        : "${meta.id}_REP${meta.replicate}.mLb${suffix}"
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

    //
    // Merged stats: merge pairtools stats per biological replicate, then per sample
    //
    if (!params.skip_pairtools_stats) {
        ch_mlb_stats = LIBRARY_HIC.out.stats
            .map { meta, stat ->
                [
                    [
                        id            : meta.id,
                        replicate     : meta.replicate,
                        cool_id       : "${meta.id}_REP${meta.replicate}.mLb"
                    ],
                    stat
                ]
            }
            .groupTuple(by: [0])
        MERGE_LIBRARY_STATS(ch_mlb_stats)

        ch_mrp_stats = MERGE_LIBRARY_STATS.out
            .map { meta, stat -> [ [ id: meta.id, cool_id: "${meta.id}.mRp" ], stat ] }
            .groupTuple(by: [0])
            .filter { _meta, stats -> stats.size() > 1 }
        MERGE_REPLICATE_STATS(ch_mrp_stats)
    }

    //
    // Merged pairs: merge Q0 pairs per biological replicate, then per sample.
    // Opt-in because pairtools merge is much more expensive than cooler merge.
    //
    if (params.merge_pairs) {
        ch_mlb_pairs = LIBRARY_HIC.out.pairs
            .filter { meta, _pairs -> meta.mapq == 0 }
            .map { meta, pairs ->
                [
                    [
                        id            : meta.id,
                        replicate     : meta.replicate,
                        cool_id       : "${meta.id}_REP${meta.replicate}.mLb"
                    ],
                    pairs
                ]
            }
            .groupTuple(by: [0])
        MERGE_LIBRARY_PAIRS(ch_mlb_pairs)

        ch_mrp_pairs = MERGE_LIBRARY_PAIRS.out
            .map { meta, pairs -> [ [ id: meta.id, cool_id: "${meta.id}.mRp" ], pairs ] }
            .groupTuple(by: [0])
            .filter { _meta, pairs -> pairs.size() > 1 }
        MERGE_REPLICATE_PAIRS(ch_mrp_pairs)
    }

    //
    // CALDER2 compartments + nested sub-domains (opt-out via --skip_calder2).
    // Needs a balanced cooler: the CALDER2 CLI dumps pixels with `--balanced`.
    //
    if (!params.skip_calder2 && !params.balance) {
        log.warn "[nf-hic] CALDER2 requires a balanced cooler; skipping because --balance false."
    }
    if (calder2_enabled) {
        //
        // Resolve the A/B phasing track: explicit --calder2_feature_track > auto gene density
        // (non-built-in genome) > [] (CALDER's built-in reference for hg19|hg38|mm9|mm10).
        // The track is a declared module input, so Nextflow stages it into the task.
        //
        def ch_feature_track = channel.empty()
        if (params.calder2_feature_track) {
            if (calder2_builtin) {
                log.warn "[nf-hic] --calder2_feature_track overrides CALDER's built-in reference for '${calder2_genome}' " +
                         "(CALDER treats the genome as 'others'; the built-in reference is not used)."
            }
            ch_feature_track = channel.value(file(params.calder2_feature_track, checkIfExists: true))
        }
        else if (autogen_track) {
            CALDER2_GENE_DENSITY(
                PREPARE_GENOME.out.chrom_sizes,
                channel.value(file(gtf, checkIfExists: true)),
                channel.value(params.calder2_resolution)
            )
            ch_feature_track = CALDER2_GENE_DENSITY.out.bed.first()
        }
        else {
            ch_feature_track = channel.value([])
        }

        MERGE_LIBRARY_CALDER2(MERGE_LIBRARY.out.mcool, channel.value(params.calder2_resolution), ch_feature_track)
        MERGE_REPLICATE_CALDER2(MERGE_REPLICATE.out.mcool, channel.value(params.calder2_resolution), ch_feature_track)
    }

    //
    // OnTAD hierarchical TADs (opt-out via --skip_ontad). Needs a balanced
    // cooler: the matrix dump reads the stored balancing weights by default.
    //
    if (!params.skip_ontad && !params.balance) {
        log.warn "[nf-hic] OnTAD requires a balanced cooler; skipping because --balance false."
    }
    if (!params.skip_ontad && params.balance) {
        def ontad_metadata = [
            minsz  : params.ontad_minsz,
            maxsz  : params.ontad_maxsz,
            lsize  : params.ontad_lsize,
            ldiff  : params.ontad_ldiff,
            penalty: params.ontad_penalty,
        ]
        MERGE_LIBRARY_ONTAD(MERGE_LIBRARY.out.mcool, PREPARE_GENOME.out.chrom_sizes, params.ontad_resolution, ontad_metadata)
        MERGE_REPLICATE_ONTAD(MERGE_REPLICATE.out.mcool, PREPARE_GENOME.out.chrom_sizes, params.ontad_resolution, ontad_metadata)
    }
}
