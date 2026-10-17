#!/usr/bin/env nextflow
/*
 * Standalone CALDER2 sub-pipeline: run CALDER2 compartments and nested
 * sub-domains on existing balanced multiresolution coolers, without repeating
 * alignment or pairtools/cooler construction.
 *
 * This is a third entry point of the nf-hic repository; the full pipeline is
 * main.nf and the standalone OnTAD sub-pipeline is ontad.nf. All share
 * modules/, subworkflows/, workflows/, conf/ and bin/.
 *
 * Usage, remote (published repository):
 *
 *   nextflow run snsinfu/nf-hic \
 *       -main-script calder2.nf \
 *       --calder2_mcool 'results/bwa/merged_library/*.mcool' \
 *       --calder2_resolution 100000 \
 *       --calder2_genome hg38
 *
 * Usage, local checkout (from the repository root):
 *
 *   nextflow run calder2.nf \
 *       --calder2_mcool 'results/bwa/merged_library/*.mcool' \
 *       --calder2_resolution 100000 \
 *       --calder2_genome hg38
 *
 * Nextflow loads the root nextflow.config automatically (including
 * conf/modules.config and conf/base.config), so no -C/-c config file is needed.
 * A bare `snsinfu/nf-hic/calder2.nf` project coordinate is rejected by Nextflow
 * ("Repository URL must not end with a script file extension"); use -main-script
 * for remote runs.
 *
 * One or more mcool paths/globs may be given. --calder2_mcool takes a single
 * path/glob or a comma-separated list (repeated flags do not accumulate in
 * Nextflow: the last one wins). A YAML list also works via -params-file. Each
 * input must be a balanced multiresolution cooler (a plain .cool cannot select
 * the resolution group), and --calder2_resolution must already exist under the
 * file's /resolutions group. Outputs are published per input base name to
 * <outdir>/calder2/<cool_id>/{sub_compartments,sub_domains}/.
 *
 * Values that start with '-' (e.g. --calder2_args '--adaptive') need the
 * --option=value form, otherwise Nextflow's CLI reads them as options.
 *
 * projectDir is the repo root here, so the shared bin/ helpers are used
 * directly (no per-subpipeline symlink).
 */

include { CALDER2_HIC } from './subworkflows/local/calder2'
include { CALDER2_GENE_DENSITY } from './modules/local/calder2_gene_density'

workflow {
    if (!params.calder2_mcool) {
        error("No mcool input specified: use --calder2_mcool <file|glob|comma-separated list>.")
    }
    if (!(params.calder2_resolution as int > 0)) {
        error("Invalid --calder2_resolution '${params.calder2_resolution}': expected a positive integer.")
    }

    //
    // --calder2_mcool is a path/glob, a comma-separated string, or a list (from
    // -params-file). Each input becomes one entry whose base name is the output
    // prefix. Materialize the list to fail early on empty globs or colliding
    // output names, then fan the entries back out as a queue channel.
    //
    def mcool_paths = (params.calder2_mcool instanceof List)
        ? params.calder2_mcool.collect { item -> item.toString() }
        : params.calder2_mcool.toString().split(',').collect { token -> token.trim() }.findAll { token -> token }

    ch_cool = channel.fromPath(mcool_paths, checkIfExists: true)
        .map { path ->
            if (path.isDirectory()) {
                error("--calder2_mcool expects a cooler file, got directory: ${path}")
            }
            def name = path.baseName
            [ [id: name, cool_id: name], path ]
        }
        .toList()
        .map { entries ->
            def duplicate_ids = entries
                .groupBy { entry -> entry[0].cool_id }
                .findAll { _name, group -> group.size() > 1 }
                .keySet()
            if (duplicate_ids) {
                error("Duplicate --calder2_mcool base name(s): ${duplicate_ids.join(', ')}. " +
                      "Each output is named <base>, so rename the inputs or process them separately.")
            }
            entries
        }
        .flatMap { entries -> entries }

    //
    // CALDER2 A/B phasing track. Precedence:
    //   1. explicit --calder2_feature_track
    //   2. auto-generated gene density (non-built-in genome; needs --gtf and --chrom_sizes)
    //   3. none -> CALDER's built-in reference compartments (hg19|hg38|mm9|mm10)
    //
    def calder2_genome  = params.calder2_genome ?: params.genome
    def calder2_builtin = ['hg19', 'hg38', 'mm9', 'mm10'].contains(calder2_genome)

    def ch_feature_track = channel.empty()
    if (params.calder2_feature_track) {
        if (calder2_builtin) {
            log.warn "[nf-hic] --calder2_feature_track overrides CALDER's built-in reference for '${calder2_genome}' " +
                     "(CALDER treats the genome as 'others'; the built-in reference is not used)."
        }
        ch_feature_track = channel.value(file(params.calder2_feature_track, checkIfExists: true))
    }
    else if (!calder2_builtin) {
        def chrom_sizes = params.chrom_sizes ?: getGenomeAttribute('chrom_sizes')
        def gtf         = params.gtf         ?: getGenomeAttribute('gtf')
        if (!gtf) {
            error("CALDER2 A/B phasing for genome '${calder2_genome}' needs a gene-density track, but no GTF annotation was found. " +
                  "Use --gtf <genes.gtf> or a --genome with a 'gtf' catalog entry.")
        }
        if (!chrom_sizes) {
            error("CALDER2 A/B phasing for genome '${calder2_genome}' needs --chrom_sizes to generate a gene-density track. " +
                  "Use --chrom_sizes <chrom.sizes> or a --genome with a 'chrom_sizes' catalog entry.")
        }
        CALDER2_GENE_DENSITY(
            channel.value(file(chrom_sizes, checkIfExists: true)),
            channel.value(file(gtf, checkIfExists: true)),
            channel.value(params.calder2_resolution)
        )
        // All inputs are value channels, so the output is a single-value channel;
        // no `.first()` needed (it would log "operator `first` is useless").
        ch_feature_track = CALDER2_GENE_DENSITY.out.bed
    }
    else {
        ch_feature_track = channel.value([])
    }

    CALDER2_HIC(ch_cool, channel.value(params.calder2_resolution), ch_feature_track, channel.value(params.save_calder2_intermediates))
}

/*
 * Get an attribute from the genome catalog, e.g. gtf/chrom_sizes (same helper
 * as main.nf). Resolve references from the catalog only when not given
 * explicitly on the CLI, so --gtf/--chrom_sizes always win.
 */
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}
