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
 * The A/B phasing track is shared by default but can be overridden per input.
 * With --calder2_feature_track, keys are the .mcool base names:
 *
 *   --calder2_feature_track 'shared.bed,sampleA=A.bed,sampleB=B.bed'
 *
 * or, from a -params-file, a map with '*' as the shared fallback:
 *
 *   calder2_feature_track:
 *     '*': shared.bed
 *     sampleA: A.bed
 *
 * Precedence: per-sample > shared > auto gene density > CALDER's built-in reference.
 * See CALDER2_FEATURE_TRACK for the exact syntax.
 *
 * Values that start with '-' (e.g. --calder2_args '--adaptive') need the
 * --option=value form, otherwise Nextflow's CLI reads them as options.
 *
 * projectDir is the repo root here, so the shared bin/ helpers are used
 * directly (no per-subpipeline symlink).
 */

include { CALDER2_HIC } from './subworkflows/local/calder2'
include { CALDER2_FEATURE_TRACK; CALDER2_ATTACH_TRACK } from './subworkflows/local/calder2_feature_track'

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
    // CALDER2 A/B phasing track. CALDER2_FEATURE_TRACK resolves the default and
    // CALDER2_ATTACH_TRACK applies this precedence per input (key = mcool base name):
    //   1. per-sample --calder2_feature_track entry
    //   2. shared --calder2_feature_track entry
    //   3. auto-generated gene density (non-built-in genome; needs --gtf and --chrom_sizes)
    //   4. none -> CALDER's built-in reference compartments (hg19|hg38|mm9|mm10)
    //
    def calder2_genome  = params.calder2_genome ?: params.genome
    def calder2_builtin = ['hg19', 'hg38', 'mm9', 'mm10'].contains(calder2_genome)
    def chrom_sizes     = params.chrom_sizes ?: getGenomeAttribute('chrom_sizes')
    def gtf             = params.gtf         ?: getGenomeAttribute('gtf')

    CALDER2_FEATURE_TRACK(
        params.calder2_feature_track,
        calder2_genome,
        calder2_builtin,
        gtf,
        chrom_sizes ? channel.value(file(chrom_sizes, checkIfExists: true)) : null,
        params.calder2_resolution
    )

    ch_known_sample_ids = ch_cool.map { meta, _cool -> meta.id }.toList()

    ch_cool_calder = CALDER2_ATTACH_TRACK(
        ch_cool,
        CALDER2_FEATURE_TRACK.out.track,
        CALDER2_FEATURE_TRACK.out.per_sample,
        ch_known_sample_ids,
        true
    )

    CALDER2_HIC(ch_cool_calder, channel.value(params.calder2_resolution), channel.value(params.save_calder2_intermediates), channel.value(params.calder2_nproc as int))
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
