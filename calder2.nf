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
 * --calder2_mcool can also be a map from sample id to path/glob (or a list of
 * paths/globs). This standalone pipeline does not infer sample ids from
 * filenames: a bare glob/string has an empty sample id and uses only the
 * shared/default feature track; use the map form to assign per-sample tracks:
 *
 *   calder2_mcool:
 *     sampleA: 'results/bwa/merged_replicate/sampleA*.mcool'
 *     sampleB: ['results/bwa/merged_library/sampleB_REP1.mLb.mcool',
 *               'results/bwa/merged_library/sampleB_REP2.mLb.mcool']
 *
 * The A/B phasing track is shared by default but can be overridden per sample.
 * With --calder2_feature_track, keys are the sample ids above. A -params-file
 * map with '*' as the shared fallback works:
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
    // --calder2_mcool accepts:
    //   - a path/glob or comma-separated string (empty sample id: default track only);
    //   - a list of paths/globs (same);
    //   - a map of sample id -> path/glob, or sample id -> [path/glob, ...].
    // A bare glob/string has an empty sample id and uses only the shared/default
    // feature track. Each file becomes one output entry named by its base name;
    // materialize the list to fail early on empty globs or colliding output names,
    // then fan the entries back out.
    //
    def mcool_patterns = []   // [sample: <id|''>, pattern: <glob>]
    if (params.calder2_mcool instanceof Map) {
        params.calder2_mcool.each { sample, spec ->
            def sampleId = sample == null ? null : sample.toString().trim()
            if (!sampleId) {
                error("Invalid --calder2_mcool map: empty sample id.")
            }
            def patterns = (spec instanceof List)
                ? spec.collect { item -> item.toString().trim() }
                : [ spec.toString().trim() ]
            patterns.each { pattern ->
                if (!pattern) {
                    error("Invalid --calder2_mcool map for '${sampleId}': empty path/glob.")
                }
                mcool_patterns << [sample: sampleId, pattern: pattern]
            }
        }
    }
    else {
        def patterns = (params.calder2_mcool instanceof List)
            ? params.calder2_mcool.collect { item -> item.toString() }
            : params.calder2_mcool.toString().split(',').collect { token -> token.trim() }
        patterns.findAll { pattern -> pattern }.each { pattern ->
            mcool_patterns << [sample: '', pattern: pattern]
        }
    }
    if (mcool_patterns.isEmpty()) {
        error("No mcool input specified: use --calder2_mcool <file|glob|comma-separated list|sample->glob map>.")
    }

    def mcool_files = []
    mcool_patterns.each { entry ->
        def matches = files(entry.pattern)
        if (!matches) {
            error("--calder2_mcool pattern '${entry.pattern}' matched no files.")
        }
        matches.each { path ->
            mcool_files << [sample: entry.sample, path: path]
        }
    }

    ch_cool = channel.fromList(mcool_files)
        .map { entry ->
            def path = entry.path
            if (path.isDirectory()) {
                error("--calder2_mcool expects a cooler file, got directory: ${path}")
            }
            def name = path.baseName
            [ [id: name, cool_id: name, sample: entry.sample], path ]
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
    // CALDER2_ATTACH_TRACK applies this precedence per input (key = --calder2_mcool map key;
    // a bare input has an empty sample id and uses only the shared/default track):
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

    ch_known_sample_ids = ch_cool.map { meta, _cool -> meta.sample }.toList()

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
