#!/usr/bin/env nextflow
/*
 * Re-run the OnTAD subworkflow (dump -> OnTAD -> assemble) on existing
 * balanced multiresolution coolers, without repeating alignment or
 * pairtools/cooler construction.
 *
 * Usage, remote (published repository):
 *
 *   nextflow run snsinfu/nf-hic \
 *       -main-script subpipeline/ontad.nf \
 *       -C subpipeline/ontad.config \
 *       --mcool 'results/bwa/merged_library/*.mcool' \
 *       --chrom_sizes _data/hg38/chrom.sizes \
 *       --ontad_resolution 10000
 *
 * Usage, local checkout (from the repository root):
 *
 *   nextflow run subpipeline/ontad.nf \
 *       -C subpipeline/ontad.config \
 *       --mcool 'results/bwa/merged_library/*.mcool' \
 *       --chrom_sizes _data/hg38/chrom.sizes \
 *       --ontad_resolution 10000
 *
 * A bare `snsinfu/nf-hic/subpipeline/ontad.nf` project coordinate is rejected
 * by Nextflow ("Repository URL must not end with a script file extension");
 * use -main-script for remote runs. The full pipeline stays `nextflow run
 * snsinfu/nf-hic` (main.nf).
 *
 * -C is used on purpose: it loads only subpipeline/ontad.config (plus
 * conf/base.config) instead of the main pipeline nextflow.config, so this is
 * a standalone subworkflow run. -c also works but additionally loads the main
 * config.
 *
 * Values that start with '-' (e.g. --ontad_args '-log2') need the
 * --option=value form, otherwise Nextflow's CLI reads them as options.
 *
 * One or more mcool paths/globs may be given. --mcool takes a single
 * path/glob or a comma-separated list (repeated --mcool flags do not
 * accumulate in Nextflow: the last one wins). A YAML list also works via
 * -params-file. Each input must be a multiresolution cooler (a plain .cool
 * cannot select the resolution group), and the requested --ontad_resolution
 * must already exist under the file's /resolutions group. The assembled
 * <cool_id>.ontad.tsv files are published flat into --outdir.
 *
 * Reuses ../subworkflows/local/ontad unchanged; the process ext.args and
 * publishDir rules live in subpipeline/ontad.config.
 *
 * Nextflow adds <main-script-dir>/bin to the task PATH, which is
 * subpipeline/bin here because projectDir becomes subpipeline/. That directory
 * is a symlink to the project-root bin/ so the same helpers are reused (and not
 * duplicated). If you later switch to a container profile, check that the
 * symlink target is mounted; a real bin/ copy does not have that caveat.
 */

include { ONTAD_HIC } from '../subworkflows/local/ontad'

workflow {
    if (!params.mcool) {
        error("No mcool input specified: use --mcool <file|glob|comma-separated list>.")
    }
    if (!params.chrom_sizes) {
        error("No chrom.sizes specified: use --chrom_sizes <chrom.sizes>.")
    }
    if (!(params.ontad_resolution as int > 0)) {
        error("Invalid --ontad_resolution '${params.ontad_resolution}': expected a positive integer.")
    }

    //
    // --mcool is a path/glob, a comma-separated string, or a list (from
    // -params-file). Each input becomes one entry whose base name is the
    // output prefix. Materialize the list to fail early on empty globs or
    // colliding output names, then fan the entries back out as a queue
    // channel.
    //
    def mcool_paths = (params.mcool instanceof List)
        ? params.mcool.collect { item -> item.toString() }
        : params.mcool.toString().split(',').collect { token -> token.trim() }.findAll { token -> token }

    ch_cool = channel.fromPath(mcool_paths, checkIfExists: true)
        .map { path ->
            if (path.isDirectory()) {
                error("--mcool expects a cooler file, got directory: ${path}")
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
                error("Duplicate --mcool base name(s): ${duplicate_ids.join(', ')}. " +
                      "Each output is named <base>.ontad.tsv, so rename the inputs or process them separately.")
            }
            entries
        }
        .flatMap { entries -> entries }

    ch_chrom_sizes = channel.fromPath(params.chrom_sizes, checkIfExists: true)

    //
    // Same metadata recorded by the main pipeline's merged levels.
    //
    def ontad_metadata = [
        minsz  : params.ontad_minsz as int,
        maxsz  : params.ontad_maxsz as int,
        lsize  : params.ontad_lsize as int,
        ldiff  : params.ontad_ldiff as double,
        penalty: params.ontad_penalty as double,
    ]

    ONTAD_HIC(ch_cool, ch_chrom_sizes, params.ontad_resolution, ontad_metadata)
}
