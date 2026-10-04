#!/usr/bin/env nextflow
/*
 * Standalone OnTAD sub-pipeline: re-run the OnTAD subworkflow
 * (dump -> OnTAD -> assemble) on existing balanced multiresolution coolers,
 * without repeating alignment or pairtools/cooler construction.
 *
 * This is the second entry point of the nf-hic repository; the full pipeline
 * is main.nf. Both share modules/, subworkflows/, workflows/ and conf/.
 *
 * Usage, remote (published repository):
 *
 *   nextflow run snsinfu/nf-hic \
 *       -main-script ontad.nf \
 *       --ontad_mcool 'results/bwa/merged_library/*.mcool' \
 *       --chrom_sizes _data/hg38/chrom.sizes \
 *       --ontad_resolution 10000
 *
 * Usage, local checkout (from the repository root):
 *
 *   nextflow run ontad.nf \
 *       --ontad_mcool 'results/bwa/merged_library/*.mcool' \
 *       --chrom_sizes _data/hg38/chrom.sizes \
 *       --ontad_resolution 10000
 *
 * Nextflow loads the root nextflow.config automatically (including
 * conf/modules.config and conf/base.config), so no -C/-c config file is
 * needed. A bare `snsinfu/nf-hic/ontad.nf` project coordinate is rejected by
 * Nextflow ("Repository URL must not end with a script file extension"); use
 * -main-script for remote runs.
 *
 * Note: CLI -C/-c config paths are resolved against the *working directory*,
 * not projectDir, so a repo-relative `-C subpipeline/ontad.config` cannot be
 * used remotely. Keeping this script at the repo root and reusing
 * nextflow.config avoids that problem entirely.
 *
 * Values that start with '-' (e.g. --ontad_args '-log2') need the
 * --option=value form, otherwise Nextflow's CLI reads them as options.
 *
 * One or more mcool paths/globs may be given. --ontad_mcool takes a single
 * path/glob or a comma-separated list (repeated --ontad_mcool flags do not
 * accumulate in Nextflow: the last one wins). A YAML list also works via
 * -params-file. Each input must be a multiresolution cooler (a plain .cool
 * cannot select the resolution group), and the requested --ontad_resolution
 * must already exist under the file's /resolutions group. The assembled
 * <cool_id>.ontad.tsv files are published to <outdir>/ontad/.
 *
 * projectDir is the repo root here, so the shared bin/ helpers are used
 * directly (no per-subpipeline symlink).
 */

include { ONTAD_HIC } from './subworkflows/local/ontad'

workflow {
    if (!params.ontad_mcool) {
        error("No mcool input specified: use --ontad_mcool <file|glob|comma-separated list>.")
    }
    if (!params.chrom_sizes) {
        error("No chrom.sizes specified: use --chrom_sizes <chrom.sizes>.")
    }
    if (!(params.ontad_resolution as int > 0)) {
        error("Invalid --ontad_resolution '${params.ontad_resolution}': expected a positive integer.")
    }

    //
    // --ontad_mcool is a path/glob, a comma-separated string, or a list (from
    // -params-file). Each input becomes one entry whose base name is the
    // output prefix. Materialize the list to fail early on empty globs or
    // colliding output names, then fan the entries back out as a queue
    // channel.
    //
    def mcool_paths = (params.ontad_mcool instanceof List)
        ? params.ontad_mcool.collect { item -> item.toString() }
        : params.ontad_mcool.toString().split(',').collect { token -> token.trim() }.findAll { token -> token }

    ch_cool = channel.fromPath(mcool_paths, checkIfExists: true)
        .map { path ->
            if (path.isDirectory()) {
                error("--ontad_mcool expects a cooler file, got directory: ${path}")
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
                error("Duplicate --ontad_mcool base name(s): ${duplicate_ids.join(', ')}. " +
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
