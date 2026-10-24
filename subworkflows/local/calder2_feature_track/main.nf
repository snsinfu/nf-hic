/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Resolve the CALDER2 A/B phasing feature track
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    One `--calder2_feature_track` value can carry a shared track, per-sample tracks, or
    both. The default is resolved once (explicit shared track > auto gene density for a
    non-built-in genome > CALDER's built-in reference), then attached to each cooler.
*/

include { CALDER2_GENE_DENSITY } from '../../../modules/local/calder2_gene_density'

//
// Parse `--calder2_feature_track` into [shared: <path|null>, perSample: [sample: path]].
//
// String form (CLI): comma-separated tokens; `key=path` is a per-sample track, a bare
// token is the shared track. Map form (-params-file/config): sample names to paths, with
// the single reserved key '*' for the shared track (literal key, not a glob). List form
// is deliberately not supported.
//
def addCalder2TrackEntry(state, key, value) {
    def k = key == null ? null : key.toString().trim()
    def v = value == null ? '' : value.toString().trim()
    if (k == '*') {
        if (state.shared != null) {
            error("--calder2_feature_track defines the shared track ('*') more than once.")
        }
        if (!v) {
            error("--calder2_feature_track: the shared track ('*') has no bed path.")
        }
        state.shared = v
    }
    else if (k == null || k == '') {
        error("--calder2_feature_track: empty sample name in an entry. Use a bare path for the shared track.")
    }
    else {
        if (state.perSample.containsKey(k)) {
            error("--calder2_feature_track: sample '${k}' is set more than once.")
        }
        if (!v) {
            error("--calder2_feature_track: sample '${k}' has no bed path.")
        }
        state.perSample[k] = v
    }
}

def parseCalder2FeatureTrack(spec) {
    def state = [shared: null, perSample: [:]]

    if (spec == null) {
        return state
    }
    if (spec instanceof Map) {
        spec.each { key, value -> addCalder2TrackEntry(state, key, value) }
    }
    else if (spec instanceof List) {
        error("Invalid --calder2_feature_track: list form is not supported. " +
              "Use a comma-separated string or a map.")
    }
    else {
        spec.toString().split(',').each { rawToken ->
            def token = rawToken.trim()
            if (token) {
                def i = token.indexOf('=')
                if (i < 0) {
                    addCalder2TrackEntry(state, '*', token)
                }
                else {
                    addCalder2TrackEntry(state, token.substring(0, i).trim(), token.substring(i + 1).trim())
                }
            }
        }
    }

    return state
}

//
// Resolve a track path and fail early with the offending key in the message.
//
def calder2TrackFile(path, label) {
    def track = file(path)
    if (!track.exists()) {
        error("--calder2_feature_track ${label} not found: ${path}")
    }
    return track
}

workflow CALDER2_FEATURE_TRACK {
    take:
    track_spec          // raw params.calder2_feature_track (string, map or null)
    genome              // effective CALDER genome name, for messages
    builtin             // boolean: hg19|hg38|mm9|mm10
    gtf                 // path or null; needed only for auto gene density
    chrom_sizes         // channel (main pipeline) or path/null (standalone)
    resolution          // integer: gene-density window / CALDER resolution

    main:
    def spec = parseCalder2FeatureTrack(track_spec)
    def perSampleFiles = spec.perSample.collectEntries { sample, path ->
        [sample, calder2TrackFile(path, "sample '${sample}'")]
    }

    if (builtin && (spec.shared || spec.perSample)) {
        log.warn "[nf-hic] --calder2_feature_track overrides CALDER's built-in reference for '${genome}' " +
                 "(CALDER treats the genome as 'others'; the built-in reference is not used)."
    }

    track = channel.empty()
    if (spec.shared) {
        track = channel.value(calder2TrackFile(spec.shared, 'shared track'))
    }
    else if (!builtin) {
        if (!gtf) {
            error("CALDER2 A/B phasing for genome '${genome}' needs a gene-density track, but no GTF annotation was found. " +
                  "Use --gtf <genes.gtf> or a --genome with a 'gtf' catalog entry.")
        }
        if (!chrom_sizes) {
            error("CALDER2 A/B phasing for genome '${genome}' needs --chrom_sizes to generate a gene-density track. " +
                  "Use --chrom_sizes <chrom.sizes> or a --genome with a 'chrom_sizes' catalog entry.")
        }
        CALDER2_GENE_DENSITY(
            chrom_sizes,
            channel.value(file(gtf, checkIfExists: true)),
            channel.value(resolution)
        )
        track = CALDER2_GENE_DENSITY.out.bed
    }
    else {
        track = channel.value([])
    }

    emit:
    track
    per_sample = channel.value(perSampleFiles)
}

workflow CALDER2_ATTACH_TRACK {
    take:
    ch_cool             // channel: [ val(meta), path(mcool) ]
    ch_default          // path(bed) or [] (value channel)
    ch_per_sample       // Map value channel: sample id -> path(bed)
    ch_known_ids        // list of valid sample ids (value channel)
    warn_unknown        // boolean: warn about unknown keys (only when ch_known_ids covers all samples)

    main:
    // Wrap the default and the id list in maps so an empty track ([]) stays a single
    // value and a List is not expanded by `combine`.
    ch_config = ch_default
        .map { track -> [fallback: track] }
        .combine(ch_per_sample)
        .combine(ch_known_ids.map { ids -> [ids: ids] })
        .map { wrapped, perSample, knownWrapper ->
            if (warn_unknown) {
                def known = knownWrapper.ids as Set
                def unknown = perSample.keySet().findAll { key -> !(key in known) }
                if (unknown) {
                    log.warn "[nf-hic] --calder2_feature_track key(s) ${unknown.join(', ')} " +
                             "do not match any sample id; they are ignored."
                }
            }
            [fallback: wrapped.fallback, perSample: perSample]
        }

    ch_cool_tracked = ch_cool
        .combine(ch_config)
        .map { meta, cool, config ->
            [meta, cool, config.perSample[meta.id] ?: config.fallback]
        }

    emit:
    ch_cool_tracked
}
