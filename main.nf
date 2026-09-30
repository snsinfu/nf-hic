#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-hic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Minimal 4DN-style Hi-C pipeline: bwa/bwa-mem2/bwa-mem3 + pairtools + cooler.
----------------------------------------------------------------------------------------
*/

include { HIC } from './workflows/hic'

workflow {
    //
    // Resolve reference paths from the genome catalog only when not explicitly given,
    // so CLI/config --fasta/--bwa_index/--bwamem2_index/--bwamem3_index/--chrom_sizes/--gtf
    // always win.
    //
    def fasta       = params.fasta       ?: getGenomeAttribute('fasta')
    def chrom_sizes = params.chrom_sizes ?: getGenomeAttribute('chrom_sizes')
    def gtf         = params.gtf         ?: getGenomeAttribute('gtf')
    def index       = resolveAlignerIndex(params.aligner)

    //
    // Index directories are aligner-specific: bwa, bwa-mem2 and bwa-mem3 indexes are not
    // interchangeable, so each aligner reads its own parameter. Warn when another aligner's
    // parameter was supplied but is ignored.
    //
    def index_params = [
        'bwa'      : [flag: '--bwa_index',     value: params.bwa_index],
        'bwa-mem2' : [flag: '--bwamem2_index', value: params.bwamem2_index],
        'bwa-mem3' : [flag: '--bwamem3_index', value: params.bwamem3_index]
    ]
    def own_index = index_params[params.aligner]
    if (own_index && !own_index.value) {
        index_params.each { name, spec ->
            if (name != params.aligner && spec.value) {
                log.warn("${spec.flag} is ignored with --aligner ${params.aligner} (bwa, bwa-mem2 and bwa-mem3 indexes are not interchangeable). Use ${own_index.flag} for a precomputed ${params.aligner} index.")
            }
        }
    }
    if (params.genome && !index) {
        log.warn("No ${params.aligner} index for genome '${params.genome}' in the catalog; it will be built from the fasta.")
    }

    HIC(fasta, index, chrom_sizes, gtf, params.aligner)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get an attribute from the genome catalog, e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

//
// Resolve the precomputed index for an aligner. bwa, bwa-mem2 and bwa-mem3 indexes are not
// interchangeable, so each aligner reads only its own parameter (falling back to its
// own catalog key). Add a case here when adding a new aligner.
//
def resolveAlignerIndex(aligner) {
    if (aligner == 'bwa') {
        return params.bwa_index     ?: getGenomeAttribute('bwa')
    }
    else if (aligner == 'bwa-mem2') {
        return params.bwamem2_index ?: getGenomeAttribute('bwamem2')
    }
    else if (aligner == 'bwa-mem3') {
        return params.bwamem3_index ?: getGenomeAttribute('bwamem3')
    }
    return null
}
