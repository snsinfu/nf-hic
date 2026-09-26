#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-hic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Minimal 4DN-style Hi-C pipeline: bwa/bwa-mem2 + pairtools + cooler.
----------------------------------------------------------------------------------------
*/

include { HIC } from './workflows/hic'

workflow {
    //
    // Resolve reference paths from the genome catalog only when not explicitly given,
    // so CLI/config --fasta/--bwa_index/--bwamem2_index/--chrom_sizes always win.
    //
    def fasta       = params.fasta       ?: getGenomeAttribute('fasta')
    def chrom_sizes = params.chrom_sizes ?: getGenomeAttribute('chrom_sizes')

    //
    // Index directories are aligner-specific: bwa and bwa-mem2 indexes are not
    // interchangeable, so each aligner reads its own parameter (and catalog key).
    //
    def index
    if (params.aligner == 'bwa') {
        index = params.bwa_index     ?: getGenomeAttribute('bwa')
    }
    else {
        index = params.bwamem2_index ?: getGenomeAttribute('bwamem2')
    }

    if (params.aligner == 'bwa-mem2' && !params.bwamem2_index && params.bwa_index) {
        log.warn("--bwa_index is ignored with --aligner bwa-mem2 (bwa and bwa-mem2 indexes are not interchangeable). Use --bwamem2_index for a precomputed bwa-mem2 index.")
    }
    if (params.aligner == 'bwa' && !params.bwa_index && params.bwamem2_index) {
        log.warn("--bwamem2_index is ignored with --aligner bwa. Use --bwa_index for a precomputed bwa index.")
    }
    if (params.genome && params.aligner == 'bwa-mem2' && !index) {
        log.warn("No bwa-mem2 index for genome '${params.genome}' in the catalog; it will be built from the fasta.")
    }

    HIC(fasta, index, chrom_sizes)
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
