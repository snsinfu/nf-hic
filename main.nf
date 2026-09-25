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
    // so CLI/config --fasta/--bwa_index/--chrom_sizes always win.
    //
    def fasta       = params.fasta       ?: getGenomeAttribute('fasta')
    def bwa_index   = params.bwa_index   ?: getGenomeAttribute(params.aligner == 'bwa' ? 'bwa' : 'bwamem2')
    def chrom_sizes = params.chrom_sizes ?: getGenomeAttribute('chrom_sizes')

    if (params.genome && params.aligner == 'bwa-mem2' && !getGenomeAttribute('bwamem2')) {
        log.warn("No bwa-mem2 index for genome '${params.genome}' in the catalog; it will be built from the fasta.")
    }

    HIC(fasta, bwa_index, chrom_sizes)
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
