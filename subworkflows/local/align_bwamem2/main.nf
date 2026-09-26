/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bwa-mem2 aligner adapter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uniform aligner interface: reads + index + fasta -> bam. One adapter per aligner
    keeps aligner-specific module signatures out of LIBRARY_HIC.
----------------------------------------------------------------------------------------
*/

include { BWAMEM2_MEM } from '../../../modules/nf-core/bwamem2/mem'

workflow ALIGN_BWAMEM2 {
    take:
    reads      // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    index      // channel: [ val(meta), path(index) ]
    fasta      // channel: path(fasta)

    main:
    // nf-core aligner modules expect the fasta as a [ meta, fasta ] tuple
    BWAMEM2_MEM(reads, index, fasta.map { fasta_path -> [ [:], fasta_path ] }, false)

    emit:
    bam = BWAMEM2_MEM.out.bam
}
