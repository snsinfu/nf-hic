/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BWA aligner adapter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uniform aligner interface: reads + index + fasta -> bam. One adapter per aligner
    keeps aligner-specific module signatures out of LIBRARY_HIC.
----------------------------------------------------------------------------------------
*/

include { BWA_MEM } from '../../../modules/nf-core/bwa/mem'

workflow ALIGN_BWA {
    take:
    reads      // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    index      // channel: [ val(meta), path(index) ]
    fasta      // channel: path(fasta)

    main:
    // nf-core aligner modules expect the fasta as a [ meta, fasta ] tuple
    BWA_MEM(reads, index, fasta.map { fasta_path -> [ [:], fasta_path ] }, false)

    emit:
    bam = BWA_MEM.out.bam
}
