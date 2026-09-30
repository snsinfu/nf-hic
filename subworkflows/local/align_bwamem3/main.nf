/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bwa-mem3 aligner adapter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uniform aligner interface: reads + index + fasta -> bam. One adapter per aligner
    keeps aligner-specific module signatures out of LIBRARY_HIC.
----------------------------------------------------------------------------------------
*/

include { BWAMEM3_MEM } from '../../../modules/nf-core/bwamem3/mem'

workflow ALIGN_BWAMEM3 {
    take:
    reads      // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    index      // channel: [ val(meta), path(index) ]
    fasta      // channel: path(fasta)

    main:
    // nf-core aligner modules expect the fasta as a [ meta, fasta ] tuple
    BWAMEM3_MEM(reads, index, fasta.map { fasta_path -> [ [:], fasta_path ] }, false)

    emit:
    // nf-core bwamem3/mem emits the alignment under `aligned` (*.bam here), not `bam`
    bam = BWAMEM3_MEM.out.aligned
}
