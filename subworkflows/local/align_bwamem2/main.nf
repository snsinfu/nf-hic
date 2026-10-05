/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bwa-mem2 aligner adapter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uniform aligner interface: reads + index + fasta -> bam. One adapter per aligner
    keeps aligner-specific module signatures out of LIBRARY_HIC.
----------------------------------------------------------------------------------------
*/

include { BWAMEM2_MEM } from '../../../modules/local/bwamem2_mem'

workflow ALIGN_BWAMEM2 {
    take:
    reads      // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    index      // channel: [ val(meta), path(index) ]
    fasta      // channel: path(fasta)

    main:
    // mem holds the index resident, including the unpacked .0123, so pass its
    // on-disk footprint to the process (computed here: a path input is a relative
    // staged name inside a dynamic directive).
    ch_index = index.map { m, idx ->
        def resident = 0L
        idx.toFile().eachFileRecurse { f -> if (f.isFile()) resident += f.length() }
        [ m + [index_bytes: resident], idx ]
    }
    // nf-core aligner modules expect the fasta as a [ meta, fasta ] tuple
    BWAMEM2_MEM(reads, ch_index, fasta.map { fasta_path -> [ [:], fasta_path ] }, false)

    emit:
    bam = BWAMEM2_MEM.out.bam
}
