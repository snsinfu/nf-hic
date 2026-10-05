/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bwa-mem3 aligner adapter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uniform aligner interface: reads + index + fasta -> bam. One adapter per aligner
    keeps aligner-specific module signatures out of LIBRARY_HIC.
----------------------------------------------------------------------------------------
*/

include { BWAMEM3_MEM } from '../../../modules/local/bwamem3_mem'

workflow ALIGN_BWAMEM3 {
    take:
    reads      // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    index      // channel: [ val(meta), path(index) ]
    fasta      // channel: path(fasta)

    main:
    // mem holds the index resident, so pass its on-disk footprint to the process
    // (computed here: a path input is a relative staged name inside a dynamic
    // directive). bwa-mem3 pac-fetches from .pac and never reads .0123.
    ch_index = index.map { m, idx ->
        def resident = 0L
        idx.toFile().eachFileRecurse { f -> if (f.isFile() && !f.name.endsWith('.0123')) resident += f.length() }
        [ m + [index_bytes: resident], idx ]
    }
    // nf-core aligner modules expect the fasta as a [ meta, fasta ] tuple
    BWAMEM3_MEM(reads, ch_index, fasta.map { fasta_path -> [ [:], fasta_path ] }, false)

    emit:
    // nf-core bwamem3/mem emits the alignment under `aligned` (*.bam here), not `bam`
    bam = BWAMEM3_MEM.out.aligned
}
