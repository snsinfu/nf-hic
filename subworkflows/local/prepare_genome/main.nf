/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prepare reference genome files (aligner index and chromosome sizes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_INDEX      } from '../../../modules/nf-core/bwa/index'
include { BWAMEM2_INDEX  } from '../../../modules/nf-core/bwamem2/index'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx'

workflow PREPARE_GENOME {
    take:
    fasta         // path: genome fasta
    bwa_index     // path: prebuilt index directory for the selected aligner (optional)
    chrom_sizes   // path: chrom.sizes override (optional)

    main:
    ch_fasta = channel.value(file(fasta, checkIfExists: true))

    //
    // Build the aligner index if one was not supplied
    //
    if (!bwa_index) {
        if (params.aligner == 'bwa') {
            BWA_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_bwa_index = BWA_INDEX.out.index
        }
        else if (params.aligner == 'bwa-mem2') {
            BWAMEM2_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_bwa_index = BWAMEM2_INDEX.out.index
        }
        else {
            error("Invalid --aligner '${params.aligner}'. Use 'bwa' or 'bwa-mem2'.")
        }
    }
    else {
        ch_bwa_index = [ [:], file(bwa_index, checkIfExists: true) ]
    }

    //
    // Chromosome sizes: use the supplied file, otherwise derive from the fasta
    //
    if (!chrom_sizes) {
        SAMTOOLS_FAIDX(ch_fasta.map { item -> [ [:], item, [] ] }, true)
        ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { item -> item[1] }
    }
    else {
        ch_chrom_sizes = channel.value(file(chrom_sizes, checkIfExists: true))
    }

    emit:
    fasta       = ch_fasta
    bwa_index   = ch_bwa_index
    chrom_sizes = ch_chrom_sizes
}
