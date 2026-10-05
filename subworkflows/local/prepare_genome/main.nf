/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prepare reference genome files (aligner index and chromosome sizes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_INDEX      } from '../../../modules/nf-core/bwa/index'
include { BWAMEM2_INDEX  } from '../../../modules/local/bwamem2_index'
include { BWAMEM3_INDEX  } from '../../../modules/local/bwamem3_index'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx'

workflow PREPARE_GENOME {
    take:
    fasta         // path: genome fasta
    index         // path: prebuilt index directory for the selected aligner (optional)
    chrom_sizes   // path: chrom.sizes override (optional)
    aligner       // string: 'bwa' | 'bwa-mem2' | 'bwa-mem3' (add a case per new aligner)

    main:
    ch_fasta = channel.value(file(fasta, checkIfExists: true))

    //
    // Build the aligner index if one was not supplied
    //
    if (!index) {
        if (aligner == 'bwa') {
            BWA_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_index = BWA_INDEX.out.index
        }
        else if (aligner == 'bwa-mem2') {
            BWAMEM2_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_index = BWAMEM2_INDEX.out.index
        }
        else if (aligner == 'bwa-mem3') {
            BWAMEM3_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_index = BWAMEM3_INDEX.out.index
        }
        else {
            error("Invalid --aligner '${aligner}'. Use 'bwa', 'bwa-mem2' or 'bwa-mem3'.")
        }
    }
    else {
        ch_index = [ [:], file(index, checkIfExists: true) ]
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
    index       = ch_index
    chrom_sizes = ch_chrom_sizes
}
