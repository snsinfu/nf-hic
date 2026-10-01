/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Call OnTAD hierarchical TADs from a balanced multiresolution cooler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OnTAD reads one dense N*N per-chromosome text matrix. The matrices are never
    published; only the assembled genome-wide TSV is.
*/

include { ONTAD_DUMP_MATRIX } from '../../../modules/local/ontad_dump_matrix'
include { ONTAD             } from '../../../modules/local/ontad'
include { ONTAD_ASSEMBLE    } from '../../../modules/local/ontad_assemble'

workflow ONTAD_HIC {
    take:
    ch_cool         // channel: [ val(meta), path(mcool) ] ; meta.cool_id is the output basename
    ch_chrom_sizes  // channel: path(chrom.sizes) ; determines which chromosomes to process
    resolution      // integer: resolution group inside the .mcool
    metadata        // map: key/value metadata recorded in the output TSV

    main:
    //
    // One entry per chromosome, as defined by chrom.sizes
    //
    ch_chrom = ch_chrom_sizes
        .splitText()
        .map { line -> line.trim().split(/\s+/)[0] }
        .filter { name -> name }

    ONTAD_DUMP_MATRIX(
        ch_cool.combine(ch_chrom).map { meta, cool, chrom -> [ meta + [chrom: chrom], cool ] },
        resolution
    )

    ONTAD(ONTAD_DUMP_MATRIX.out.matrix)

    //
    // Group the per-chromosome .tad files of each cooler and attach the cooler
    // so the assembler can recover the bins
    //
    ONTAD_ASSEMBLE(
        ONTAD.out.tad
            .map { meta, tad -> [ [id: meta.id, cool_id: meta.cool_id], tad ] }
            .groupTuple(by: [0])
            .join(ch_cool.map { meta, cool -> [ [id: meta.id, cool_id: meta.cool_id], cool ] })
            .map { key, tads, cool -> [ key, cool, resolution, metadata, tads ] }
    )

    emit:
    ONTAD_ASSEMBLE.out.tads
}
