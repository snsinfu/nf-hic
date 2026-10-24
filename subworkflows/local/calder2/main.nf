/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Call CALDER2 compartments and nested sub-domains from a balanced multiresolution cooler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALDER2 } from '../../../modules/local/calder2'

workflow CALDER2_HIC {
    take:
    ch_cool            // channel: [ val(meta), path(mcool), path(bed) | [] ] ; meta.cool_id is the output basename
    resolution         // integer: resolution group inside the .mcool
    save_intermediates // boolean: keep and publish the full intermediate_data/ folder
    nproc              // integer: >0 overrides CALDER --nproc; 0 = task.cpus

    main:
    CALDER2(ch_cool, resolution, save_intermediates, nproc)

    emit:
    CALDER2.out.output_folder
}
