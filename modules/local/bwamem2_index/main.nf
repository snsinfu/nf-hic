process BWAMEM2_INDEX {
    tag "$fasta"
    // Local fork of nf-core/modules bwamem2/index @ efec54255f9baad3ea032173d75031929883bed8.
    // Only the memory directive differs.
    memory {
        // Expected memory footprint of the index build: ~28 B per genome base
        // (*4 when the FASTA is gzipped, *1.1 headroom above the build peak).
        def bases = fasta.size() * (fasta.name.endsWith('.gz') ? 4 : 1)
        def expected = (28L * bases).B * 1.1
        def baseLimit = [expected, 64.GB].max()
        (baseLimit * task.attempt)
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem2"), emit: index
    tuple val("${task.process}"), val('bwamem2'), eval('bwa-mem2 version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_bwamem2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir bwamem2
    bwa-mem2 \\
        index \\
        $args \\
        -p bwamem2/${prefix} \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"

    """
    mkdir bwamem2
    touch bwamem2/${prefix}.0123
    touch bwamem2/${prefix}.ann
    touch bwamem2/${prefix}.pac
    touch bwamem2/${prefix}.amb
    touch bwamem2/${prefix}.bwt.2bit.64
    """
}
