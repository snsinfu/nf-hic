process ONTAD {
    tag "$meta.cool_id ${meta.chrom}"
    label 'process_high'

    container 'ghcr.io/snsinfu/anlin00007-ontad:v1.4-p3'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${meta.cool_id}.${meta.chrom}.tad"), emit: tad
    // WARN: OnTAD has no --version flag; keep in sync with the container tag.
    tuple val("${task.process}"), val('OnTAD'), val('1.4-p3'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    OnTAD ${matrix} \\
        -o ${meta.cool_id}.${meta.chrom} \\
        $args
    """

    stub:
    """
    touch ${meta.cool_id}.${meta.chrom}.tad
    """
}
