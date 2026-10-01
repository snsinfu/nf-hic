process ONTAD_DUMP_MATRIX {
    tag "$meta.cool_id ${meta.chrom}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'quay.io/biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cool)
    val resolution

    output:
    tuple val(meta), path("${meta.cool_id}.${meta.chrom}.matrix"), emit: matrix
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ontad_dump_matrix.py \\
        --region '${meta.chrom}' \\
        --bin-size ${resolution} \\
        --output ${meta.cool_id}.${meta.chrom}.matrix \\
        $args \\
        ${cool}
    """

    stub:
    """
    touch ${meta.cool_id}.${meta.chrom}.matrix
    """
}
