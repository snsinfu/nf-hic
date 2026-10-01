process ONTAD_ASSEMBLE {
    tag "$meta.cool_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'quay.io/biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cool), val(resolution), val(metadata), path(tad_files)

    output:
    tuple val(meta), path("${meta.cool_id}.ontad.tsv"), emit: tads
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def metadata_json = groovy.json.JsonOutput.toJson(metadata ?: [:])
    """
    ontad_assemble.py \\
        --cool ${cool} \\
        --bin-size ${resolution} \\
        --prefix ${meta.cool_id} \\
        --output ${meta.cool_id}.ontad.tsv \\
        --metadata '${metadata_json}' \\
        $args \\
        ${tad_files}
    """

    stub:
    """
    touch ${meta.cool_id}.ontad.tsv
    """
}
