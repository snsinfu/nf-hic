process CALDER2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-calder2:0.7--r43hdfd78af_1' :
        'quay.io/biocontainers/r-calder2:0.7--r43hdfd78af_1' }"

    input:
    tuple val(meta), path(cool)
    val resolution
    path feature_track
    val save_intermediates

    output:
    tuple val(meta), path("${prefix}/")                    , emit: output_folder
    tuple val(meta), path("${prefix}/intermediate_data/")  , emit: intermediate_data_folder, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('calder'), val('0.7'), emit: versions_calder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = resolution ? "::/resolutions/$resolution" : ""
    def cpus = task.cpus ?: 1
    // 4-column bed (chr, start, end, score), headerless; absent (`[]`) for built-in genomes
    def track_arg = feature_track ? "--feature_track ${feature_track}" : ''
    def save_arg = save_intermediates ? '--save-intermediates' : ''
    def resolution_arg = resolution ?: ''
    """
    mkdir -p ${prefix}/logs

    # getting binsize as mandatory input for calder
    binsize="\$(cooler info --field bin-size $cool$suffix)"

    # CALDER2 can exit 0 having produced nothing: build_comp_table_opt and
    # CALDER_sub_domains are wrapped in try() and the CLI ignores system()
    # failures. Record the exit code and the log instead of aborting, so a
    # failed analysis is a diagnosable result and does not stop the pipeline.
    # --keep_intermediate retains the logs/dumps for the post-processing below.
    rc=0
    calder --input $cool$suffix \\
        --outpath ${prefix} \\
        --nproc $cpus \\
        --type cool \\
        --bin_size "\${binsize}" \\
        --keep_intermediate \\
        $args \\
        $track_arg 2>&1 | tee ${prefix}/logs/calder.log || rc=\${PIPESTATUS[0]}
    echo "\${rc}" > ${prefix}/logs/calder.exit_code

    # Map CALDER's normalized chr names back to the cooler names and extract the
    # compact logs/ set; the large intermediate_data/ is dropped unless asked.
    calder2_postprocess.py \\
        --cool ${cool} \\
        --resolution "${resolution_arg}" \\
        --outdir ${prefix} \\
        $save_arg
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/sub_compartments
    mkdir -p ${prefix}/sub_domains
    mkdir -p ${prefix}/logs
    if [ "${save_intermediates}" == "true" ]; then
        mkdir -p ${prefix}/intermediate_data
    fi

    touch ${prefix}/sub_compartments/all_sub_compartments.bed
    touch ${prefix}/sub_compartments/all_sub_compartments.tsv
    touch ${prefix}/sub_compartments/cor_with_ref.ALL.txt
    touch ${prefix}/sub_compartments/cor_with_ref.pdf
    touch ${prefix}/sub_compartments/cor_with_ref.txt

    touch ${prefix}/sub_domains/all_nested_boundaries.bed

    touch ${prefix}/logs/calder.log
    touch ${prefix}/logs/calder.exit_code
    """
}
