process CALDER2_GENE_DENSITY {
    tag "gene_density"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path chrom_sizes
    path gtf
    val  window

    output:
    path("*.bed"), emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "calder2_gene_density"
    def win = window ?: 100000
    """
    # One gene body per gene: prefer 'gene' features, else collapse features by
    # gene_id (supports GTFs without 'gene' rows, e.g. UCSC ncbiRefSeq).
    gtf_to_gene_bodies.sh ${gtf} > genes.tsv

    bedtools makewindows -g ${chrom_sizes} -w ${win} > windows.tsv

    # Headerless 4-column bed: chr, start, end, gene count
    bedtools intersect -a windows.tsv -b genes.tsv -c > ${prefix}.bed

    awk -v FS='\\t' '\$4 > 0 { n++ } END { if (n == 0) { print "ERROR: no windows overlap any gene; check that the GTF and chrom.sizes use identical chromosome names" > "/dev/stderr"; exit 1 } }' ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "calder2_gene_density"
    """
    touch ${prefix}.bed
    """
}
