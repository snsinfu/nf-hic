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
    def cat_cmd = gtf.toString().endsWith('.gz') ? "zcat" : "cat"
    """
    # Gene bodies (one interval per gene); skip GTF comment lines.
    ${cat_cmd} ${gtf} | awk -v FS='\\t' -v OFS='\\t' '!/^#/ && \$3 == "gene" { print \$1, \$4-1, \$5 }' \\
        | sort -k1,1 -k2,2n > genes.tsv
    test -s genes.tsv || { echo "ERROR: no 'gene' features found in ${gtf}; a gene-level GTF is required for the CALDER2 feature track" >&2; exit 1; }

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
