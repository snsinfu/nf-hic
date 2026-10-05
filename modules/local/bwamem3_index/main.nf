process BWAMEM3_INDEX {
    tag "${meta.id}"
    // Local fork of nf-core/modules bwamem3/index @ efec54255f9baad3ea032173d75031929883bed8.
    // Only the memory directive differs. bwa-mem3 0.8.0 index peak: 8 bytes per
    // base of doubled text (= 16 B/genome-base) below ~1.07 Gbp (32-bit SA), else
    // 12 bytes per doubled base (= 24 B/genome-base). Source:
    // https://bwa-mem3.readthedocs.io/en/latest/user-guide/indexing.html
    // The tool budgets --max-memory as detected RAM less min(max(2GB, 5%), 50%)
    // and rejects the build up front when that is below the requirement, so
    // request the requirement plus headroom.
    memory {
        def bases = fasta.size() * (fasta.name.endsWith('.gz') ? 4 : 1)
        def need  = (bases < 1070000000L ? 16L : 24L) * bases
        def headroom = Math.max(2.5 * 1024 ** 3, 0.10 * need)
        ((need + headroom) * task.attempt).toLong().B
    }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/09/09c0e84ac22ba722bae6f75c5c450d72dfdc06e7224ad2c0bd042849632ca70e/data'
        : 'community.wave.seqera.io/library/bwa-mem3:0.8.0--a430e1e1df49a92e'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem3"), emit: index
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir bwamem3
    bwa-mem3 \\
        index \\
        $args \\
        -t ${task.cpus} \\
        -p bwamem3/${prefix} \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwamem3
    touch bwamem3/${prefix}.0123
    touch bwamem3/${prefix}.amb
    touch bwamem3/${prefix}.ann
    touch bwamem3/${prefix}.bwt.2bit.64
    touch bwamem3/${prefix}.pac
    """
}
