nextflow.enable.dsl = 2

params.interval = "chr20:10000000-11000000"
params.threads = 2
params.outdir = "results/HG002_nextflow"

process PREPARE_GIAB {
    container "bioseq:latest"
    cache false

    input:
    val interval
    val threads

    output:
    path "prepared"

    script:
    """
    python -m bioseq.giab_data \\
      --interval '${interval}' \\
      --threads ${threads} \\
      --output-dir prepared
    """
}

process RUN_GIAB_BENCHMARK {
    container "bioseq:latest"
    publishDir params.outdir, mode: "copy"

    input:
    path prepared
    val interval
    val threads
    val start
    val end

    output:
    path "comparison"

    script:
    """
    python -m bioseq.full_benchmark \\
      --read1 '${prepared}/HG002_R1.fastq' \\
      --read2 '${prepared}/HG002_R2.fastq' \\
      --reference-fasta '${prepared}/chr20.fa' \\
      --known-sites-vcf '${prepared}/Mills.hg38.region.vcf.gz' \\
      --reference-sdf '${prepared}/chr20.sdf' \\
      --truth-vcf '${prepared}/HG002.truth.region.vcf.gz' \\
      --confident-regions '${prepared}/HG002_chr20_${start}_${end}.confident.bed' \\
      --interval '${interval}' \\
      --threads ${threads} \\
      --output-dir comparison
    """
}

workflow {
    def interval_match = params.interval =~ /^chr20:(\d+)-(\d+)$/
    if (!interval_match.matches()) {
        error "The HG002 preparation workflow supports intervals in chr20:start-end format."
    }

    def threads = params.threads as Integer
    if (threads < 1) {
        error "threads must be a positive integer."
    }

    def start = interval_match[0][1]
    def end = interval_match[0][2]
    PREPARE_GIAB(params.interval, threads)
    RUN_GIAB_BENCHMARK(PREPARE_GIAB.out, params.interval, threads, start, end)
}
