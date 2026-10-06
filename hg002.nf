nextflow.enable.dsl = 2

params.interval = "chr20:10000000-11000000"
params.threads = 2
params.outdir = "results/HG002_nextflow"
params.prepared_dir = null

process PREPARE_GIAB {
    container "seqcheckflow:latest"
    cache false

    input:
    val interval
    val threads

    output:
    path "prepared"

    script:
    """
    python -m seqcheckflow.giab_data \\
      --interval '${interval}' \\
      --threads ${threads} \\
      --output-dir prepared
    """
}

process FASTQ_QC {
    container "seqcheckflow:latest"
    publishDir "${params.outdir}/qc", mode: "copy"

    input:
    tuple path(read1), path(read2)

    output:
    path "read1_qc"
    path "read2_qc"

    script:
    """
    python -m seqcheckflow.pipeline ${read1} --results-dir read1_qc
    python -m seqcheckflow.pipeline ${read2} --results-dir read2_qc
    """
}

process SETUP_REFERENCE {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    path prepared

    output:
    path "reference_bundle"

    script:
    """
    mkdir reference_bundle
    cp ${prepared}/chr20.fa reference_bundle/reference.fa
    cp -a ${prepared}/chr20.sdf reference_bundle/reference.sdf
    cp ${prepared}/Mills.hg38.region.vcf.gz* reference_bundle/
    cp ${prepared}/HG002.truth.region.vcf.gz reference_bundle/truth.vcf.gz
    cp ${prepared}/HG002_chr20_*.confident.bed reference_bundle/confident_regions.bed
    bwa index reference_bundle/reference.fa
    samtools faidx reference_bundle/reference.fa
    gatk CreateSequenceDictionary \\
      -R reference_bundle/reference.fa \\
      -O reference_bundle/reference.dict
    """
}

process BWA_ALIGN {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    tuple path(read1), path(read2)
    path reference_bundle

    output:
    path "alignment.sam"

    script:
    def read_group = "@RG\\tID:HG002\\tSM:HG002\\tPL:ILLUMINA"
    """
    bwa mem \\
      -R '${read_group}' \\
      -t ${task.cpus} \\
      ${reference_bundle}/reference.fa \\
      ${read1} ${read2} > alignment.sam
    """
}

process SAM_TO_BAM {
    container "seqcheckflow:latest"

    input:
    path alignment_sam

    output:
    path "alignment.bam"

    script:
    """
    samtools view -b -o alignment.bam ${alignment_sam}
    """
}

process NAME_SORT {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    path alignment_bam

    output:
    path "alignment.name_sorted.bam"

    script:
    """
    samtools sort -@ ${task.cpus} -n -O BAM \\
      -o alignment.name_sorted.bam ${alignment_bam}
    """
}

process FIXMATE {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    path name_sorted_bam

    output:
    path "alignment.fixmate.bam"

    script:
    """
    samtools fixmate -@ ${task.cpus} -m -O BAM \\
      ${name_sorted_bam} alignment.fixmate.bam
    """
}

process COORDINATE_SORT {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    path fixmate_bam

    output:
    path "alignment.sorted.bam"

    script:
    """
    samtools sort -@ ${task.cpus} -O BAM \\
      -o alignment.sorted.bam ${fixmate_bam}
    """
}

process MARK_DUPLICATES {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    path coordinate_sorted_bam

    output:
    path "marked"

    script:
    """
    mkdir marked
    samtools markdup -@ ${task.cpus} \\
      ${coordinate_sorted_bam} marked/alignment.deduplicated.bam
    samtools index -@ ${task.cpus} \\
      marked/alignment.deduplicated.bam marked/alignment.deduplicated.bam.bai
    samtools coverage -r '${params.interval}' \\
      marked/alignment.deduplicated.bam > marked/alignment_coverage.txt
    """
}

process BASE_RECALIBRATOR {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    path marked
    path reference_bundle

    output:
    path "bqsr.table"

    script:
    """
    gatk BaseRecalibrator \\
      -R ${reference_bundle}/reference.fa \\
      -I ${marked}/alignment.deduplicated.bam \\
      --known-sites ${reference_bundle}/Mills.hg38.region.vcf.gz \\
      -L '${params.interval}' \\
      -O bqsr.table
    """
}

process APPLY_BQSR {
    container "seqcheckflow:latest"
    cpus params.threads

    input:
    path marked
    path bqsr_table
    path reference_bundle

    output:
    path "recalibrated"

    script:
    """
    mkdir recalibrated
    gatk ApplyBQSR \\
      -R ${reference_bundle}/reference.fa \\
      -I ${marked}/alignment.deduplicated.bam \\
      --bqsr-recal-file ${bqsr_table} \\
      -L '${params.interval}' \\
      -O recalibrated/alignment.recalibrated.bam
    samtools index -@ ${task.cpus} \\
      recalibrated/alignment.recalibrated.bam \\
      recalibrated/alignment.recalibrated.bam.bai
    """
}

process HAPLOTYPE_CALLER {
    container "seqcheckflow:latest"
    cpus params.threads
    publishDir "${params.outdir}/variants", mode: "copy"

    input:
    path recalibrated
    path reference_bundle

    output:
    path "calls.vcf.gz"

    script:
    """
    gatk HaplotypeCaller \\
      -R ${reference_bundle}/reference.fa \\
      -I ${recalibrated}/alignment.recalibrated.bam \\
      -O calls.vcf.gz \\
      -L '${params.interval}' \\
      --native-pair-hmm-threads ${task.cpus}
    """
}

process VCF_QC {
    container "seqcheckflow:latest"
    publishDir "${params.outdir}/qc", mode: "copy"

    input:
    path calls_vcf

    output:
    path "vcf_qc"

    script:
    """
    python -m seqcheckflow.pipeline ${calls_vcf} --results-dir vcf_qc
    """
}

process SCORE_VARIANTS {
    container "seqcheckflow:latest"
    publishDir "${params.outdir}/truth_evaluation", mode: "copy"

    input:
    path calls_vcf
    path reference_bundle

    output:
    path "rtg_eval"

    script:
    """
    rtg vcfeval \\
      --baseline ${reference_bundle}/truth.vcf.gz \\
      --calls ${calls_vcf} \\
      --template ${reference_bundle}/reference.sdf \\
      --evaluation-regions ${reference_bundle}/confident_regions.bed \\
      --output rtg_eval \\
      --no-roc
    test -s rtg_eval/summary.txt
    """
}

workflow {
    def interval_match = params.interval =~ /^chr20:(\d+)-(\d+)$/
    if (!interval_match.matches()) {
        error "The HG002 workflow supports intervals in chr20:start-end format."
    }

    def threads = params.threads as Integer
    if (threads < 1) {
        error "threads must be a positive integer."
    }

    if (params.prepared_dir) {
        prepared = Channel.fromPath(params.prepared_dir, checkIfExists: true)
    } else {
        prepared = PREPARE_GIAB(params.interval, threads)
    }
    FASTQ_QC(prepared.map { data -> tuple(file("${data}/HG002_R1.fastq"),
                                           file("${data}/HG002_R2.fastq")) })
    reference = SETUP_REFERENCE(prepared)
    reference_bundle = reference

    alignment = BWA_ALIGN(
        prepared.map { data -> tuple(file("${data}/HG002_R1.fastq"),
                                     file("${data}/HG002_R2.fastq")) },
        reference_bundle
    )
    converted = SAM_TO_BAM(alignment)
    name_sorted = NAME_SORT(converted)
    fixmate = FIXMATE(name_sorted)
    coordinate_sorted = COORDINATE_SORT(fixmate)
    marked = MARK_DUPLICATES(coordinate_sorted)
    recalibration = BASE_RECALIBRATOR(marked, reference_bundle)
    recalibrated = APPLY_BQSR(marked, recalibration, reference_bundle)
    calls = HAPLOTYPE_CALLER(recalibrated, reference_bundle)

    VCF_QC(calls)
    SCORE_VARIANTS(calls, reference_bundle)
}
