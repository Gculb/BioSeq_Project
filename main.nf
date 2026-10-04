nextflow.enable.dsl = 2

params.outdir = "results/rnaseq_airway_raw"
params.samplesheet = "${projectDir}/analysis/airway_samples.tsv"
params.threads = 2
params.ensembl_release = "112"

process DOWNLOAD_READS {
    container "bioseq-rnaseq:latest"
    cpus 1
    maxForks 1

    input:
    tuple val(meta), val(run)

    output:
    tuple val(meta), path("read_1.fastq.gz"), path("read_2.fastq.gz")

    script:
    """
    curl --retry 3 --retry-delay 5 -fsSL \\
      'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp&format=tsv' \\
      -o run_info.tsv
    urls=\$(awk -F '\\t' 'NR == 2 { print \$2 }' run_info.tsv | tr ';' '\\n')
    read1=\$(printf '%s\\n' "\$urls" | grep -E '_1\\.fastq\\.gz' | head -n 1)
    read2=\$(printf '%s\\n' "\$urls" | grep -E '_2\\.fastq\\.gz' | head -n 1)
    if [ -z "\$read1" ] || [ -z "\$read2" ]; then
      echo "ENA did not provide paired FASTQ files for ${run}" >&2
      exit 1
    fi
    curl --retry 3 --retry-delay 5 -fsSL "https://\$read1" -o read_1.fastq.gz
    curl --retry 3 --retry-delay 5 -fsSL "https://\$read2" -o read_2.fastq.gz
    gzip -t read_1.fastq.gz
    gzip -t read_2.fastq.gz
    """
}

process FASTP {
    container "bioseq-rnaseq:latest"
    cpus params.threads
    maxForks 1
    publishDir "${params.outdir}/qc", mode: "copy", pattern: "*fastp.*"

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta),
          path("${meta.sample}_R1.trimmed.fastq.gz"),
          path("${meta.sample}_R2.trimmed.fastq.gz"),
          emit: trimmed
    path "${meta.sample}.fastp.html", emit: html
    path "${meta.sample}.fastp.json", emit: json

    script:
    """
    fastp \\
      --in1 ${read1} \\
      --in2 ${read2} \\
      --out1 ${meta.sample}_R1.trimmed.fastq.gz \\
      --out2 ${meta.sample}_R2.trimmed.fastq.gz \\
      --thread ${task.cpus} \\
      --html ${meta.sample}.fastp.html \\
      --json ${meta.sample}.fastp.json
    """
}

process SALMON_INDEX {
    container "bioseq-rnaseq:latest"
    cpus params.threads

    output:
    path "salmon_index", emit: index
    path "tx2gene.tsv", emit: tx2gene

    script:
    """
    base='https://ftp.ensembl.org/pub/release-${params.ensembl_release}/fasta/homo_sapiens/cdna'
    annotation='https://ftp.ensembl.org/pub/release-${params.ensembl_release}/gtf/homo_sapiens'
    curl --retry 3 --retry-delay 5 -fsSL \\
      "\${base}/Homo_sapiens.GRCh38.cdna.all.fa.gz" -o transcripts.fa.gz
    curl --retry 3 --retry-delay 5 -fsSL \\
      "\${annotation}/Homo_sapiens.GRCh38.${params.ensembl_release}.gtf.gz" -o annotation.gtf.gz
    gunzip -f transcripts.fa.gz annotation.gtf.gz
    salmon index -t transcripts.fa -i salmon_index -k 31
    awk -F '\\t' '
      \$3 == "transcript" {
        gene = \$9
        transcript = \$9
        sub(/.*gene_id "/, "", gene)
        sub(/".*/, "", gene)
        sub(/.*transcript_id "/, "", transcript)
        sub(/".*/, "", transcript)
        if (gene != "" && transcript != "") print transcript "\\t" gene
      }
    ' annotation.gtf | sort -u > tx2gene.tsv
    test -s tx2gene.tsv
    """
}

process SALMON_QUANT {
    container "bioseq-rnaseq:latest"
    cpus params.threads
    maxForks 1
    publishDir "${params.outdir}/quant", mode: "copy", pattern: "*.quant.sf"

    input:
    tuple val(meta), path(read1), path(read2), path(salmon_index)

    output:
    tuple val(meta),
          path("${meta.sample}.quant.sf"),
          path("${meta.sample}_salmon"),
          emit: quant

    script:
    """
    salmon quant \\
      -i ${salmon_index} \\
      -l A \\
      -1 ${read1} \\
      -2 ${read2} \\
      --validateMappings \\
      --seqBias \\
      --gcBias \\
      -p ${task.cpus} \\
      -o quant
    cp quant/quant.sf ${meta.sample}.quant.sf
    mkdir -p ${meta.sample}_salmon/aux_info
    cp quant/aux_info/meta_info.json ${meta.sample}_salmon/aux_info/
    if [ -f quant/aux_info/lib_format_counts.json ]; then
      cp quant/aux_info/lib_format_counts.json ${meta.sample}_salmon/aux_info/
    fi
    if [ -f quant/aux_info/flenDist.txt ]; then
      cp quant/aux_info/flenDist.txt ${meta.sample}_salmon/aux_info/
    fi
    """
}

process MULTIQC {
    container "multiqc/multiqc:v1.27@sha256:b4c30167e87ecb120925792c20c4c98d06076dde5f6ba9af48dda792fa549225"
    publishDir "${params.outdir}/multiqc", mode: "copy"

    input:
    path fastp_reports
    path salmon_reports

    output:
    path "multiqc_report.html"

    script:
    """
    mkdir -p multiqc_input
    cp -L *.fastp.json multiqc_input/
    cp -LR *_salmon multiqc_input/
    multiqc multiqc_input --outdir multiqc_output --filename multiqc_report.html --force --no-data-dir
    cp multiqc_output/multiqc_report.html .
    """
}

process DESEQ2 {
    container "bioseq-rnaseq:latest"
    publishDir "${params.outdir}/results", mode: "copy"

    input:
    tuple val(sample_meta), path(quant_files)
    path tx2gene
    path analysis_script

    output:
    path "deseq2_results.csv"
    path "analysis_summary.txt"
    path "ma_plot.png"
    path "volcano_plot.png"
    path "pca_plot.png"
    path "session_info.txt"
    script:
    def sample_manifest = (0..<sample_meta.size()).collect { index ->
        def meta = sample_meta[index]
        "${meta.sample}\t${meta.cell}\t${meta.dex}\t${quant_files[index]}"
    }.join("\n")
    """
    printf 'sample\\tcell\\tdex\\tquant_file\\n%s\\n' '${sample_manifest}' > quant_manifest.tsv
    Rscript ${analysis_script} quant_manifest.tsv ${tx2gene}
    """
}

workflow {
    def samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
    def salmon_samples = samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                [sample: row.sample, cell: row.cell, dex: row.dex],
                row.run
            )
        }

    DOWNLOAD_READS(salmon_samples)
    FASTP(DOWNLOAD_READS.out)
    SALMON_INDEX()
    salmon_index = SALMON_INDEX.out.index
    tx2gene = SALMON_INDEX.out.tx2gene
    SALMON_QUANT(FASTP.out.trimmed.combine(salmon_index))
    MULTIQC(
        FASTP.out.json.collect(),
        SALMON_QUANT.out.quant
            .map { meta, quant_file, salmon_report -> salmon_report }
            .collect()
    )

    quant_results = SALMON_QUANT.out.quant
        .map { meta, quant_file, salmon_report -> tuple(meta, quant_file) }
        .collect(flat: false)
        .map { results ->
            tuple(
                results.collect { meta, quant_file -> meta },
                results.collect { meta, quant_file -> quant_file }
            )
        }
    DESEQ2(
        quant_results,
        tx2gene,
        file("${projectDir}/analysis/rnaseq_airway_deseq2.R")
    )
}
