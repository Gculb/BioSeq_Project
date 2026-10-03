nextflow.enable.dsl = 2

params.outdir = "results/rnaseq_airway"
params.analysis_script = "${projectDir}/analysis/rnaseq_airway_deseq2.R"

process DIFFERENTIAL_EXPRESSION {
    container "bioseq-rnaseq:latest"
    publishDir params.outdir, mode: "copy"

    input:
    path analysis_script

    output:
    path "deseq2_results.csv"
    path "analysis_summary.txt"
    path "ma_plot.png"
    path "volcano_plot.png"
    path "session_info.txt"

    script:
    """
    Rscript ${analysis_script}
    """
}

workflow {
    DIFFERENTIAL_EXPRESSION(file(params.analysis_script))
}
