rule PicardCollectMultipleMetrics:
    input:
        bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam",
        ref=config["genome"]
    output:
        # Through the output file extensions the different tools for the metrics can be selected
        # so that it is not necessary to specify them under params with the "PROGRAM" option.
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html.
        multiext(
            "qc/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".quality_distribution_metrics",
            ".quality_distribution.pdf",
            ".quality_by_cycle_metrics",
            ".quality_by_cycle.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
            ".rna_metrics",
            ".bait_bias_detail_metrics",
            ".bait_bias_summary_metrics",
            ".error_summary_metrics",
            ".pre_adapter_detail_metrics",
            ".pre_adapter_summary_metrics",
            ".quality_yield_metrics",
        ),
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=4096,
    log:
        "logs/{sample}.PicardCollectMultipleMetrics.log",
    params:
        # optional parameters
        # REF_FLAT is required if RnaSeqMetrics are used
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE --REF_FLAT resources/refFlat.txt.gz",
    wrapper:
        "v3.13.1/bio/picard/collectmultiplemetrics"

rule PicardLibraryComplexity:
    input:
        config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.bam"
    output:
        "qc/{sample}.est_lib_complex_metrics.txt"
    threads: 1
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/{sample}.PicardLibraryComplexity.log"
    shell:
        "picard EstimateLibraryComplexity I={input} O={output} 2>{logs}"