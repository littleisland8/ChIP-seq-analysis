rule ScaleFactor:
	input:
		"qc/{sample}.bwa.flagstat"
	output:
		"results/{sample}.scalefactor.txt"
	log:
		"logs/{sample}.ScaleFactor.log"
	threads: 1
	shell:
		"""
		SCALE_FACTOR=$(grep '[0-9] mapped (' {input} | awk '{{print 10000000/$1}}') && echo ${{SCALE_FACTOR}} > {output} 2>{log}
		"""

rule BedtoolsGenomeCov:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam",
		ScaleFactor="results/{sample}.scalefactor.txt"
	output:
		"results/{sample}.bedGraph"
	log:
		"logs/{sample}.BedtoolsGenomeCov.log"
	threads: 1
	conda:
		"../envs/bedtools.yaml"
	shell:
		"""
		bedtools genomecov -ibam {input.bam} -scale {input.ScaleFactor} -bg -pc -fs 300 | sort -T '.' -k1,1 -k2,2n > {output} 2>{log} #-scale {input.ScaleFactor}
		"""

rule ChromSizes:
	input:
		"/resources/Transcriptomestic/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
	output:
		"resources/chromo_sizes.txt"
	log:
		"logs/ChromSizes.log"
	threads:1
	shell:
		"cat {input} | cut -f1,2 > {output} 2>{log}"

rule BedGraphToBigWig:
	input:
		bedgraph="results/{sample}.bedGraph",
		chromo_sizes="resources/chromo_sizes.txt"
	output:
		"results/{sample}.bigWig"
	log:
		"logs/{sample}.BedGraphToBigWig.log"
	threads: 1
	conda:
		"../envs/BedGraphToBigWig.yaml"
	shell:
		"bedGraphToBigWig {input.bedgraph} {input.chromo_sizes} {output} 2>{log}"

rule GenerateBigWig:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam.bai"
	output:
		"results/{sample}.bwa.bw",
	log:
		"logs/{sample}.GenerateBigWig.log",
	threads: 4  
	conda:
		"../envs/deeptools.yaml"
	shell:
		"bamCoverage -b {input.bam} -o {output} --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --exactScaling --binSize 1 --numberOfProcessors {threads} -v 2>{log}"