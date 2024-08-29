rule MergeBed10InputNarrow:
	input:
		expand(f"peaks/narrow/10_input/{{sample}}_peaks.narrowPeak", sample={**config["controls_macs2"],**config["lactate_macs2"]}.values())
	output:
		"results/narrow/10_input/total_peaks.narrowPeak"
	threads: 1
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/MergeBed10InputNarrow.log"
	params:
		chrs="resources/classic.chrs.txt"
	shell:
		"cat {input} |cut -f 1,2,3 |sortBed |mergeBed -i stdin | grep -w -f {params.chrs} |sortBed -i stdin > {output} 2>{log}"

rule multiCovBams10InputNarrow:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["Ab"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai",sample=config["Ab"].values()),
		bed="results/narrow/10_input/total_peaks.narrowPeak"
	output:
		"results/narrow/10_input/total.multicov.bed"
	threads: 1 
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/multiCovBams10InputNarrow.log"
	shell:
		"multiBamCov -bams {input.bam} -bed {input.bed} > {output} 2>{log}"

rule multiCovBams10InputParsedNarrow:
	input:
		bed="results/narrow/10_input/total.multicov.bed"
	output:
		"results/narrow/10_input/total.multicov.parsed.bed"
	threads: 1 
	log:
		"logs/multiCovBams10InputParsedNarrow.log"
	params:
		sample=expand("{sample}",sample=config["Ab"].values())
	shell:
		'''
		cat <(echo -e "chrom start end" {params.sample}|tr " " "\t") {input} > {output} 2>{log}
		'''

rule MergeBedLAvsCtrNarrow:
	input:
		expand(f"peaks/narrow/LAvsCtr/{{sample}}_peaks.narrowPeak", sample=config["lactate_macs2"].values())
	output:
		"results/narrow/LAvsCtr/total_peaks.narrowPeak"
	threads: 1
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/MergeBedLAvsCtrNarrow.log"
	params:
		chrs="resources/classic.chrs.txt"
	shell:
		"cat {input} |cut -f 1,2,3 |sortBed |mergeBed |grep -w -f {params.chrs} |sortBed -i stdin > {output} 2>{log}"

rule multiCovBamsLAvsCtrNarrow:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["Ab"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai",sample=config["Ab"].values()),
		bed="results/narrow/LAvsCtr/total_peaks.narrowPeak"
	output:
		"results/narrow/LAvsCtr/total.multicov.bed"
	threads: 1 
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/multiCovBamsLAvsCtrNarrow.log"
	shell:
		"multiBamCov -bams {input.bam} -bed {input.bed} > {output} 2>{log}"

rule multiCovBamsLAvsCtrParsedNarrow:
	input:
		bed="results/narrow/LAvsCtr/total.multicov.bed"
	output:
		"results/narrow/LAvsCtr/total.multicov.parsed.bed"
	threads: 1 
	log:
		"logs/multiCovBamsLAvsCtrParsedNarrow.log"
	params:
		sample=expand("{sample}",sample=config["Ab"].values()),
		chrs="resources/classic.chrs.txt"
	shell:
		'''
		cat <(echo -e "chrom start end" {params.sample}|tr " " "\t") {input} > {output} 2>{log}
		'''

rule MergeBed10InputBroad:
	input:
		expand(f"peaks/broad/10_input/{{sample}}_peaks.broadPeak", sample={**config["controls_macs2"],**config["lactate_macs2"]}.values())
	output:
		"results/broad/10_input/total_peaks.broadPeak"
	threads: 1
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/MergeBed10InputBroad.log"
	params:
		chrs="resources/classic.chrs.txt"
	shell:
		"cat {input} |cut -f 1,2,3 |sortBed |mergeBed -i stdin | grep -w -f {params.chrs} |sortBed -i stdin > {output} 2>{log}"

rule multiCovBams10InputBroad:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["Ab"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai",sample=config["Ab"].values()),
		bed="results/broad/10_input/total_peaks.broadPeak"
	output:
		"results/broad/10_input/total.multicov.bed"
	threads: 1 
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/multiCovBams10InputBroad.log"
	shell:
		"multiBamCov -bams {input.bam} -bed {input.bed} > {output} 2>{log}"

rule multiCovBams10InputParsedBroad:
	input:
		bed="results/broad/10_input/total.multicov.bed"
	output:
		"results/broad/10_input/total.multicov.parsed.bed"
	threads: 1 
	log:
		"logs/multiCovBams10InputParsedBroad.log"
	params:
		sample=expand("{sample}",sample=config["Ab"].values())
	shell:
		'''
		cat <(echo -e "chrom start end" {params.sample} |tr " " "\t") {input} > {output} 2>{log}
		'''

rule MergeBedLAvsCtrBroad:
	input:
		expand(f"peaks/broad/LAvsCtr/{{sample}}_peaks.broadPeak", sample=config["lactate_macs2"].values())
	output:
		"results/broad/LAvsCtr/total_peaks.broadPeak"
	threads: 1
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/MergeBedLAvsCtrBroad.log"
	params:
		chrs="resources/classic.chrs.txt"
	shell:
		"cat {input} |cut -f 1,2,3 |sortBed |mergeBed |grep -w -f {params.chrs} |sortBed -i stdin > {output} 2>{log}"

rule multiCovBamsLAvsCtrBroad:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["Ab"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai",sample=config["Ab"].values()),
		bed="results/broad/LAvsCtr/total_peaks.broadPeak"
	output:
		"results/broad/LAvsCtr/total.multicov.bed"
	threads: 1 
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/multiCovBamsLAvsCtrBroad.log"
	shell:
		"multiBamCov -bams {input.bam} -bed {input.bed} > {output} 2>{log}"

rule multiCovBamsLAvsCtrParsedBroad:
	input:
		bed="results/broad/LAvsCtr/total.multicov.bed"
	output:
		"results/broad/LAvsCtr/total.multicov.parsed.bed"
	threads: 1 
	log:
		"logs/multiCovBamsLAvsCtrParsedBroad.log"
	params:
		sample=expand("{sample}",sample=config["Ab"].values()),
		chrs="resources/classic.chrs.txt"
	shell:
		'''
		cat <(echo -e "chrom start end" {params.sample}|tr " " "\t") {input} > {output} 2>{log}
		'''