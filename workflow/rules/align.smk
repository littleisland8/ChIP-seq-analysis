rule BwaIndex:
    input:
        config["genome"]
    output:
        idx=multiext(config["genome"], ".bwt", ".pac", ".ann", ".amb", ".sa"),
    threads: 1
    log:
        "logs/BwaIndex.log",
    wrapper:
        "v3.13.1/bio/bwa/index"

rule BwaAlign:
	input:
		R1=config["pipedir"] + "/" + "data/{sample}_1.tr.fq.gz", 
		R2=config["pipedir"] + "/" + "data/{sample}_2.tr.fq.gz",
		ref=config["genome"]
	output:
		config["pipedir"] + "/" + "alignments/{sample}.bwa.srt.bam"
	log:
		"logs/{sample}.BwaAlign.log"
	message:
		"Mapping reads with bwa-mem2 {wildcards.sample}"
	threads: 5
	conda:
		"../envs/align.yaml"
	shell:
		"bwa mem -t {threads} {input.ref} {input.R1} {input.R2} |samtools sort -@ {threads} -o {output} --output-fmt BAM -T /scratch/tmp_sromagnoli/tmp5lxzdu46 2>{log}"

rule markDuplicates:
	input:
		bams=config["pipedir"] + "/" + "alignments/{sample}.bwa.srt.bam"
	output:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bai",
		metrics="qc/{sample}.MarkDuplicates.metrics.txt"
	threads: 5
	message:
		"Mark and remove PCR-optical duplicates"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
	resources:
		mem_mb=4096    
	log:
		"logs/{sample}.PicardRemoveDuplicates.log"
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule SamtoolsFilter:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bai"
	output:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.bam.bai"
	threads: 5
	conda:
		"../envs/filterbam.yaml"
	log:
		"logs/{sample}.SamtoolsFilter.log"
	params:
		script="workflow/scripts/bamtools.filter.json"
	shell:
		 "samtools view -F 0x004 -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -F 0x0400 -q 1 -b {input.bam} | bamtools filter -out {output.bam} -script {params.script} 2>{log} && samtools index -@ {threads} {output.bam}"

rule SamtoolsAddReadGroup:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.bam.bai"
	output:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam.bai"
	threads: 5
	conda:
		"../envs/samtools.yaml"
	log:
		"logs/{sample}.SamtoolsAddReadGroup.log"
	params:
		RG=r"'@RG\tID:{sample}\tSM:{sample}'"
	shell:
		 "samtools addreplacerg -r {params.RG} -@ {threads} -O BAM -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>{log}"

rule SamtoolsStats:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam.bai"
	output:
		"qc/{sample}.bwa.stats"
	threads: 5
	conda:
		"../envs/samtools.yaml"
	log:
		"logs/{sample}.SamtoolsStats.log"
	params:
		ref=config["genome"]
	shell:
		 "samtools stats -@ {threads} --reference {params.ref} {input.bam} > {output} 2>{log}"

rule SamtoolsFlagStats:
	input:
		metrics="qc/{sample}.MarkDuplicates.metrics.txt",
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam.bai"
	output:
		"qc/{sample}.bwa.flagstat"
	threads: 5
	conda:
		"../envs/samtools.yaml"
	log:
		"logs/{sample}.SamtoolsFlagStats.log"
	shell:
		 "samtools flagstat -@ {threads} {input.bam} > {output} 2>{log}"

rule SamtoolsStatsRaw:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bai"
	output:
		"qc/{sample}.bwa.raw.stats"
	threads: 5
	conda:
		"../envs/samtools.yaml"
	log:
		"logs/{sample}.SamtoolsStatsRaw.log"
	params:
		ref=config["genome"]
	shell:
		 "samtools stats -@ {threads} --reference {params.ref} {input.bam} > {output} 2>{log}"

rule SamtoolsFlagStatsRaw:
	input:
		metrics="qc/{sample}.MarkDuplicates.metrics.txt",
		bam=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.bwa.mdup.bai"
	output:
		"qc/{sample}.bwa.raw.flagstat"
	threads: 5
	conda:
		"../envs/samtools.yaml"
	log:
		"logs/{sample}.SamtoolsFlagStatsRaw.log"
	shell:
		 "samtools flagstat -@ {threads} {input.bam} > {output} 2>{log}"