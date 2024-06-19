rule CreateSequenceDictionary:
	input:
		config["genome"]
	output:
		config["genome"].split(".fa")[0] + ".dict"
	log:
		"logs/CreateSequenceDictionary.log"
	threads: 1
	conda:
		"../envs/gatk4.yaml"
	shell:
		"java -jar picard.jar CreateSequenceDictionary -R {input} -O {output}"

rule bwa_align:
	input:
		reads=["data/{sample}.R1.tr.fastq.gz", "data{sample}.R2.tr.fastq.gz"],
		idx=multiext(config['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
	output:
		temp("alignments/{sample}.bwa.srt.bam")
	log:
		"logs/{sample}.bwa_align.log"
	message:
		"Mapping reads with bwa-mem2 on tumor {wildcards.sample}"
	benchmark:
		"SnakeWES/benchmarks/{sample}.bwaTumor.txt"
	params:
		#extra=r"-R '@RG\tID:{sample}\tSM:{sample}_tumor'",
		sort="samtools",  
		sort_order="coordinate"
	threads: 1
	wrapper:
		"v3.3.3/bio/bwa-mem2/mem"

rule PicarRemoveDuplicates:
	input:
		bams="alignments/{sample}.bwa.srt.bam"
	output:
		bam=temp("alignments/{sample}.bwa.dd.bam"),
		bai=temp("alignments/{sample}.bwa.dd.bai"),
		metrics="data/{sample}.bwa.metrics.picard.txt"
	threads: 1
	log:
		"logs/{sample}.PicarRemoveDuplicates.log"
	message:
		"Mark and remove PCR-optical duplicates - Picard"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
	resources:
		mem_mb=4096    
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule SamtoolsRemoveDuplicates:
	input:
		bams="alignments/{sample}.bwa.dd.bam",
		bai="alignments/{sample}.bwa.dd.bai"
	output:
		bam=temp("alignments/{sample}.bwa.dd.samtools.bam"),
		bai=temp("alignments/{sample}.bwa.dd.samtools.bam.bai"),
		metrics="data/{sample}.bwa.metrics.samtools.txt"
	threads: 1
	envs:
		"../envs/samtools.yaml"
	log:
		"logs/{sample}.SamtoolsRemoveDuplicates.log"
	message:
		"Mark and remove PCR-optical duplicates - Samtools"
	shell:
		 "samtools markdup -@ {threads} -r -f {input.metrics} {input.bams} {output.bam} && samtools index -@ {threads} {output.bam}"

rule SamtoolsFilter:
	input:
		bam="alignments/{sample}.bwa.dd.samtools.bam",
		bai="alignments/{sample}.bwa.dd.samtools.bam.bai"
	output:
		bam="alignments/{sample}.bwa.filtered.bam",
		bai="alignments/{sample}.bwa.filtered.bam.bai"
	threads: 1
	envs:
		"../envs/filterbam.yaml"
	log:
		"logs/{sample}.SamtoolsFilter.log"
	params:
		script="workflow/scripts/bamtools.filter.json"
	shell:
		 "samtools view -F 0x004 -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -F 0x0400 -b {input.bam} | bamtools filter -out {output} -script {params.script}"