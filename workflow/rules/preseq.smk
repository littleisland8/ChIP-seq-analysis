rule Preseq:
	input:
		config["pipedir"] + "/" + "alignments/{sample}.bwa.filtered.rg.bam"
	output:
		"qc/{sample}.lc_extrap.txt"
	threads: 5
	log:
		"logs/{sample}.Preseq.log"
	conda:
		"../envs/preseq.yaml"
	shell:
		"preseq lc_extrap -bam -verbose -pe -output {output} {input} 2>{log}"