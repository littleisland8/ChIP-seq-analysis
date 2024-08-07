rule fastp: ## aggiungere flag per i threads
	input:
		sample=[config["datadir"] + "/" + "{sample}_1.fq.gz", config["datadir"] + "/" + "{sample}_2.fq.gz"]
	output:
		trimmed=[config["pipedir"] + "/" + "data/{sample}_1.tr.fq.gz", config["pipedir"] + "/" + "data/{sample}_2.tr.fq.gz"],
		json=config["pipedir"] + "/" + "data/{sample}.json",
		failed=config["pipedir"] + "/" + "data/{sample}.failedreads.txt",
		html=config["pipedir"] + "/" + "data/{sample}.html",
		unpaired1=config["pipedir"] + "/" + "data/{sample}.u1.fq.gz",
		unpaired2=config["pipedir"] + "/" + "data/{sample}.u2.fq.gz"
	threads: 1
	log:
		"logs/{sample}.fastp.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra="--length_required 30 --detect_adapter_for_pe --disable_quality_filtering"
	wrapper:
		"v3.3.3/bio/fastp"