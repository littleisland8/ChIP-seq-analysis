rule Macs2NarrowpeakLAIgG:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/narrow/IgG/{{sample}}_peaks.narrowPeak", sample=config["lactate_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["lactate_IgG"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --outdir peaks/narrow/IgG 2>logs/${{final}}.Macs2NarrowpeakLAIgG.log ; done
		"""

rule Macs2NarrowpeakLA10:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/narrow/10_input/{{sample}}_peaks.narrowPeak", sample=config["lactate_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["lactate_10"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --outdir peaks/narrow/10_input 2>logs/${{final}}.Macs2NarrowpeakLA10.log ; done
		"""

rule Macs2NarrowpeakCtrIgG:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/narrow/IgG/{{sample}}_peaks.narrowPeak", sample=config["controls_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["controls_IgG"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --outdir peaks/narrow/IgG 2>logs/${{final}}.Macs2NarrowpeakCtrIgG.log ; done
		"""

rule Macs2NarrowpeakCtr10:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/narrow/10_input/{{sample}}_peaks.narrowPeak", sample=config["controls_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["controls_10"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --outdir peaks/narrow/10_input 2>logs/${{final}}.Macs2NarrowpeakCtr10.log ; done
		"""

rule Macs2BroadpeakLAIgG:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/broad/IgG/{{sample}}_peaks.broadPeak", sample=config["lactate_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["lactate_IgG"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --broad --outdir peaks/broad/IgG 2>logs/${{final}}.Macs2BroadpeakLAIgG.log ; done
		"""

rule Macs2BroadpeakLA10:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/broad/10_input/{{sample}}_peaks.broadPeak", sample=config["lactate_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["lactate_10"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --broad --outdir peaks/broad/10_input 2>logs/${{final}}.Macs2BroadpeakLA10.log ; done
		"""

rule Macs2BroadpeakCtrIgG:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/broad/IgG/{{sample}}_peaks.broadPeak", sample=config["controls_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["controls_IgG"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --broad --outdir peaks/broad/IgG 2>logs/${{final}}.Macs2BroadpeakCtrIgG.log ; done
		"""

rule Macs2BroadpeakCtr10:
	input:
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam", sample=config["samples"].values()),
		bai=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.bwa.filtered.rg.bam.bai", sample=config["samples"].values())
	output:
		expand(f"peaks/broad/10_input/{{sample}}_peaks.broadPeak", sample=config["controls_macs2"].values())
	threads: 1
	conda: 
		"../envs/macs2.yaml"
	params:
		paired=' '.join(list(config["controls_10"].values())),
		aligndir=config["pipedir"] + "/" + "alignments"
	shell:
		"""
		for sample in {params.paired}; do t=$(echo ${{sample}} |cut -d "." -f 1) && c=$(echo ${{sample}} |cut -d "." -f 2) && final=$(echo ${{sample}} |cut -d "." -f1 |sed -e "s/[0-9]//") && macs2 callpeak -t {params.aligndir}/${{t}}.bwa.filtered.bam -c {params.aligndir}/${{c}}.bwa.filtered.bam -f BAM -g hs -n ${{final}} --broad --outdir peaks/broad/10_input 2>logs/${{final}}.Macs2BroadpeakCtr10.log ; done
		"""