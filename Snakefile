contaminated_samples = []
with open("data/contaminated_samples.txt", "r") as infile:
    for i,line in enumerate(infile):
        if i > 0:
            line = line.strip()
            contaminated_samples.append(line)

BATCHES = []
SAMPLES = []

with open("data/sequenced_isolates.txt", "r") as infile:
    for i,line in enumerate(infile):
        if i > 0:
            line = line.strip().split()
            if line[1] not in contaminated_samples:
                BATCHES.append(line[0])
                SAMPLES.append(line[1])

localrules: all

rule all:
    input:
        expand("data/annotations/{batch}/{sample}/{sample}.gff", zip, batch = BATCHES, sample = SAMPLES),
        "data/gubbins/core_alignment.final_tree.tre"

def get_path(wildcards):
    if wildcards.batch == "batch_2":
        return "/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/{batch}/results/{sample}/LSARP_Results/Assembly/{sample}.genome.fa"
    else:
        return "/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/{batch}/final_results/{sample}/Assembly/{sample}.genome.fa"
    
 
rule annotation:
    input:
        fasta=get_path
    params:
        name="{sample}",
        batch="{batch}",
        partition="cpu2019"
    output:
        gff="data/annotations/{batch}/{sample}/{sample}.gff"
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt * 8000,
        time=lambda wildcards, attempt: attempt * 30
    log:
        "logs/annotation/{batch}/{sample}.log"
    shell:
        """
        mkdir -p annotations/{params.batch}
        prokka --force --outdir data/annotations/{params.batch}/{params.name} --prefix {params.name} --locustag {params.name} --genus Staphylococcus --species aureus --strain {params.name} --usegenus --cpus 8 {input.fasta}
        """

rule roary:
    input:
        expand("data/annotations/{batch}/{sample}/{sample}.gff", zip, batch = BATCHES, sample = SAMPLES)
    output:
        "data/roary/gene_presence_absence.csv",
        "data/roary/core_gene_alignment.aln"
    params:
        partition="cpu2019"
    resources:
        cpus=12,
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        time=lambda wildcards, attempt: attempt * 1200
    log:
        "logs/roary.log"
    shell:
        """
        roary -p 12 -z -e -n -v -s -i 95 -f ./data/roary {input}
        """

rule gubbins:
    input:
        "data/roary/210207-tm__roary/core_gene_alignment.aln" 
    output:
        "data/gubbins/core_alignment.final_tree.tre"
    params:
        partition="cpu2019"
    resources:
        cpus=12,
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        time=lambda wildcards, attempt: attempt * 1200
    log:
        "logs/gubbins.log"
    shell:
        """
        mkdir -p data/gubbins
        run_gubbins.py --threads {resources.cpus} --prefix data/gubbins/core_alignment {input}
        """

rule unitigs:

rule lmm:

rule
