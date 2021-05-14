contaminated_samples = []
with open("data/contaminated_samples.txt", "r") as infile:
    for i, line in enumerate(infile):
        if i > 0:
            line = line.strip()
            contaminated_samples.append(line)

BATCHES = []
SAMPLES = []

with open("data/sequenced_isolates.txt", "r") as infile:
    for i, line in enumerate(infile):
        if i > 0:
            line = line.strip().split()
            if line[1] not in contaminated_samples:
                BATCHES.append(line[0])
                SAMPLES.append(line[1])


localrules:
    all,
    create_unitig_input,
    filter_significant,


rule all:
    input:
        "data/CLX/unitig_significance_annotated.txt",
        "data/CLX/unitigs_TCH1516_chromosome_position.txt",
        "data/CLX/unitig_annotated_TCH1516.txt",


def get_path(wildcards):
    if wildcards.batch == "batch_2":
        return "/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/{batch}/results/{sample}/LSARP_Results/Assembly/{sample}.genome.fa"
    else:
        return "/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/{batch}/final_results/{sample}/Assembly/{sample}.genome.fa"


rule annotation:
    input:
        fasta=get_path,
    params:
        name="{sample}",
        batch="{batch}",
    output:
        gff="data/annotations/{batch}/{sample}/{sample}.gff",
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt * 8000,
        time=lambda wildcards, attempt: attempt * 30,
    log:
        "logs/annotation/{batch}/{sample}.log",
    conda:
        "conda_envs/prokka.yml"
    shell:
        """
        mkdir -p annotations/{params.batch}
        prokka --force --outdir data/annotations/{params.batch}/{params.name} --prefix {params.name} --locustag {params.name} --genus Staphylococcus --species aureus --strain {params.name} --usegenus --cpus 8 {input.fasta}
        """


rule roary:
    input:
        expand(
            "data/annotations/{batch}/{sample}/{sample}.gff",
            zip,
            batch=BATCHES,
            sample=SAMPLES,
        ),
    output:
        "data/roary/roary/gene_presence_absence.csv",
        "data/roary/roary/core_gene_alignment.aln",
    params:
    resources:
        cpus=12,
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        time=lambda wildcards, attempt: attempt * 1200,
    log:
        "logs/roary.log",
    conda:
        "conda_envs/roary.yml"
    shell:
        """
        roary -p 12 -z -e -n -v -s -i 95 -f ./data/roary/roary {input}
        """


rule gubbins:
    input:
        "data/roary/roary/core_gene_alignment.aln",
    output:
        "data/gubbins/gubbins/core_alignment.final_tree.tre",
    params:
    resources:
        cpus=12,
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        time=lambda wildcards, attempt: attempt * 1200,
    log:
        "logs/gubbins.log",
    conda:
        "conda_envs/gubbins.yml"
    shell:
        """
        mkdir -p data/gubbins
        run_gubbins.py --threads {resources.cpus} --prefix data/gubbins/gubbins/core_alignment {input}
        """


rule create_unitig_input:
    input:
        "data/sequenced_isolates.txt",
        "data/contaminated_samples.txt",
    output:
        strain_list="data/unitigs/strain_list.txt"
    run:
        with open(output.strain_list, "w") as outfile:
            outfile.write("wgs_id\tpath\n")
            for sample, batch in zip(SAMPLES, BATCHES):
                if batch == "batch_2":
                    path = f"/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/{batch}/results/{sample}/LSARP_Results/Assembly/{sample}.genome.fa"
                else:
                    path = f"/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/{batch}/final_results/{sample}/Assembly/{sample}.genome.fa"
                outfile.write(f"{sample}\t{path}\n")


rule unitigs:
    input:
        strain_list="data/unitigs/strain_list.txt",
    output:
        "data/unitigs/s_aureus_unitigs/unitigs.txt",
    params:
        directory="data/unitigs/s_aureus_unitigs/",
    resources:
        cpus=4,
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        time=lambda wildcards, attempt: attempt * 480,
    log:
        "logs/unitig-counter.log",
    conda:
        "conda_envs/unitig-counter.yml"
    shell:
        """
        unitig-counter -strains {input.strain_list} -output unitigs/ -nb-cores {resources.cpus}
        mv unitigs/* {params.directory}
        rmdir unitigs/
        """


rule similarity_matrix:
    input:
        tree="data/gubbins/gubbins/core_alignment.final_tree.tre",
    output:
        matrix="data/gubbins/gubbins/similarity_matrix.txt",
    params:
    resources:
        cpus=1,
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        time=lambda wildcards, attempt: attempt * 30,
    log:
        "logs/similarity-matrix.log",
    conda:
        "conda_envs/pyseer.yml"
    shell:
        "python software/pyseer/scripts/phylogeny_distance.py --lmm {input.tree} > {output.matrix}"


rule lmm_gwas:
    input:
        pheno="data/{phenotype}/{phenotype}.txt",
        unitigs="data/unitigs/s_aureus_unitigs/unitigs.txt",
        similarity="data/gubbins/gubbins/similarity_matrix.txt",
    output:
        patterns="data/{phenotype}/unitig_patterns.txt",
        significance="data/{phenotype}/unitig_significance.txt",
    params:
    resources:
        cpus=1,
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        time=lambda wildcards, attempt: attempt * 360,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        pyseer --lmm --phenotypes {input.pheno} --similarity {input.similarity} --uncompressed --kmers {input.unitigs} --phenotype-column {wildcards.phenotype} --covariates {input.pheno} --use-covariates 2 --output-patterns {output.patterns} --cpu 1 > {output.significance}
        """


rule pyseer_count_patterns:
    input:
        "data/{phenotype}/unitig_patterns.txt",
    output:
        "data/{phenotype}/significance_limits.txt",
    params:
    resources:
        cpus=1,
        mem_mb=1000,
        time=10,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        python software/pyseer/scripts/count_patterns.py {input} > {output}
        """


rule filter_significant:
    input:
        script="scripts/filter_significant_unitigs.R",
        limit="data/{phenotype}/significance_limits.txt",
        unitig_significance="data/{phenotype}/unitig_significance.txt",
    output:
        "data/{phenotype}/unitig_significance_filtered.txt",
    shell:
        """
        module load R
        Rscript scripts/filter_significant_unitigs.R {input.limit} {input.unitig_significance} {output}
        """

rule annotate:
    input:
        unitig_filtered="data/{phenotype}/unitig_significance_filtered.txt",
        reference="data/reference/references.txt"
    output:
        "data/{phenotype}/unitig_significance_annotated.txt"
    resources:
        cpus=1,
        mem_mb=1000,
        time=10,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        annotate_hits_pyseer {input.unitig_filtered} {input.reference} {output}
        """

rule manhattan_input:
    input:
        unitigs="data/{phenotype}/unitig_significance.txt",
        reference_genome="data/reference/s_aureus_{reference}.fasta"
    output:
        "data/{phenotype}/unitigs_{reference}_position.txt"
    resources:
        cpus=1,
        mem_mb=1000,
        time=10
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        phandango_mapper {input.unitigs} {input.reference_genome} {output}
        """

rule annotate_onereference_allunitigs:
    input:
        significance="data/{phenotype}/unitig_significance.txt",
        reference_file="data/reference/references_{reference}.txt"
    output:
        "data/{phenotype}/unitig_annotated_{reference}.txt"
    resources:
        cpus=1,
        mem_mb=1000,
        time=10,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        annotate_hits_pyseer {input.significance} {input.reference_file} {output}
        """
