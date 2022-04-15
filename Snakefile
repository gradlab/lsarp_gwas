PHENOTYPE_FILE = "data/test/test.txt"
SAMPLES = []

with open(PHENOTYPE_FILE, "r") as infile:
    for i, line in enumerate(infile):
        if i > 0:
            line = line.strip().split()
            SAMPLES.append(line[0])


localrules:
    all,
    create_unitig_input,
    filter_significant,
    single_reference_file,
    manhattan_plot


rule all:
    input:
        "data/test/unitig_significance_annotated.txt",
        "data/test/unitigs_TCH1516_chromosome_position.txt",
        "data/test/manhattan_plot_TCH1516_chromosome.pdf"


def get_path(wildcards):
    return "/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/results/{sample}/Assembly/{sample}.genome.fa"


rule annotation:
    input:
        fasta=get_path,
    params:
        name="{sample}",
    output:
        gff="data/annotations/{sample}/{sample}.gff",
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt * 8000,
        time=lambda wildcards, attempt: attempt * 60,
    log:
        "logs/annotation/{sample}.log",
    conda:
        "conda_envs/prokka.yml"
    shell:
        """
        mkdir -p data/annotations/
        prokka --force --outdir data/annotations/{params.name} --prefix {params.name} --locustag {params.name} --genus Staphylococcus --species aureus --strain {params.name} --usegenus --cpus 8 {input.fasta}
        """


rule roary:
    input:
        expand(
            "data/annotations/{sample}/{sample}.gff",
            sample=SAMPLES,
        ),
    output:
        "data/roary/gene_presence_absence.csv",
        "data/roary/core_gene_alignment.aln",
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
        rm -rf data/roary
        roary -p 12 -z -e -n -v -s -i 95 -f ./data/roary {input}
        """


rule gubbins:
    input:
        "data/roary/core_gene_alignment.aln",
    output:
        "data/gubbins/core_alignment.final_tree.tre",
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
        run_gubbins.py --threads {resources.cpus} --prefix data/gubbins/core_alignment {input}
        """


rule create_unitig_input:
    input:
        "data/{phenotype}/{phenotype}.txt",
    output:
        strain_list="data/unitigs/strain_list_{phenotype}.txt"
    run:
        with open(output.strain_list, "w") as outfile:
            outfile.write("wgs_id\tpath\n")
            for sample in SAMPLES:
                path = f"/bulk/LSARP/genomics/pipeline/Staphylococcus_aureus/results/{sample}/Assembly/{sample}.genome.fa"
                outfile.write(f"{sample}\t{path}\n")


rule unitigs:
    input:
        strain_list="data/unitigs/strain_list_{phenotype}.txt",
    output:
        "data/unitigs/s_aureus_{phenotype}/unitigs.txt",
    params:
        directory="data/unitigs/s_aureus_{phenotype}/",
    resources:
        cpus=4,
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        time=lambda wildcards, attempt: attempt * 480,
    log:
        "logs/unitig-counter-{phenotype}.log",
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
        tree="data/gubbins/core_alignment.final_tree.tre",
    output:
        matrix="data/gubbins/similarity_matrix.txt",
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
        unitigs="data/unitigs/s_aureus_{phenotype}/unitigs.txt",
        similarity="data/gubbins/similarity_matrix.txt",
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
        pyseer --lmm --phenotypes {input.pheno} --similarity {input.similarity} --uncompressed --kmers {input.unitigs} --phenotype-column {wildcards.phenotype} --output-patterns {output.patterns} --cpu 1 > {output.significance}
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
        reference="reference/s_aureus/references.txt",
    output:
        annotated="data/{phenotype}/unitig_significance_annotated.txt",
        unannotated_fasta="data/{phenotype}/remaining_kmers.fa",
        unannotated_text="data/{phenotype}/remaining_kmers.txt"
    params:
        out_dir="data/{phenotype}/"
    shadow: "shallow"
    resources:
        cpus=1,
        mem_mb=1000,
        time=10,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        annotate_hits_pyseer {input.unitig_filtered} {input.reference} {output.annotated}
        mv remaining_kmers.fa {params.out_dir}
        mv remaining_kmers.txt {params.out_dir}
        """

rule phandango_input:
    input:
        unitigs="data/{phenotype}/unitig_significance.txt",
        reference_genome="reference/s_aureus/s_aureus_{reference}.fasta"
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

rule single_reference_file:
    input:
        reference_file="reference/s_aureus/references.txt"
    output:
        temp("reference/s_aureus/single_reference_{reference}_tmp.txt")
    params:
        reference="{reference}"
    shell:
        """
        grep {params.reference} {input.reference_file} > {output}
        """

rule annotate_onereference_allunitigs:
    input:
        significance="data/{phenotype}/unitig_significance.txt",
        reference_file="reference/s_aureus/single_reference_{reference}_tmp.txt"
    output:
        annotated_singleref="data/{phenotype}/unitig_annotated_{reference}.txt",
    shadow: "shallow"
    resources:
        cpus=1,
        mem_mb=1000,
        time=10,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        annotate_hits_pyseer {input.significance} {input.reference_file} {output.annotated_singleref}
        """

rule manhattan_plot:
    input:
        unitigs="data/{phenotype}/unitig_annotated_{reference}.txt",
        limit="data/{phenotype}/significance_limits.txt",
    output:
        plot="data/{phenotype}/manhattan_plot_{reference}.pdf"
    params:
        color="'#1f78b4'"
    shell:
        """
        module load R
        Rscript scripts/manhattan_plot.R {input.limit} {input.unitigs} {params.color} {output.plot}
        """
