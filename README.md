# lsarp_gwas

Description of GWAS pipeline developed for LSARP project. Pipeline scripts (but not data) are additionally tracked via git and backed up to a GitHub repository.  

## Pipeline Overview

This pipeline generates a table of unitigs (unique sequences that are present in > 1% of < 99% of the input genomes) that are significantly associated with a phenotype after controlling for population structure using a linear mixed model (LMM) implemented in the software pyseer. 

The general pipeline is as follows:

```
de novo assemblies -> annotations -> core genome alignment -> phylogeny -> similarity matrix
        |
        V
      unitigs

similiarity matrix + unitigs + phenotype -> LMM -> significant unitigs -> annotated unitigs -> phylogeny annotation
                                                                                   |
                                                                                   V        
                                                                            Manhattan plot

```
Details of the software used can be found by reading the Snakefile.

## Software Requirements

This pipeline is implemented in snakemake. Other required software is installed via conda as the pipeline runs. 

### conda

Links to conda installation

### Snakemake

Links to installation

Example slurm config

### pyseer

Explanation about pyseer accessory scripts and getting them from github

## Pipeline Components

### Snakefile

### start_snakemake.sh

### conda_envs/

### scripts/

## Input Details

### sequenced_isolates.txt

### contaminanted_isolates.txt

### phenotype file

## Output Details

### unitigs

### annotations

### roary

### gubbins

### gwas output (directory named by phenotype)

## Interpretation of Significant Unitigs

