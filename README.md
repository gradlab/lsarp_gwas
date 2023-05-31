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

### Core genome alignment and phylogeny
The core genome is identified using Roary v 3.13 [2], using annotations from Prokka [1], with a 95% identity cutoff and aligned with MAFFT [3]. To understand the population structure of our dataset, we generate a maximum likelihood phylogeny [4], which is based on the core genome alignment.
### Unitig presence matrix
To summarize the genomic variation in our sample, unitigs (unique sequences of variable length) are identified and counted using unitig-counter (https://github.com/johnlees/unitig-counter), which uses a compressed de Bruijn graph and is based on the approach in DBGWAS [5].
### Genome wide association study (GWAS) using linear mixed model (LMM)
To identify genetic variants associated with a phenotype of interest, we perform bacterial GWAS as implemented in pyseer v 1.3.6 [6]. We use LMM to control for population structure by including a similarity matrix generated from the phylogeny in the model. Unitig significance is determined with a likelihood ratio test and a Bonferroni-corrected p-value threshold based on the number of unique unitig presence/absence patterns.
### Unitig annotation and visualization 
Unitigs are annotated by mapping to a panel of closed, reference genomes representing each of the major clades or clonal complexes present in the dataset. Additionally, all unitigs are mapped to a single reference genome containing the gene content associated with the most significant unitigs for visualization (i.e. Manhattan plots). The presence of significant unitigs is displayed across the phylogeny using ITOL [7].
### References
1. 	Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014; 30:2068–2069. 
2. 	Page AJ, Cummins CA, Hunt M, et al. Roary: rapid large-scale prokaryote pan genome analysis. Bioinformatics 2015; 31:3691–3693. 
3. 	Katoh K, Standley DM. MAFFT: iterative refinement and additional methods. Methods Mol Biol Clifton NJ 2014; 1079:131–146. 
4. 	Croucher NJ, Page AJ, Connor TR, et al. Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. Nucleic Acids Res 2015; 43:e15–e15. 
5. 	Jaillard M, Lima L, Tournoud M, et al. A fast and agnostic method for bacterial genome-wide association studies: Bridging the gap between k-mers and genetic events. PLOS Genet 2018; 14:e1007758. 
6. 	Lees JA, Galardini M, Bentley SD, Weiser JN, Corander J. pyseer: a comprehensive tool for microbial pangenome-wide association studies. Bioinforma Oxf Engl 2018; 34:4310–4312. 
7. 	Letunic I, Bork P. Interactive Tree Of Life (iTOL) v4: recent updates and new developments. Nucleic Acids Res.

## Quick Start Guide
1. Install snakemake and activate environment
2. Clone this repository to your working directory
3. Edit `slurm/config.yaml`to reflect the path to your conda installation and ensure that `slurm/slurm-status.py` is excecutable
4. Make a directory called `software/`
5. Clone the pyseer repository to `software/` (`git clone https://github.com/mgalardini/pyseer.git`)
6. **Optional test:** skip directly to step 10 and run test files.
7. Make a directory within data for your specific phenotype `data/phenotype_name/`
8. Save or update input files in `data/phenotype_name/phenotype_name.txt`
9. Edit `Snakefile` to list desired output files under `rule all`. (e.g. replace 'test' in the example file names with your phenotype name)
10. Submit pipeline job with `sbatch start_snakemake.sh` (add email to `start_snakemake.sh` if you would like to receive email updates on your job)

## Software Requirements

This pipeline is implemented in snakemake. Additionally, there are some helper scripts that are not packaged with pyseer via conda, so the pyseer git repository should be cloned to a software/ directory. Other required software is installed via conda as the pipeline runs. 


### conda

Miniconda3 can be installed from https://docs.conda.io/en/latest/miniconda.html.

If you have having trouble with automatic software install via conda, try setting channel priority to flexible using the following command: `conda config --set channel_priority flexible`

### Snakemake

Snakemake can be installed via conda. Instructions here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

This pipeline uses a snakemake profile to interact with the slurm job submission system on ARC. The configuration file is in `slurm` in addition to a script that detects the job status from slurm. Jobs that fail due to requested time or memory limitations will be retried with the requested time/memory increased for a total of 3 tries (you can change by editing `restart-times` in config). 

This config file can be stored in `$HOME/.config/snakemake/slurm/config.yaml` if you would like to use it for other projects or in your working directory:

The config file should be edited to reflext your conda installation directory (see below for an example if miniconda3 was installed). Also check to see if `slurm-status.py` is executable. If not, use `chmod +x slurm-status.py` to change this.

```
restart-times: 3
cluster: "sbatch --parsable -t {resources.time} --nodes 1 --ntasks 1 --cpus-per-task {resources.cpus} --mem={resources.mem_mb} -o logs/slurm/{rule}_{wildcards}.out -e logs/slurm/{rule}_{wildcards}.err"
default-resources: [cpus=1, mem_mb=1000, time=60]
max-jobs-per-second: 1
max-status-checks-per-second: 10
local-cores: 1
use-conda: true
conda-prefix: /home/ARC_USERNAME/miniconda3/
jobs: 100
rerun-incomplete: true
cluster-status: "slurm-status.py"
```
### pyseer

The majority of pyseer components are installed via conda as part of the snakemake pipeline. However, the pipeline does require some helper scripts. To ensure the helper scripts can be found by the pipeline, run the following in the directory where you will be running the pipeline:

```
mkdir -p software/
cd software/
git clone https://github.com/mgalardini/pyseer.git
```

## Pipeline Components

### Snakefile
The Snakefile list rules describing required input and output for all files generated by the pipeline as well as the software needed to generate those files.
In general, the rules do not need to be edited. The exceptions are:
1. You must list desired output files for your phenotype under `rule all`. The files listed currently are examples. Please replace these file names with those listing your phenotype and chosen reference. If you don't want to choose a reference right away, list `[phenotype]/unitig_significance_annotated.txt`.
2. You can change the color of significant unitigs in the manhattan plot by editing the hex code listed under `params` in `rule manhattan_plot` (line 292 in the example Snakefile)

### start_snakemake.sh
Job submission script for running the pipeline (submit with `sbatch`).

### conda_envs/
Files describing conda environments for required pipeline software. 

### reference/
Directories for each species containing finished reference genomes and annotations to be used for
unitig mapping.

### scripts/

#### fastbaps.R
Example script for running fastbaps to assign BAPS groups using an alignment and phylogeny as input. 
(Not currently used by GWAS pipeline)

#### filter_significant_unitigs.R
Gets the appropriate Bonferroni corrected significance threshold from pyseer output and filters pyseer unitig output.

#### find_closest_reference_sig_unitigs.R
Finds the most common reference for the top 10% of significant unitigs. Note that this script will output the NCBI RefSeq accession number. This will be listed in the headers of the fasta. For example, if the reference is "NC_010079.1", the corresponding file is `s_aureus_TCH1516_chromosome.fasta`. You can easily figure this out by using `grep -H "NC_010079.1" *.fasta` in the `reference` directory. 

usage: `Rscript find_closest_reference_sig_unitigs.R unitigs_annotated`

example: `Rscript scripts/find_closest_reference_sig_unitigs.R  data/CLX/unitig_significance_annotated.txt`

#### get_mlst.py
Summarizes MLST output from pipeline results for all samples.
(Not currently used by GWAS pipeline)

#### manhattan_plot.R
Creates a Manhattan plot for unitigs mapped to a single reference.

usage: `Rscript manhattan_plot.R significance_limit unitigs_annotated_singleReference color outfile_name`

example: `Rscript scripts/manhattan_plot.R data/CLX/significance_limits.txt data/CLX/unitig_annotated_TCH1516_chromosome.txt "#D3D3D3" manhattan_plot_TCH1516_chromosome.pdf`

#### unitig_to_itol.py
Makes an annotation file that can be dragged and dropped into ITOL showing the presence and absence of a specific
unitig.

usage: `unitig_to_itol.py [-h] unitig_file unitig unitig_name color`

example: `scripts/unitig_to_itol.py data/unitigs/s_aureus_unitigs/unitigs.txt ACTGACTGACTGACTGACTG significant_unitig_1 "#D3D3D3"`

## Input Details

### phenotype file
The tab-delimited file with phenotypic information for the GWAS. Binary phenotypes should be represented by 0 or 1. This file should be stored in `data/phenotype1/phenotype1.txt` (with the actual phenotype name, not the word "phenotype1")

Format:

```
BI_NBR         phenotype1
BI_16_0013     1.36990875
BI_16_0017     0.80548875
BI_16_0028     1.11414875
```
### Directory Structure when you start

```
GWAS_directory/
    Snakefile
    start_snakemake.sh
    conda_envs/
    data/
        phenotype1/
            phenotype1.txt
    reference/
    scripts/
    software/
        pyseer/

```

## Output Details

### unitigs

Unitigs generated using unitig-caller (https://github.com/johnlees/unitig-counter)

### annotations

Annotations generated by Prokka (https://github.com/tseemann/prokka)

### roary

Pan-genome analysis by Roary (https://sanger-pathogens.github.io/Roary/)

### gubbins

Recombination detection and phylogeny by Gubbins (https://sanger-pathogens.github.io/gubbins/)

### gwas output (directory named by phenotype)

#### manhattan_plot_[reference].pdf
PDF of a Manhattan plot with unitigs mapped to a single reference. The color of significant points is dark blue by default, but that can be changed by editing the params
section of the manhattan_plot rule in the Snakefile.

#### remaining_kmers.fa
A fasta file of any significant unitigs that did not map to the panel of reference genomes.

#### remaining_kmers.txt
A text file of any significant unitigs that did not map to the panel of reference genomes.

#### significance_limits.txt
File describing appropriate Bonferroni corrected significance threshold for unitig p-values

#### unitig_patterns.txt
File used to calculate the number of unique unitig presence/absence patterns in the dataset (used for Bonferroni correction)

#### unitig_significance_filtered.txt
File containing information for significant unitigs only

#### unitig_significance.txt
File containing information for all tested unitigs

#### unitig_significance_annotated.txt
File containing information for significant unitigs annotated with a panel of reference genomes (mapped in order according to the file reference/[species]/references.txt)

#### unitigs_[reference]\_position.txt
File that can be used as an input for phandango (https://jameshadfield.github.io/phandango/#/) along with the .gff file for the reference chosen for an interactive
Manhattan plot in your browser.

#### unitigs_annotated_[reference].txt
All unitigs annotated according to a single reference

## Interpretation of Significant Unitigs

A major componenent of this pipeline is pyseer, which is a microbial GWAS tool using linear mixed models. A good first step in understanding GWAS results is to read through 
the documentation and tutorials for this software: https://pyseer.readthedocs.io/en/master/.

To visualize significant unitigs, you can use the interactive tools phandango and ITOL for an interactive Manhattan plot and annotated phylogeny, respectively. Additionally,
the pipeline can create a Manhattan plot figure.

### Manhattan plot
The pipeline will automatically generate a Manhattan plot for a single reference genome if `manhattan_plot_[reference].pdf` with `[reference]` replaced with the name of your chosen reference genome.

### phandango
For phandango, you will need the output file from the pipeline called unitigs_[reference]\_position.txt with [reference] replaced with the name of your chosen reference genome.
If you would like to use this file, make sure that the file name is listed under `rule all` in the Snakefile.

### ITOL
To make an annotation file for a specific unitig, you can use the `unitig_to_itol.py` script. This is not run automatically as part of the pipeline because there are 
potentially hundreds of significant unitigs to choose from. The output from this script can be dropped into ITOL (https://itol.embl.de/) along with the phylogeny (`data/gubbins/gubbins/core_alignment.final_tree.tre`) to view the distribution of the significant unitig across samples.

### Choosing a reference for annotation
A panel of reference genomes representing several lineages of _S. aureus_ is available in the `reference/s_aureus` directory. They are named using a standardized naming scheme, so if for example you would like to use `s_aureus_MW2_chromosome.fasta` as the reference, the pipeline will recognize output file names with `MW2_chromosome` in the name. The test run uses the TCH1516_chromosome as the reference, which is a USA300 strain.

If you aren't sure which single reference genome you want to use, you can take the output of significant unitigs mapped to the whole panel (`unitig_significance_annotated.txt`) and check which reference was chosen for the most significant unitigs using `find_closest_reference_sig_unitigs.R`.

## Alternate versions of the pipeline

If you already have a phylogeny for the samples in your dataset, you can skip a few steps in the pipeline. To use this version of the pipeline, copy the Snakefile in the `alternate_snakefile` directory to the main directory and edit the variable "TREE_FILE" at the top.
