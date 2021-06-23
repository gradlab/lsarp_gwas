#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH -p cpu2019
#SBATCH -t 3-00:00
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err
#SBATCH --mail-type=END
#SBATCH --mail-user=

mkdir -p logs/slurm
snakemake --profile slurm --rerun-incomplete --latency-wait 90
