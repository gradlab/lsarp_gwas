md -p ~/.config/snakemake/slurm

sed 's/ARC_USERNAME/'$USER'/g' config/slurm/config.yaml > ~/.config/snakemake/slurm/config.yaml

cd ~

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

chmod +x Miniconda3-latest-Linux-x86_64.sh
    
./Miniconda3-latest-Linux-x86_64.sh
    
conda create -c conda-forge -c bioconda -n gwas pyseer mamba snakemake-minimal pip ipykernel

    