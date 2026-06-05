#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -J run_vtga
#SBATCH --output=slurm-%x.%j.out
#SBATCH --error=slurm-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your_email_address>

eval "$(conda shell.bash hook)"
conda activate snakemake

vtga.py \
    --reads_dir /path/to/raw_reads/directory \
    --genome_info /path/to/genome/information/table \
    -o output_dir 
