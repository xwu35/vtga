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
    --long_reads /path/to/nanopore/long_reads \
    --short_r1 /path/to/illumina/forward_reads \
    --short_r2 /path/to/illumina/reverse_reads \
    -o output_dir 
