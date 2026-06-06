#!/usr/bin/env python3

import sys
import os
import pandas as pd
import subprocess
import click
from datetime import datetime

version = "1.1.0"
@click.version_option(version, "--version", "-v")

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n vtga.py --reads_dir <reads directory> '
    '--genome_info <genome information table> -o <output directory>'
)
@click.option(
    '--reads_dir',
    required=True,
    type=click.Path(dir_okay=True, exists=True, resolve_path=True),
    help='Reads directory'
)
@click.option(
    '--genome_info',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Genome information table (tab separated).' 
    ' The input table for hybrid-read assemblies must contain four columns: genome, short_R1, short_R2, long_reads.'
    ' For short-read only assemblies, the long column is not required.'
    ' For long-read only assemblies, the short_R1 and short_R2 columns are not required.'
)
@click.option("-o", 
    '--output_dir',
    default="OUTPUT",
    type=click.Path(dir_okay=True, resolve_path=True),
    show_default=True,
    help=('Output directory')
)
@click.option(
    '--reads_type',
    default='hybrid',
    type=str,
    show_default=True,
    help=('Reads type; available options are: short, long, hybrid')
)
@click.option(
    '--count',
    default=4,
    type=int,
    show_default=True,
    help=('Number of subgenomed read sets. This option only applies when long reads are provided')
)
@click.option(
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce'
)
@click.option(
    '--conda_envs',
    default='',
    show_default=True,
    help='Directory to store conda environments.'
    ' By default, the "conda_env" directory within vtga is used'
)
@click.option(
    '--profile',
    default='slurm',
    show_default=True,
    help='Snakemake profile for cluster execution'
)

def run_genomeassembly(reads_dir, genome_info, 
          output_dir, reads_type, count, dryrun, conda_envs, profile):

    # get snakefile and conda_envs path
    script_dir=os.path.dirname(os.path.abspath(__file__))
    snakefile=os.path.join(script_dir, "workflow", "Snakefile")
    default_envs=os.path.join(script_dir, "conda_envs")

    df = pd.read_table(genome_info)

    if reads_type=="hybrid":
        required_columns = {"genome", "short_R1", "short_R2", "long_reads"}
        missing_cols = required_columns - set(df.columns)
        if missing_cols:
            raise ValueError(f"Genome information table is missing required columns: {missing_cols}")

    if reads_type=="short":
        required_columns = {"genome", "short_R1", "short_R2"}
        missing_cols = required_columns - set(df.columns)
        if missing_cols:
            raise ValueError(f"Sample information table is missing required columns: {missing_cols}")
    
    if reads_type=="long":
        required_columns = {"genome", "long_reads"}
        missing_cols = required_columns - set(df.columns)
        if missing_cols:
            raise ValueError(f"Sample information table is missing required columns: {missing_cols}")
        
    # write run log if it is not a dry run
    if not dryrun:
        os.makedirs(output_dir, exist_ok=True)
        logfile = os.path.join(output_dir, f"{os.path.basename(output_dir)}_run.log")
        with open(logfile, "w") as log:
            log.write("================VTGenomeAssembly run log==============\n")
            log.write(f"Start time: {datetime.now()}\n")
            log.write(f"VTGenomeAssembly version: {version}\n")
            log.write(f"Reads directory: {reads_dir}\n")
            log.write(f"Genome information table: {genome_info}")

    cmd = (
        'snakemake --snakefile {snakefile} '
        '--use-conda --conda-frontend mamba '
        '--conda-prefix {envs} '
        '--profile {profile} --rerun-incomplete ' 
        '--printshellcmds --nolock --show-failed-logs '
        '{dryrun} '
        '--config reads_dir={reads_dir} genome_info={genome} '
        'results_dir={results} reads_type={type} count={c} conda_envs={envs}'
        ).format(
            snakefile=snakefile,
            envs=default_envs if conda_envs=='' else conda_envs,
            profile=profile,
            dryrun='--dryrun' if dryrun else '',
            reads_dir=reads_dir,
            genome=genome_info,
            results=output_dir,
            type=reads_type,
            c=count
            )

    # run snakemake with command-line config
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError:
        print("Snakemake failed. see log for details.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    run_genomeassembly()