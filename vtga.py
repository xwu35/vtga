#!/usr/bin/env python3

import sys
import os
import subprocess
import click
from datetime import datetime

version = "1.0.0"
@click.version_option(version, "--version", "-v")

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n vtga.py --long_reads <ONT reads> '
    '--short_r1 <Illumina R1> --short_r2 <Illumina R2> -o <output directory>'
)
@click.option(
    '--long_reads',
    required=False,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Long-read FASTQ'
)
@click.option(
    '--short_r1',
    required=False,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Short-read R1 FASTQ'
)
@click.option(
    '--short_r2',
    required=False,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Short-read R2 FASTQ'
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
    help=('Number of subsampled read sets. This option only applies when long reads are provided')
)
@click.option(
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce'
)
@click.option(
    '--profile',
    default='slurm',
    show_default=True,
    help='Snakemake profile for cluster execution'
)
def run_genomeassembly(long_reads, short_r1, short_r2, 
          output_dir, reads_type, count, dryrun, profile):

    # get snakefile and conda_envs path
    script_dir=os.path.dirname(os.path.abspath(__file__))
    snakefile=os.path.join(script_dir, "workflow", "Snakefile")
    conda_envs=os.path.join(script_dir, "conda_envs")

    if reads_type=="hybrid":
        if long_reads is None:
            raise click.UsageError("--long_reads is required when both short and long reads are used.")
        if short_r1 is None:
            raise click.UsageError("--short_r1 is required when both short and long reads are used.")
        if short_r2 is None:
            raise click.UsageError("--short_r2 is required when both short and long reads are used.")
        
    if reads_type=="short":
        if short_r1 is None:
            raise click.UsageError("--short_r1 is required when short reads assembly is selected.")
        if short_r2 is None:
            raise click.UsageError("--short_r2 is required when short reads assembly is selected.")
    
    if reads_type=="long":
        if long_reads is None:
            raise click.UsageError("--long_reads is required when long reads assembly is selected.")
        
    # write run log if it is not a dry run
    if not dryrun:
        os.makedirs(output_dir, exist_ok=True)
        logfile = os.path.join(output_dir, f"{os.path.basename(output_dir)}_run.log")
        with open(logfile, "w") as log:
            log.write("================VTGenomeAssembly run log==============\n")
            log.write(f"Start time: {datetime.now()}\n")
            log.write(f"VTGenomeAssembly version: {version}\n")
            log.write(f"Long reads: {long_reads}\n")
            log.write(f"Short reads: {short_r1}, {short_r2}")

    cmd = (
        'snakemake --snakefile {snakefile} '
        '--use-conda --conda-frontend mamba '
        '--conda-prefix {conda_envs} '
        '--profile {profile} --rerun-incomplete ' 
        '--printshellcmds --nolock --show-failed-logs '
        '{dryrun} '
        '--config long_reads={long} '
        'short_r1={r1} short_r2={r2} '
        'results_dir={results} reads_type={type} count={c}'
        ).format(
            snakefile=snakefile,
            conda_envs=conda_envs,
            profile=profile,
            dryrun='--dryrun' if dryrun else '',
            long=long_reads,
            r1=short_r1,
            r2=short_r2,
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
