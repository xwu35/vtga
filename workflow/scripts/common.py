#!/usr/bin/env python

import os
import pandas as pd

# function to get variables and sequences path
def parse_genomes_and_sequences(reads_type, metadata, raw_reads_dir=""):
    
    df = pd.read_table(metadata).set_index("genome", drop=False)
    genomes = df.index.tolist()

    if reads_type=="hybrid":
      short_R1_map = {}
      short_R2_map = {}
      long_reads_map = {}
      for genome in genomes:
        genome_R1 = df.loc[genome, "short_R1"]
        genome_R2 = df.loc[genome, "short_R2"]
        genome_long = df.loc[genome, "long_reads"]
        short_R1_map[genome] = os.path.join(raw_reads_dir, genome_R1)
        short_R2_map[genome] = os.path.join(raw_reads_dir, genome_R2)
        long_reads_map[genome] = os.path.join(raw_reads_dir, genome_long)
      return genomes, short_R1_map, short_R2_map, long_reads_map

    if reads_type=="short":
      short_R1_map = {}
      short_R2_map = {}
      for genome in genomes:
        genome_R1 = df.loc[genome, "short_R1"]
        genome_R2 = df.loc[genome, "short_R2"]
        short_R1_map[genome] = os.path.join(raw_reads_dir, genome_R1)
        short_R2_map[genome] = os.path.join(raw_reads_dir, genome_R2)
      return genomes, short_R1_map, short_R2_map
    
    if reads_type=="long":
      long_reads_map = {}
      for genome in genomes:
        genome_long = df.loc[genome, "long_reads"]
        long_reads_map[genome] = os.path.join(raw_reads_dir, genome_long)
      return genomes, long_reads_map