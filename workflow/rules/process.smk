rule trim_short_reads:
    input: 
        r1=R1,
        r2=R2
    output:
        r1_filt=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_R1_trimmed.fastq.gz"),
        r2_filt=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_R2_trimmed.fastq.gz"),
        html=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_fastp.html"),
        json=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_fastp.json")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "fastp.yml")
    shell:
        """
        fastp --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1_filt} \
            --out2 {output.r2_filt} \
            --unpaired1 /dev/null \
            --unpaired2 /dev/null \
            -h {output.html} \
            -j {output.json} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --trim_poly_g # By default trimming is automatically enabled for Illumina NextSeq/NovaSeq data, but it did not, so force polyG tail trimming. 
        """

rule shovill_short_reads_assembly:
    input: 
        r1=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_R1_trimmed.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_R2_trimmed.fastq.gz")
    output:
        genome=os.path.join(RESULTS_DIR, "short_reads_only", "genome_sequences", "contigs.fa")
    params:
        dir=os.path.join(RESULTS_DIR, "short_reads_only", "genome_sequences"),
        gsize=config["shovill"]["genome_size"],
        assembler=config["shovill"]["assembler"],
        length=config["shovill"]["length"],
        coverage=config["shovill"]["coverage"]
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "shovill.yml")
    shell:
        """
        shovill --outdir {params.dir} \
            --R1 {input.r1} \
            --R2 {input.r2} \
            --cpus {threads} \
            --assembler {params.assembler} \
            --gsize {params.gsize} \
            --minlen {params.length} \
            --mincov {params.coverage} \
            --force 
        """

rule checkm2_shovill_assembled_genome:
    input: 
        done=os.path.join(dir["db"], "CheckM2_database", ".checkm2_db_done"),
        genome=os.path.join(RESULTS_DIR, "short_reads_only", "genome_sequences", "contigs.fa")
    output:
        report=os.path.join(RESULTS_DIR, "short_reads_only", "genome_quality", "quality_report.tsv")
    params:
        dir=os.path.join(RESULTS_DIR, "short_reads_only", "genome_quality"),
        db=os.path.join(dir["db"], "CheckM2_database", "uniref100.KO.1.dmnd")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "checkm2.yml")
    shell:
        """
        checkm2 predict --threads {threads} \
            --input {input.genome} \
            --output-directory {params.dir} \
            --database_path {params.db} \
            --force 
        """

rule trim_long_reads:
    input: 
        long=LONG
    output:
        long_filt=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "long_reads_trimmed.fastq.gz")
    params:
        length=config["filtlong"]["length"],
        percent=config["filtlong"]["percent"]
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "filtlong.yml")
    shell:
        """
        filtlong --min_length {params.length} --keep_percent {params.percent} {input.long} | gzip > {output.long_filt}
        """

rule calculate_genome_size:
    input:
        os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "long_reads_trimmed.fastq.gz")
    output:
        size=os.path.join(RESULTS_DIR, "intermediate", "estimated_genome_size.txt")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "lrge.yml")
    shell:
        """
        lrge -t {threads} {input} > {output.size}
        """

rule autocycler_subsample:
    """
    Step 1: subsample the long-read set into multiple files
    """
    input: 
        long_filt=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "long_reads_trimmed.fastq.gz"),
        genome_size=os.path.join(RESULTS_DIR, "intermediate", "estimated_genome_size.txt")
    output:
        outdir=directory(os.path.join(RESULTS_DIR, "intermediate", "subsampled_reads")),
        sub_seq=temp(reads_list),
        version=os.path.join(RESULTS_DIR, "long_reads_only", "software_versions.txt")
    params:
        num_sub=SUBREADS_COUNT
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        genome_size=$(cat {input.genome_size}) 

        autocycler subsample --reads {input.long_filt} --out_dir {output.outdir} --genome_size "$genome_size" --count {params.num_sub}

        autocycler --version >> {output.version}
        flye -v | sed 's/^/flye /' >> {output.version}
        plassembler --version >> {output.version}
        """

rule autocycler_assembly_flye:
    """
    Step 2: assemble each subsampled file using flye
    """
    input: 
        sub_seq=os.path.join(RESULTS_DIR, "intermediate", "subsampled_reads", "sample_{number}.fastq"),
        genome_size=os.path.join(RESULTS_DIR, "intermediate", "estimated_genome_size.txt")
    output:
        multiext(os.path.join(RESULTS_DIR, "intermediate", "assemblies", "flye_{number}"), ".fasta", ".gfa", ".log")
    params:
        prefix=os.path.join(RESULTS_DIR, "intermediate", "assemblies", "flye_{number}")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        genome_size=$(cat {input.genome_size})

        autocycler helper flye --reads {input.sub_seq} --out_prefix {params.prefix} --threads {threads} --genome_size "$genome_size"
        """

rule autocycler_assembly_plassembler:
    """
    Step 2: assemble each subsampled file using plassembler (for plasmids)
    """
    input: 
        done=os.path.join(dir["conda_envs"], ".plassembler_db_done"),
        sub_seq=os.path.join(RESULTS_DIR, "intermediate", "subsampled_reads", "sample_{number}.fastq"),
        genome_size=os.path.join(RESULTS_DIR, "intermediate", "estimated_genome_size.txt")
    output:
        multiext(os.path.join(RESULTS_DIR, "intermediate", "assemblies", "plassembler_{number}"), ".fasta", ".gfa", ".log")
    params:
        prefix=os.path.join(RESULTS_DIR, "intermediate", "assemblies", "plassembler_{number}")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        genome_size=$(cat {input.genome_size})

        autocycler helper plassembler \
            --reads {input.sub_seq} \
            --out_prefix {params.prefix} \
            --threads {threads} \
            --genome_size "$genome_size"
        """

rule autocycler_compress_assemblies:
    """
    Step 3: compress the input assemblies into a unitig graph
    """
    input: 
        flye=flye_list,
        plassembler=plassembler_list
    output:
        gfa=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "input_assemblies.gfa")
    params:
        indir=os.path.join(RESULTS_DIR, "intermediate", "assemblies"),
        outdir=directory(os.path.join(RESULTS_DIR, "intermediate", "autocycler_out"))
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        autocycler compress \
            -i {params.indir} \
            -a {params.outdir} \
            -t {threads}
        """

rule autocycler_cluster:
    """
    Step 4: cluster the input contigs into putative genomic sequences
    """
    input: 
        gfa=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "input_assemblies.gfa")
    output:
        os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "clustering", "clustering.tsv")
    params:
        indir=directory(os.path.join(RESULTS_DIR, "intermediate", "autocycler_out")),
        outdir=directory(os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "clustering"))
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        autocycler cluster -a {params.indir} 
        """

rule autocycler_trim_and_resolve:
    """
    Steps 5 and 6: trim and resolve each QC-pass cluster
    """
    input: 
       os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "clustering", "clustering.tsv")
    output:
        done=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "clustering", "qc_pass", ".done")
    params:
        dir=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "clustering", "qc_pass", "cluster_*")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        # the number of clusters and cluster names constantly change, so keep using the loop instead of submit a job for each cluster
        for c in {params.dir}; do
            autocycler trim -c "$c"
            autocycler resolve -c "$c"
        done  &&
        touch {output.done}
        """

rule autocycler_combine_resolved_clusters:
    """
    Step 7: combine resolved clusters into a final assembly
    """
    input: 
        done=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "clustering", "qc_pass", ".done")
    output:
        gfa=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "consensus_assembly.gfa"),
        fasta=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "consensus_assembly.fasta")
    params:
        autocycler_dir=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out"),
        final_gfa=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "clustering", "qc_pass", "cluster_*", "5_final.gfa")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        autocycler combine \
            -a {params.autocycler_dir} \
            -i {params.final_gfa}
        """

rule medaka_polishing:
    """
    polish genome assemblies with long reads
    """
    input: 
        long_filt=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "long_reads_trimmed.fastq.gz"),
        fasta=os.path.join(RESULTS_DIR, "intermediate", "autocycler_out", "consensus_assembly.fasta")
    output:
        genome=os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "consensus.fasta")
    params:
        mode=config["medaka"]["mode"],
        dir=os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "medaka.yml")
    shell:
        """
        medaka_consensus \
            -i {input.long_filt} \
            -d {input.fasta} \
            -o {params.dir} \
            -t {threads} \
            {params.mode} 
        """

rule extract_chromosome_genome_after_medaka:
    input: 
        genome=os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "consensus.fasta")
    output:
        chr=os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "chromosome_genome.fasta")
    params:
        script=os.path.join(dir["scripts"], "extract_longest_sequence.py")
    threads:
        config["resources"]["small_cpu"]
    shell:
        """
        {params.script} -i {input.genome} -o {output.chr}
        """

rule checkm2_medaka_polished_chromosome_genome:
    input: 
        done=os.path.join(dir["db"], "CheckM2_database", ".checkm2_db_done"),
        chr=os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "chromosome_genome.fasta")
    output:
        report=os.path.join(RESULTS_DIR, "long_reads_only", "genome_quality", "quality_report.tsv")
    params:
        dir=os.path.join(RESULTS_DIR, "long_reads_only", "genome_quality"),
        db=os.path.join(dir["db"], "CheckM2_database", "uniref100.KO.1.dmnd")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "checkm2.yml")
    shell:
        """
        checkm2 predict --threads {threads} \
            --input {input.chr} \
            --output-directory {params.dir} \
            --database_path {params.db} \
            --force 
        """

rule polypolish_polishing_alignment:
    """
    polishing genome assemblies with short reads. Step 1: align illumina reads to the assembly
    """
    input: 
        r1_filt=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_R1_trimmed.fastq.gz"),
        r2_filt=os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_R2_trimmed.fastq.gz"),
        fasta=os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "consensus.fasta")
    output:
        index=temp(multiext(os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "consensus.fasta"), ".amb", ".ann", ".bwt", ".pac", ".sa")),
        r1_aligned=temp(os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_R1.sam")),
        r2_aligned=temp(os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_R2.sam"))
    threads:
        config["resources"]["big_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1_filt} > {output.r1_aligned}
        bwa mem -t {threads} -a {input.fasta} {input.r2_filt} > {output.r2_aligned}
        """

rule polypolish_polishing_filter_alignment:
    """
    Step 2: run Polypolish's insert size filter to exclude alignments with incongruous insert sizes
    """
    input: 
        index=multiext(os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "consensus.fasta"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
        r1_aligned=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_R1.sam"),
        r2_aligned=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_R2.sam")
    output:
        r1_trimmed=temp(os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_trimmed_R1.sam")),
        r2_trimmed=temp(os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_trimmed_R2.sam"))
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "polypolish.yml")
    shell:
        """
        polypolish filter --in1 {input.r1_aligned} \
            --in2 {input.r2_aligned} \
            --out1 {output.r1_trimmed} \
            --out2 {output.r2_trimmed}
        """
        
rule polypolish_polishing_final:
    """
    Step 3: run Polypolish with the trimmed alignments
    """
    input: 
        fasta=os.path.join(RESULTS_DIR, "long_reads_only", "genome_sequences", "consensus.fasta"),
        r1_trimmed=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_trimmed_R1.sam"),
        r2_trimmed=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "alignments_trimmed_R2.sam")
    output:
        result=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_sequences", "short_reads_polished.fasta")
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "polypolish.yml")
    shell:
        """
        polypolish polish {input.fasta} {input.r1_trimmed} {input.r2_trimmed} > {output.result}
        """

rule extract_chromosome_genome_after_polypolish:
    input: 
        genome=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_sequences", "short_reads_polished.fasta")
    output:
        chr=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_sequences", "chromosome_genome.fasta")
    params:
        script=os.path.join(dir["scripts"], "extract_longest_sequence.py")
    threads:
        config["resources"]["small_cpu"]
    shell:
        """
        {params.script} -i {input.genome} -o {output.chr}
        """

rule checkm2_polypolish_polished_genome:
    input: 
        done=os.path.join(dir["db"], "CheckM2_database", ".checkm2_db_done"),
        chr=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_sequences", "chromosome_genome.fasta")
    output:
        os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_quality", "quality_report.tsv")
    params:
        dir=os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_quality"),
        db=os.path.join(dir["db"], "CheckM2_database", "uniref100.KO.1.dmnd")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "checkm2.yml")
    shell:
        """
        checkm2 predict --threads {threads} \
            --input {input.chr} \
            --output-directory {params.dir} \
            --database_path {params.db} \
            --force 
        """