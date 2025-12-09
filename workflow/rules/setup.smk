#------------ SET UP THE DIRECTORIES
dir = dict()
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")
dir["db"]      = os.path.join(workflow.basedir, "..", "db")
dir["conda_envs"] = os.path.join(workflow.basedir, "..", "conda_envs")

#------------ GET LIST FOR RULES
# get the list of subsampled reads
reads_list=[os.path.join(RESULTS_DIR, "intermediate", "subsampled_reads", "sample_" + n + ".fastq").format(n) for n in NUMBER]

# get the list of fasta files from flye assembly
flye_list=[os.path.join(RESULTS_DIR, "intermediate", "assemblies", "flye_" + n + ".fasta").format(n) for n in NUMBER] 

# get the list of fasta files from plassembler assembly
plassembler_list=[os.path.join(RESULTS_DIR, "intermediate", "assemblies", "plassembler_" + n + ".fasta").format(n) for n in NUMBER]

#------------ SET UP THE OUTPUT
# only short reads provided
short_only_input = [
    os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_fastp.html"),
    os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_fastp.json"),
    os.path.join(RESULTS_DIR, "short_reads_only", "genome_sequences", "contigs.fa"),
    os.path.join(RESULTS_DIR, "short_reads_only", "genome_quality", "quality_report.tsv")
]

# only long reads are provided
long_only_input = [
    os.path.join(RESULTS_DIR, "long_reads_only", "genome_quality", "quality_report.tsv"),
    os.path.join(RESULTS_DIR, "long_reads_only", "software_versions.txt"),
]

# both short and long reads are provided
hybrid_input = [
    os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_fastp.html"),
    os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", "short_reads_fastp.json") 
] + long_only_input + [
    os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_quality", "quality_report.tsv")
]

#------------ DOWNLOAD DATABASES
rule download_plassembler_db:
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    output:
        os.path.join(dir["conda_envs"], ".plassembler_db_done")
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    shell:
        """
        plassembler download -d "$CONDA_PREFIX"/plassembler_db &&
        touch {output}
        """

rule download_checkm2_db:
    conda:
        os.path.join(dir["env"], "checkm2.yml")
    output:
        os.path.join(dir["db"], "CheckM2_database", ".checkm2_db_done")
    params:
        dir=dir["db"]
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    shell:
        """
        checkm2 database --download --path {params.dir} &&
        touch {output}
        """