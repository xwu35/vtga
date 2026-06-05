#------------ SET UP THE DIRECTORIES
dir = dict()
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")
dir["db"]      = os.path.join(workflow.basedir, "..", "db")

#------------ GET LIST FOR RULES
flye_list={}
plassembler_list={}
fastp_result=[]
short_assembly=[]
long_quality_report=[]
hybrid_quality_report=[]

for genome in GENOME:
    # get the list of fasta files from flye assembly 
    flye_list[genome]=[os.path.join(RESULTS_DIR, "intermediate", "autocycler_output", "assemblies", genome, "flye_" + n + ".fasta").format(n) for n in NUMBER]
    # get the list of marker files from plassembler assembly
    plassembler_list[genome]=[os.path.join(RESULTS_DIR, "intermediate", "autocycler_output", "assemblies", genome, ".finished_plassembler_" + n).format(n) for n in NUMBER]
    # output
    fastp_result.append(os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", genome + "_sr_fastp.html"))
    fastp_result.append(os.path.join(RESULTS_DIR, "intermediate", "trimmed_reads", genome + "_sr_fastp.json"))
    short_assembly.append(os.path.join(RESULTS_DIR, "short_reads_only", "genome_sequences", genome + "_contigs.gfa"))
    short_assembly.append(os.path.join(RESULTS_DIR, "short_reads_only", "genome_quality", genome, "quality_report.tsv"))
    long_quality_report.append(os.path.join(RESULTS_DIR, "long_reads_only", "genome_quality", genome, "quality_report.tsv"))
    hybrid_quality_report.append(os.path.join(RESULTS_DIR, "long_reads_plus_short_reads_polished", "genome_quality", genome, "quality_report.tsv"))

#------------ SET UP THE OUTPUT
# only short reads provided
short_only_input = fastp_result + short_assembly

# only long reads are provided
long_only_input = long_quality_report + [
    os.path.join(RESULTS_DIR, "long_reads_only", "software_versions.txt"),
]

# both short and long reads are provided
hybrid_input = fastp_result + long_only_input + hybrid_quality_report

#------------ DOWNLOAD DATABASES
rule download_plassembler_db:
    conda:
        os.path.join(dir["env"], "autocycler.yml")
    output:
        os.path.join(CONDA_ENVS, ".plassembler_db_done")
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