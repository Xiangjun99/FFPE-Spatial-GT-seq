############################################
# Bulk DNA Snakemake pipeline
############################################

############################################
# Paths
############################################

#change the path
FASTQ_DIR = "/path/to/FASTQ_files"
bowtie2_index = "/path/to/GRCh38_noalt_decoy_as" #for Human
bowtie2_index = "/path/to/mm10" #for mouse

bowtie2_path = "bowtie2"
samtools_path = "samtools"
sambamba_path = "sambamba"

############################################
# Samples
############################################

(samples,) = glob_wildcards(
    f"{FASTQ_DIR}/{{sample}}_R1.fastq.gz"
)

samples = sorted(samples)

############################################
# Rules
############################################

rule all:
    input:
        expand("marked/{sample}.bam", sample=samples)


############################################
# Mapping
############################################

rule bowtie2:
    input:
        r1 = f"{FASTQ_DIR}/{{sample}}_R1.fastq.gz",
        r2 = f"{FASTQ_DIR}/{{sample}}_R2.fastq.gz"
    output:
        temp("mapped/{sample}.bam")
    log:
        "logs/bowtie2/{sample}.log"
    threads: 8
    shell:
        """
        ({bowtie2_path} \
            -x {bowtie2_index} \
            -p {threads} \
            -1 {input.r1} \
            -2 {input.r2} | \
         {samtools_path} view -Sb -@ {threads} > {output}) \
         2> {log}
        """


############################################
# Sort BAM
############################################

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        temp("sort/{sample}.bam")
    threads: 4
    shell:
        """
        {samtools_path} sort {input} -@ {threads} -o {output}
        """


############################################
# Index BAM
############################################

rule index:
    input:
        "sort/{sample}.bam"
    output:
        temp("sort/{sample}.bam.bai")
    shell:
        """
        {samtools_path} index {input}
        """


############################################
# Mark duplicates
############################################

rule sambamba_markdup:
    input:
        "sort/{sample}.bam"
    output:
        "marked/{sample}.bam"
    threads: 4
    shell:
        """
        {sambamba_path} markdup -t {threads} {input} {output}
        """
