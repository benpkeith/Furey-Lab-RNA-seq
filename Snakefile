# Ben Keith
# AstraZeneca RNA-seq pipeline
import os
shell.prefix("module load python/2.7.12; ")

configfile: "config.yaml"

# Target to run whole workflow:
rule all:
    input:
        expand("results/rsem/{sample}.genome.sorted.bam.bai", \
                                        sample = config["samples"]),
        "results/multiqc/multiqc.html"

###################################-
#### Fetch and split SRA files ####
###################################

rule fetch_SRA:
    output:
        "raw/fastq/{sample}.sra"
    log: "logs/prefetch/{sample}.log"
    shell:
        """
        module load sratoolkit/2.10.1
        prefetch -L 5 -o raw/fastq/{wildcards.sample}.sra {wildcards.sample} \
        >{log}
        """

rule split_SRA:
    input: "raw/fastq/{sample}.sra"
    output:
        "raw/fastq/{sample}_1.fastq.gz",
        "raw/fastq/{sample}_2.fastq.gz"
    params:
        fastq_dir = "raw/fastq/",
        other_flags = "--split-files"
    log: "logs/fastq_dump/{sample}.log"
    shell:
        """
        module load sratoolkit/2.10.1
        fastq-dump {params.other_flags} \
                   -L 5 -O {params.fastq_dir} \
                   {input} >{log}
        gzip raw/fastq/{wildcards.sample}_1.fastq
        gzip raw/fastq/{wildcards.sample}_2.fastq
        """

##################################################
#### Alignment and quantification - STAR/rsem ####
##################################################

rule rsem:
    input:
        fastq1 = "raw/fastq/{sample}_1.fastq.gz",
        fastq2 = "raw/fastq/{sample}_2.fastq.gz"
    output:
        "results/rsem/{sample}.genes.results",
        directory("temp/rsem/{sample}"),
        "results/rsem/{sample}.genome.bam"
    params:
        flags = "--paired-end --star --output-genome-bam \
                --keep-intermediate-files  --star-gzipped-read-file",
        tempFolder = "temp/rsem/{sample}",
        outPrefix = "results/rsem/{sample}",
        threads = config["rsem"]["threads"],
        reference = config["rsem"]["reference"]
    log:
        "logs/rsem/{sample}.log"
    shell:
        """
        module load star/2.7.0e
        module load rsem/1.2.31
        rsem-calculate-expression {params.flags} --num-threads {params.threads} \
        --temporary-folder {params.tempFolder} \
        {input.fastq1} {input.fastq2} {params.reference} \
         {params.outPrefix} >{log}
        """

rule rsem_bam_sort:
    input:
        results = "results/rsem/{sample}.genes.results",
        bam = "results/rsem/{sample}.genome.bam"
    output:
        "results/rsem/{sample}.genome.sorted.bam"
    log:
        "logs/samtools_sort/{sample}.log"
    shell:
        """
        module load samtools/1.9
        samtools sort {input.bam} -o {output} >{log}
        """

rule rsem_bam_index:
    input:
        "results/rsem/{sample}.genome.sorted.bam"
    output:
        "results/rsem/{sample}.genome.sorted.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    shell:
        """
        module load samtools/1.9
        samtools index {input}
        """

#########################
#### Quality control ####
#########################

rule fastqc:
    input:
        "raw/fastq/{sample}_{pair}.fastq.gz"
    output:
        html = "results/fastqc/{sample}_{pair}_fastqc.html",
        zip = "temp/fastqc/{sample}_{pair}_fastqc.zip"
    shell:
        """
        module load fastqc/0.11.8
        fastqc {input} -q -o .
        mv {wildcards.sample}_{wildcards.pair}_fastqc.html {output.html}
        mv {wildcards.sample}_{wildcards.pair}_fastqc.zip {output.zip}
        """

rule rseqc_tin:
    input:
        "results/rsem/{sample}.genome.sorted.bam.bai"
    output:
        "results/rseqc/{sample}.genome.sorted.summary.txt"
    params:
        model = config["rseqc"]["modelFile"]
    log:
        "logs/rseqc/{sample}.log"
    shell:
        """
        tin.py -r {params.model} -i \
        results/rsem/{wildcards.sample}.genome.sorted.bam >{log}
        mv {wildcards.sample}*genome* results/rseqc
        """

rule rseqc_geneBodyCoverage:
    input:
        expand("results/rseqc/{sample}.genome.sorted.summary.txt",
                sample = config["samples"])
    output:
        touch("temp/rseqc_coverage_done.flag")
    params:
        model = config["rseqc"]["modelFile"],
        analysis = config["analysis"]["name"]
    log:
        "logs/rseqc/geneBodyCoverage.log"
    shell:
        """
        geneBody_coverage.py -r {params.model} \
        -i results/rsem/ -o {params.analysis} >{log}
        mv {params.analysis}* results/rseqc
        mv log.txt logs/rseqc
        """

rule multiqc:
    input:
        expand("temp/fastqc/{sample}_{pair}_fastqc.zip", \
                sample = config["samples"], pair = ["1","2"]),
        "temp/rseqc_coverage_done.flag"
    output:
        html = "results/multiqc/multiqc.html",
        data = directory("results/multiqc/multiqc_data")
    params:
        configFile = config["multiqc"]["configFile"]
    log:
        "logs/multiqc/multiqc.log"
    shell:
        """
        module load multiqc/1.7
        multiqc -n {output.html} -c {params.configFile} . > {log}
        """
