# Ben Keith
# last updated 2020.03.25
# Furey Lab Pipeline 2020
# Snakemake 1.0

import os
shell.prefix("module load python/2.7.12; ")

configfile: "project_config.yaml"

#genome variables
organism = config["analysis"]["organism"]
genomeBuild = config["analysis"]["genomeBuild"]

# Target to run whole workflow:
rule all:
    input:
        expand("results/rsem/{sample}.genome.sorted.bam.bai", \
                                        sample = config["samples"]),
        "results/multiqc/multiqc.html"

###################################-
#### Fetch and split SRA files ####
###################################


rule prefetch:
    output:
        "raw/fastq/{sample}.sra"
    log:
        "results/{sample}/logs/prefetch.log"
    shell:
        """
        module load sratoolkit/2.10.1
        prefetch -L 5 -o raw/fastq/{wildcards.sample}.sra {wildcards.sample}
        """

rule fastq_dump:
    input:
        "raw/fastq/{sample}.sra"
    output:
        "raw/fastq/{sample}_1.fastq.gz",
        "raw/fastq/{sample}_2.fastq.gz"
    params:
        fastq_dir = "raw/fastq/",
        other_flags = "--split-files"
    log:
        "results/{sample}/logs/fastq_dump.log"
    shell:
        """
        module load sratoolkit/2.10.1
        fastq-dump {params.other_flags} -L 5 -O {params.fastq_dir} {input}
        gzip raw/fastq/{wildcards.sample}_1.fastq
        gzip raw/fastq/{wildcards.sample}_2.fastq
        """

##################################
#### Preprocessing - cutadapt ####
##################################

rule cutadapt:
    input:
        fastq1 = "data/fastq/{sample}_1.fastq.gz",
        fastq2 = "data/fastq/{sample}_2.fastq.gz"
    output:
        trimmed1 = "data/fastq/{sample}_1.fastq.trimmed.gz",
        trimmed2 = "data/fastq/{sample}_2.fastq.trimmed.gz"
    params:
        a = config["cutadapt"]["a"]
        A = config["cutadapt"]["A"]
        qualityCutoff = config["cutadapt"]["qualityCutoff"]
        minimumLength = config["cutadapt"]["minimumLength"]
    log:
        "results/{sample}/logs/cutadapt.log"
    shell:
        """
        module load cutadapt/2.9
        cutadapt -a {params.a} -A {params.A} \
        --quality-cutoff {params.qualityCutoff} --minimum-length {params.minimumLength} \
        -o {output.trimmed1} -p {output.trimmed2} \
        {input.fastq1} {input.fastq2}
        """

##################################################
#### Alignment and quantification - STAR/rsem ####
##################################################

# if config["quantification"] == "rsem":
#     rule rsem:
#         input:
#             fastq1 = "raw/fastq/{sample}_1.fastq.gz",
#             fastq2 = "raw/fastq/{sample}_2.fastq.gz"
#         output:
#             "results/rsem/{sample}.genes.results",
#             directory("temp/rsem/{sample}"),
#             "results/rsem/{sample}.genome.bam"
#         params:
#             otherFlags = config["rsem"]["otherFlags"],
#             tempFolder = "temp/rsem/{sample}",
#             outPrefix = "results/rsem/{sample}",
#             threads = config["rsem"]["threads"],
#             reference = config["rsem"]["reference"]
#         log:
#             "results/{sample}/logs/rsem/{sample}.log"
#         shell:
#             """
#             module load star/2.7.0e
#             module load rsem/1.2.31
#             rsem-calculate-expression {params.otherFlags} --num-threads {params.threads} \
#               --temporary-folder {params.tempFolder} \
#               {input.fastq1} {input.fastq2} {params.reference} \
#               {params.outPrefix} >{log}
#             """
#
#     rule rsem_bam_sort:
#         input:
#             results = "results/rsem/{sample}.genes.results",
#             bam = "results/rsem/{sample}.genome.bam"
#         output:
#             "results/rsem/{sample}.genome.sorted.bam"
#         log:
#             "results/{sample}/logs/samtools_sort/{sample}.log"
#         shell:
#             """
#             module load samtools/1.9
#             samtools sort {input.bam} -o {output} >{log}
#             """
#
#     rule rsem_bam_index:
#         input:
#             "results/rsem/{sample}.genome.sorted.bam"
#         output:
#             "results/rsem/{sample}.genome.sorted.bam.bai"
#         log:
#             "results/{sample}/logs/samtools_index/{sample}.log"
#         shell:
#             """
#             module load samtools/1.9
#             samtools index {input}
#             """

####################################################
#### Alignment and quantification - STAR/salmon ####
####################################################

if config["quantification"] == "salmon":
    rule star:
        input:
            fastq1 = "data/fastq/{sample}_1.fastq.trimmed.gz",
            fastq2 = "data/fastq/{sample}_2.fastq.trimmed.gz"
        output:
            "results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            "results/star/{sample}/{sample}.SJ.out.tab"
        params:
            theads  = config["star"]["threads"]
            genomeLoad = config["star"]["genomeLoad"]
            outSAMtype = config["star"]["outSAMtype"]
            quantMode = config["star"]["quantMode"]
            readFilesCommand = config["star"]["readFilesCommand"]
            genomeDir = config[organism][genomeBuild]["starIndex"]
            sjdbOverhang = config["star"]["sjdbOverhang"]
            outFileNamePrefix = "results/{sample}/star/{sample}."
        log:
            "results/{sample}/logs/star/star.log"
        shell:
            """
            module load star/2.7.3a
            star --runThreadN {params.threads} --sjdbOverhang {params.sjdbOverhang} \
            --genomeDir {params.genomeDir} --readFilesCommand {params.readFilesCommand} \
            --outSAMtype {params.outSAMtype} --outSAMunmapped Within \
            --outFileNamePrefix {params.outFileNamePrefix} \
            --readFilesIn {input.fastq1} {input.fastq2}
            mv {params.outFileNamePrefix}.Log.final.out results/{sample}/logs/star/
            mv {params.outFileNamePrefix}.Log.progress.out results/{sample}/logs/star/
            mv {params.outFileNamePrefix}.Log.out results/{sample}/logs/star/
            """

    rule salmon
        input:
            fastq1 = "data/fastq/{sample}_1.fastq.trimmed.gz",
            fastq2 = "data/fastq/{sample}_2.fastq.trimmed.gz"
        output:
            "results/{sample}/salmon/quant.sf",
            "results/{sample}/salmon/{sample}.aligned.sam"
        params:
            index = config[organism][genomeBuild]["salmonIndex"]
            libType = config["salmon"]["libType"]
            threads = config["salmon"]["threads"]
            numBootstraps = config["salmon"]["numBootstraps"]
            otherFlags = config["rsem"]["otherFlags"]
            outDir = "results/{sample}/salmon"
        log:
            "results/{sample}/logs/salmon.log"
        shell:
            """
            module load salmon/1.1.0
            salmon quant --libType {params.libType} {params.otherFlags} \
            --numBootstraps={params.numBootstraps} --threads {params.threads} \
            --writeMappings={params.outDir}/{wildcards.sample}.aligned.sam \
            -i {params.index} -o {params.outDir} \
            -1 {input.fastq1} -2 {input.fastq2}
            mv {params.outDir}/logs/salmon_quant.log results/{sample}/logs
            rm -r {params.outDir}/logs
            """

#############################
#### Post-quantification ####
#############################

# rule deconvolution
#     input:
#     output:
#     log:
#     shell:
#
# rule matrixGeneration
#     input:
#     output:
#     log:
#     shell:
#
# rule bam_sort
#         input:
#             results = "results/STAR/{sample}.genes.results",
#             bam = "results/STAR/{sample}.genome.bam"
#         output:
#             "results/STAR/{sample}.genome.sorted.bam"
#         log:
#             "results/{sample}/logs/samtools_sort/{sample}.log"
#         shell:
#             """
#             module load samtools/1.9
#             samtools sort {input.bam} -o {output} >{log}
#             """
#
# rule bam_index
#         input:
#             "results/STAR/{sample}.genome.sorted.bam"
#         output:
#             "results/STAR/{sample}.genome.sorted.bam.bai"
#         log:
#             "results/{sample}/logs/samtools_index/{sample}.log"
#         shell:
#             """
#             module load samtools/1.9
#             samtools index {input}
#             """
#
# rule bamCoverage
#     input:
#     output:
#     log:
#     shell:

#########################
#### Quality control ####
#########################

# rule fastqc:
#     input:
#         "raw/fastq/{sample}_{pair}.fastq.gz"
#     output:
#         html = "results/fastqc/{sample}_{pair}_fastqc.html",
#         zip = "temp/fastqc/{sample}_{pair}_fastqc.zip"
#     shell:
#         """
#         module load fastqc/0.11.8
#         fastqc {input} -q -o .
#         mv {wildcards.sample}_{wildcards.pair}_fastqc.html {output.html}
#         mv {wildcards.sample}_{wildcards.pair}_fastqc.zip {output.zip}
#         """
#
# rule multiqc_raw:
#     input:
#         expand("temp/fastqc/{sample}_{pair}_fastqc.zip", \
#                 sample = config["samples"], pair = ["1","2"]),
#     output:
#         html = "results/multiqc/multiqc_raw.html",
#     params:
#         configFile = config["multiqc"]["configFile"]
#     log:
#         "results/{sample}/logs/multiqc/multiqc.log"
#     shell:
#         """
#         module load multiqc/1.7
#         multiqc -n {output.html} -c {params.configFile} results/fastqc \
#         temp/fastqc > {log}
#         """
#
# #rule virusSeq
#
# rule rseqc_tin:
#     input:
#         "results/rsem/{sample}.genome.sorted.bam.bai"
#     output:
#         "results/rseqc/{sample}.genome.sorted.summary.txt"
#     params:
#         model = config["rseqc"]["modelFile"]
#     log:
#         "results/{sample}/logs/rseqc/{sample}.log"
#     shell:
#         """
#         tin.py -r {params.model} -i \
#         results/rsem/{wildcards.sample}.genome.sorted.bam >{log}
#         mv {wildcards.sample}*genome* results/rseqc
#         """
#
# rule rseqc_geneBodyCoverage:
#     input:
#         expand("results/rseqc/{sample}.genome.sorted.summary.txt",
#                 sample = config["samples"])
#     output:
#         touch("temp/rseqc_coverage_done.flag")
#     params:
#         model = config["rseqc"]["modelFile"],
#         analysis = config["analysis"]["name"]
#     log:
#         "results/{sample}/logs/rseqc/geneBodyCoverage.log"
#     shell:
#         """
#         geneBody_coverage.py -r {params.model} \
#           -i results/rsem/ -o {params.analysis} >{log}
#         mv {params.analysis}* results/rseqc
#         mv log.txt results/{sample}/logs/rseqc
#         """
#
# rule rseqc_readDistribution
#     input:
#     output:
#     log:
#     shell:
#
# rule rseqc_junctionSaturation
#     input:
#     output:
#     log:
#     shell:
#
# rule rseqc_rRNAcontamination
#     input:
#     output:
#     log:
#     shell:
#
# rule multiqc:
#     input:
#         expand("temp/fastqc/{sample}_{pair}_fastqc.zip", \
#                 sample = config["samples"], pair = ["1","2"]),
#         "temp/rseqc_coverage_done.flag"
#     output:
#         html = "results/multiqc/multiqc.html",
#         data = directory("results/multiqc/multiqc_data")
#     params:
#         configFile = config["multiqc"]["configFile"]
#     log:
#         "results/{sample}/logs/multiqc/multiqc.log"
#     shell:
#         """
#         module load multiqc/1.7
#         multiqc -n {output.html} -c {params.configFile} . > {log}
#         """
