# Ben Keith
# last updated 2020.05.07
# Furey Lab Pipeline 2020
# Snakemake 1.0

import os
shell.prefix("module load python/2.7.12; ")

configfile: "project_config.yaml"
[os.makedirs("results/" + str(i) + "/logs", exist_ok=True) for i in config["samples"]]

#genome variables
organism = config["analysis"]["organism"]
genomeBuild = config["analysis"]["genomeBuild"]

# Target to run whole workflow:
rule all:
    input:
        "results/multiqc/multiqc.html",
        expand("results/{sample}/logs/cleanup.log", sample=config["samples"])

# StarSalmon flag.s
rule quantification:
    input:
        expand("results/{sample}/salmon/quant.sf", sample=config["samples"]),
        expand("results/{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
                                                   sample=config["samples"]),
        expand("results/{sample}/salmon/{sample}.aligned.sorted.bam",
                                                   sample=config["samples"])
    output:
        touch("temp/starSalmon_run.flag")

rule QC:
    input:
        "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.rnaseq.html",
        "results/{sample}/QC/qualimap/{sample}.bamqc/{sample}.bamqc.html",
        "results/{sample}/QC/rseqc/{sample}.Aligned.sortedByCoord.out.summary.txt",
        "results/{sample}/QC/rseqc/{sample}.junctionSaturation_plot.pdf",
        "results/{sample}/{sample}.salmon/quant.sf"
    output:
        touch("temp/QC_complete.flag")

###################################-
#### Fetch and split SRA files ####
###################################


rule prefetch:
    output:
        "temp/{sample}/fastq/{sample}.sra"
    log:
        "results/{sample}/logs/prefetch.log"
    shell:
        """
        module load sratoolkit/2.10.1
        prefetch -L 5 -o {output} {wildcards.sample} > {log}
        """

rule fastq_dump:
    input:
        "temp/{sample}/fastq/{sample}.sra"
    output:
        "results/{sample}/fastq/{sample}_1.fastq.gz",
        "results/{sample}/fastq/{sample}_2.fastq.gz"
    params:
        fastq_dir = "results/{sample}/fastq/",
        other_flags = "--split-files"
    log:
        "results/{sample}/logs/fastq_dump.log"
    shell:
        """
        module load sratoolkit/2.10.1
        fastq-dump {params.other_flags} -L 5 -O {params.fastq_dir} {input} \
          > {log}
        gzip results/{wildcards.sample}/fastq/{wildcards.sample}_1.fastq
        gzip results/{wildcards.sample}/fastq/{wildcards.sample}_2.fastq
        """

##################################
#### Preprocessing - cutadapt ####
##################################

rule cutadapt:
    input:
        fastq1 = "results/{sample}/fastq/{sample}_1.fastq.gz",
        fastq2 = "results/{sample}/fastq/{sample}_2.fastq.gz"
    output:
        trimmed1 = "temp/{sample}/fastq/{sample}_1.fastq.trimmed.gz",
        trimmed2 = "temp/{sample}/fastq/{sample}_2.fastq.trimmed.gz"
    params:
        a = config["cutadapt"]["a"],
        A = config["cutadapt"]["A"],
        qualityCutoff = config["cutadapt"]["qualityCutoff"],
        minimumLength = config["cutadapt"]["minimumLength"],
        threads = config["cutadapt"]["threads"]
    log:
        "results/{sample}/logs/cutadapt.log"
    shell:
        """
        module load cutadapt/2.9
        cutadapt -a {params.a} -A {params.A} \
          --quality-cutoff {params.qualityCutoff} \
          --cores {params.threads} \
          --minimum-length {params.minimumLength} \
          -o {output.trimmed1} -p {output.trimmed2} \
          {input.fastq1} {input.fastq2} \
          > {log}
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
            fastq1 = "temp/{sample}/fastq/{sample}_1.fastq.trimmed.gz",
            fastq2 = "temp/{sample}/fastq/{sample}_2.fastq.trimmed.gz"
        output:
            "temp/{sample}.Aligned.sortedByCoord.out.bam",
            "results/{sample}/star/{sample}.SJ.out.tab"
        params:
            threads  = config["star"]["threads"],
            genomeLoad = config["star"]["genomeLoad"],
            outSAMtype = config["star"]["outSAMtype"],
            quantMode = config["star"]["quantMode"],
            readFilesCommand = config["star"]["readFilesCommand"],
            genomeDir = config[genomeBuild]["starIndex"],
            sjdbOverhang = config["star"]["sjdbOverhang"],
            outFileNamePrefix = "results/{sample}/star/{sample}.",
            featureFile = config[genomeBuild]["featureFile"]
        log:
            "results/{sample}/logs/star/star.log"
        shell:
            """
            mkdir -p results/{wildcards.sample}/star/
            module load star/2.7.3a
            star --runThreadN {params.threads} --sjdbOverhang {params.sjdbOverhang} \
              --genomeDir {params.genomeDir} --readFilesCommand {params.readFilesCommand} \
              --outSAMtype {params.outSAMtype} --outSAMunmapped Within \
              --quantMode {params.quantMode} \
              --outFileNamePrefix {params.outFileNamePrefix} \
              --readFilesIn {input.fastq1} {input.fastq2} \
              > {log}
            samtools index \
              results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam

            cp results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam* \
              temp
            cp {params.featureFile} temp; gunzip temp/*gtf.gz

            mv {params.outFileNamePrefix}Log.final.out \
              results/{wildcards.sample}/logs/star/
            mv {params.outFileNamePrefix}Log.progress.out \
              results/{wildcards.sample}/logs/star/
            mv {params.outFileNamePrefix}Log.out \
              results/{wildcards.sample}/logs/star/
            """

    rule salmon:
        input:
            fastq1 = "temp/{sample}/fastq/{sample}_1.fastq.trimmed.gz",
            fastq2 = "temp/{sample}/fastq/{sample}_2.fastq.trimmed.gz"
        output:
            "results/{sample}/{sample}.salmon/quant.sf",
            "temp/{sample}/salmon/{sample}.aligned.sam"
        params:
            index = config[genomeBuild]["salmonIndex"],
            libType = config["salmon"]["libType"],
            threads = config["salmon"]["threads"],
            numBootstraps = config["salmon"]["numBootstraps"],
            otherFlags = config["salmon"]["otherFlags"],
            outDir = "results/{sample}/{sample}.salmon",
            tmpDir = "temp/{sample}/salmon"
        log:
            "results/{sample}/logs/salmon.log"
        shell:
            """
            mkdir -p {params.outDir}
            mkdir -p {params.tmpDir}
            module load salmon/1.1.0
            salmon quant --libType {params.libType} {params.otherFlags} \
              --numBootstraps={params.numBootstraps} --threads {params.threads} \
              --writeMappings={params.tmpDir}/{wildcards.sample}.aligned.sam \
              -i {params.index} -o {params.outDir} \
              -1 {input.fastq1} -2 {input.fastq2} \
              > {log}
            mv {params.outDir}/logs/salmon_quant.log \
              results/{wildcards.sample}/logs
            rm -r {params.outDir}/logs
            """

    rule salmon_sam:
        input:
            "temp/{sample}/salmon/{sample}.aligned.sam"
        output:
            "results/{sample}/salmon/{sample}.aligned.sorted.bam",
            "results/{sample}/salmon/{sample}.aligned.sorted.bam.bai"
        params:
            outDir = "results/{sample}/salmon",
            tmpDir = "temp/{sample}/salmon"
        log:
            "results/{sample}/logs/salmon_sam.log"
        shell:
            """
            module load samtools/1.9
            samtools view -S -b {params.tmpDir}/{wildcards.sample}.aligned.sam \
              -o {params.tmpDir}/{wildcards.sample}.aligned.bam > {log}
            samtools sort {params.tmpDir}/{wildcards.sample}.aligned.bam \
              -o {params.outDir}/{wildcards.sample}.aligned.sorted.bam >> {log}
            samtools index {params.outDir}/{wildcards.sample}.aligned.sorted.bam \
              >> {log}
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

# rule bamCoverage
#     input:
#     output:
#     log:
#     shell:

#########################
#### Quality control ####
#########################

rule fastqc:
    input:
        "results/{sample}/fastq/{sample}_{pair}.fastq.gz"
    output:
        html = "results/{sample}/QC/fastqc/{sample}_{pair}_fastqc.html",
        zip = "results/{sample}/QC/fastqc/{sample}_{pair}_fastqc.zip"
    log:
        "results/{sample}/logs/fastqc/fastqc_{pair}.log"
    shell:
        """
        module load fastqc/0.11.8
        fastqc {input} -q -o . > {log}
        mv {wildcards.sample}_{wildcards.pair}_fastqc.html \
          results/{wildcards.sample}/QC/fastqc/
        mv {wildcards.sample}_{wildcards.pair}_fastqc.zip \
          results/{wildcards.sample}/QC/fastqc/
        """

os.makedirs("results/multiqc", exist_ok=True)
rule multiqc_raw:
    input:
        expand("results/{sample}/QC/fastqc/{sample}_{pair}_fastqc.zip", \
                sample = config["samples"], pair = ["1","2"]),
    output:
        html = "results/multiqc/multiqc_raw.html",
    params:
        configFile = config[genomeBuild]["multiqcConfig"]
    log:
        "results/multiqc/multiqc_raw.log"
    shell:
        """
        module load multiqc/1.7
        multiqc -n {output.html} -c {params.configFile} results/*/QC/fastqc \
          temp/fastqc > {log}
        mv results/multiqc/multiqc_data results/multiqc/multiqc_raw_data
        """

rule fastQscreen:
    input:
        "results/{sample}/fastq/{sample}_{pair}.fastq.gz"
    output:
        "results/{sample}/QC/fastQscreen/{sample}_{pair}_screen.txt"
    params:
        aligner = config["fastQscreen"]["aligner"],
        ourDir = "results/{sample}/QC/fastQscreen",
        subset = config["fastQscreen"]["subset"],
        conf = config["fastQscreen"]["conf"]
    log:
        "results/{sample}/logs/fastQscreen.log"
    shell:
        """
        /proj/fureylab/bin/fastq_screen --conf {params.conf} \
          --aligner {params.aligner} --outdir {params.outDir} \
          --subset {params.subset} {input} > {log}
        """

rule rseqc:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/rseqc/{sample}.Aligned.sortedByCoord.out.summary.txt",
        "results/{sample}/QC/rseqc/{sample}.junctionSaturation_plot.pdf"
    params:
        model = config[genomeBuild]["rseqcModel"]
    log:
        "results/{sample}/logs/rseqc.log"
    shell:
        """
        module load rseqc/3.0.1
        mkdir -p results/{wildcards.sample}/QC/rseqc
        tin.py -i {input} -r {params.model}
        junction_saturation.py -i {input} -r {params.model} \
          -o {wildcards.sample}
        mv {wildcards.sample}* results/{wildcards.sample}/QC/rseqc
        """

rule QM_bamqc:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/qualimap/{sample}.bamqc/{sample}.bamqc.html",
        "results/{sample}/QC/qualimap/{sample}.bamqc/{sample}.genomeCoverage.txt"
    params:
        outDir = "results/{sample}/QC/qualimap/{sample}.bamqc",
        outCoverage = "results/{sample}/QC/qualimap/{sample}.bamqc/{sample}.genomeCoverage.txt",
        outHTML = "results/{sample}/QC/qualimap/{sample}.bamqc/{sample}.bamqc.html",
        seqProtocol = "strand-specific-reverse",
        featureFile = "temp/*gtf",
        genomeGC = organism,
        javaMemSize = "4G",
        otherFlags = "--collect-overlap-pairs"
    log:
        "results/{sample}/logs/QM_bamqc.log"
    shell:
        """
        module load qualimap/2.2.1
        mkdir -p {params.outDir}
        qualimap bamqc -bam {input} -outdir {params.outDir} \
          -oc {params.outCoverage} -outfile {params.outHTML} -outformat HTML \
          --sequencing-protocol {params.seqProtocol} \
          --feature-file {params.featureFile} --genome-gc-distr {params.organism} \
          --java-mem-size={params.javaMemSize} {params.otherFlags}
        """

rule QM_rnaseq:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.computedCounts.txt",
        "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.rnaseq.html"
    params:
        outDir = "results/{sample}/QC/qualimap/{sample}.rnaseq",
        outCounts = "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.computedCounts.txt",
        outHTML = "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.rnaseq.html",
        seqProtocol = "strand-specific-reverse",
        featureFile = "temp/*gtf",
        javaMemSize = "4G",
        otherFlags = "--paired --sorted"
    log:
        "results/{sample}/logs/QM_rnaseq.log"
    shell:
        """
        module load qualimap/2.2.1
        mkdir -p {params.outDir}
        qualimap rnaseq -bam {input} -outdir {params.outDir} -outformat HTML \
        -oc {params.outCounts} -outfile {params.outHTML} \
        --sequencing-protocol {params.seqProtocol} -gtf {params.featureFile} \
        --java-mem-size={params.otherFlags} {params.otherFlags}
        """

os.makedirs("results/multiqc", exist_ok=True)
rule multiqc:
    input:
        "temp/QC_complete.flag"
    output:
        html = "results/multiqc/multiqc.html",
        data = directory("results/multiqc/multiqc_data")
    params:
        configFile = config[genomeBuild]["multiqcConfig"]
    log:
        "results/multiqc/multiqc.log"
    shell:
        """
        module load multiqc/1.7
        multiqc -n {output.html} -c {params.configFile} . > {log}
        """

##################
#### CLEAN UP ####
##################

# For the purposes of multiqc, some directories had to be sample specific
# This rule fixes those names.
rule name_clean:
    input:
        "results/multiqc/multiqc.html",
        expand("results/{sample}}/{sample}.salmon", sample=config["samples"]),
        "results/{sample}/QC/qualimap/{sample}.rnaseq",
        "results/{sample}/QC/qualimap/{sample}.bamqc"
    output:
        directory("results/{sample}}/salmon"),
        directory("results/{sample}/QC/qualimap/rnaseq"),
        directory("results/{sample}/QC/qualimap/bamqc")
    log:
        "results/{sample}/logs/cleanup.log"
    shell:
        """
        mv results/{sample}}/{sample}.salmon \
          results/{sample}}/salmon
        mv results/{sample}/QC/qualimap/{sample}.rnaseq \
          results/{sample}/QC/qualimap/rnaseq
        mv results/{sample}/QC/qualimap/{sample}.bamqc \
          results/{sample}/QC/qualimap/bamqc
        """
