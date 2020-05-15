# Ben Keith
# last updated 2020.05.07
# Furey Lab Pipeline 2020
# Snakemake 1.0

########################
#### Initial set up ####
########################

import os
import glob
shell.prefix("module load python/2.7.12; ")

configfile: "project_config.yaml"
# [os.makedirs("results/" + str(i) + "/logs", exist_ok=True) for i in config["samples"]]

#genome variables
configFilename = "project_config.yaml"
analysisName = config["analysis"]["name"]
organism = config["analysis"]["organism"]
genomeBuild = config["analysis"]["genomeBuild"]

#sample/path variables

if config["useSRA"]:
    samples = config["samples"]
    [os.makedirs("results/" + str(sample) + "/logs", exist_ok=True) for sample in samples]
else:
    samplePaths = config["samples"]
    samplePaths = [i[:-1] for i in samplePaths if i.endswith('/')]
    samples = [i.rsplit('/', 1)[1] for i in samplePaths]
    # fastqPaths = [str(i + '/fastq') for i in paths]
    #print(samples)
    [os.makedirs("results/" + str(sample) + "/logs", exist_ok=True) for sample in samples]
    [os.makedirs("results/" + str(sample) + "/fastq", exist_ok=True) for sample in samples]

    for path in samplePaths:
        samp = str(path.rsplit('/', 1)[1])
        if config["moveOutFiles"]:
            os.makedirs(str(path) + "/snakemakeRNA", exist_ok=True)
            print("Moving final files...")
            cmd = "cp -rf results/" + str(samp) + "/* " + str(path) + "/snakemakeRNA"
            os.system(cmd)
            print("moved sample " + str(samp))
            print("Files moved! Exiting...")

        for file in glob.glob(str(path + "/fastq/*gz")):

            if "_R1_" in file or file.endswith("_1.f*q.gz"):
                cmd = "ln -s " + str(file) + " results/" + \
                  samp + "/fastq/" + samp + "_1.fastq.gz" +  " >/dev/null 2>&1"
                #print(cmd)
                os.system(cmd)
            else:
                cmd = "ln -s " + str(file) + " results/" + \
                  samp + "/fastq/" + samp + "_2.fastq.gz" +  " >/dev/null 2>&1"
                #print(cmd)
                os.system(cmd)

    # paths = config["paths"]
    # fastqPaths = [str(i + '/fastq') for i in paths]
    # paths = [i[:-1] for i in paths if i.endswith('/')]
    # samples = [i.rsplit('/', 1)[1] for i in paths]
    #
    # [os.makedirs("results/" + str(i) + "/logs", exist_ok=True) for i in samples]
    # [os.makedirs("results/" + str(i) + "/fastq", exist_ok=True) for i in samples]
    #
    # paths = [i.rsplit('/', 1)[0] for i in paths]
    # for path in fastqPaths:
    #     for file in glob.glob(str(path + "/*gz")):
    #         os.system('ln -s file ')

# print(fastqPaths)
# print(paths)
# print(samples)


# [os.makedirs("results/" + str(i) + "/logs", exist_ok=True) for i in samples]

##################
#### Pipeline ####
##################

####################
#### all target ####
####################

# Target to run whole workflow:
rule all:
    input:
        "results/multiqc/multiqc.html",
        "results/counts/" + analysisName + ".txi.rds",
        expand("results/{sample}/logs/cleanup.log", sample=samples)

##########################
#### Checkpoint rules ####
##########################

# StarSalmon flag.s
rule alignmentQuantification:
    input:
        "results/{sample}/{sample}.salmon/quant.sf",
        "results/{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        "results/{sample}/{sample}.salmon/{sample}.aligned.sorted.bam"
    output:
        touch("temp/starSalmon_run.flag")

rule QC:
    input:
        expand("results/{sample}/QC/qualimap/{sample}.rnaseq/qualimapReport.html",
            sample=samples),
        expand("results/{sample}/QC/qualimap/{sample}.bamqc/qualimapReport.html",
            sample=samples),
        expand("results/{sample}/QC/rseqc/{sample}.Aligned.sortedByCoord.out.summary.txt",
            sample=samples),
        expand("results/{sample}/QC/rseqc/{sample}.junctionSaturation_plot.pdf",
            sample=samples),
        expand("results/{sample}/{sample}.salmon/quant.sf",
            sample=samples),
        expand("results/{sample}/QC/fastQscreen/{sample}_{pair}_screen.txt",
            sample=samples, pair=["1","2"]),
        "results/multiqc_raw/multiqc_raw.html"
    output:
        touch("temp/QC_complete.flag")

#################################################
####### useSRA: Fetch and split SRA files #######
#### else: Create links to local fastq files ####
#################################################

if config["useSRA"]:
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
            fastq_dir = "fastq/",
            other_flags = "--split-files"
        log:
            "results/{sample}/logs/fastq_dump.log"
        shell:
            """
            module load sratoolkit/2.10.1
            fastq-dump {params.other_flags} -L 5 -O {params.fastq_dir} {input} \
              > {log}
            gzip fastq/{wildcards.sample}_1.fastq
            gzip fastq/{wildcards.sample}_2.fastq
            mkdir results/{wildcards.sample}/fastq
            ln -s fastq/{wildcards.sample}_1.fastq.gz \
              results/{wildcards.sample}/fastq/{wildcards.sample}_1.fastq.gz
            ln -s fastq/{wildcards.sample}_2.fastq.gz \
              results/{wildcards.sample}/fastq/{wildcards.sample}_2.fastq.gz
            """

# if config["fastqR1"]:
#     rule fastq_links:
#         input:
#             R1 = expand("{path}/{sample}/fastq/{prepend}_R1_{sampleInfo}.f{fastq}q.gz",
#                      path = paths, sample = samples, allow_missing=True),
#             R2 = expand("{path}/{sample}/fastq/{prepend}_R2_{sampleInfo}.f{fastq}q.gz",
#                      path = paths, sample = samples, allow_missing=True)
#         output:
#             "results/{sample}/fastq/{sample}_1.fastq.gz",
#             "results/{sample}/fastq/{sample}_2.fastq.gz"
#         log:
#             "results/{sample}/logs/fastq_links.log"
#         shell:
#             """
#             ln -s {input.R1} \
#               results/{wildcards.sample}/fastq/{wildcards.sample}_1.fastq.gz > {log}
#             ln -s {input.R2} \
#               results/{wildcards.sample}/fastq/{wildcards.sample}_2.fastq.gz >> {log}
#               """
#
# else:
#     rule fastq_links:
#         input:
#             R1 = "{path}/{sample}/fastq/*1.f*q*.gz",
#             R2 = "{path}/{sample}/fastq/*2.f*q*.gz"
#         output:
#             "results/{sample}/fastq/{sample}_1.fastq.gz",
#             "results/{sample}/fastq/{sample}_2.fastq.gz"
#         log:
#             "results/{sample}/logs/fastq_links.log"
#         shell:
#             """
#             ln -s {input.R1} \
#               results/{wildcards.sample}/fastq/{wildcards.sample}_1.fastq.gz > {log}
#             ln -s {input.R2} \
#               results/{wildcards.sample}/fastq/{wildcards.sample}_2.fastq.gz >> {log}
#             """

# else:
#     rule fastq_links:
#         output:
#             "results/{sample}/fastq/{sample}_1.fastq.gz",
#             "results/{sample}/fastq/{sample}_2.fastq.gz"
#         log:
#             "results/{sample}/logs/fastq_links.log"
#         run:

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

            module load samtools/1.9
            samtools index \
              results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam

            cp results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam \
              temp
            samtools index \
              temp/{wildcards.sample}.Aligned.sortedByCoord.out.bam
            rm temp/*gtf*; cp {params.featureFile} temp; gunzip temp/*gtf.gz \
              >/dev/null 2>&1

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
            "results/{sample}/{sample}.salmon/{sample}.aligned.sorted.bam",
            "results/{sample}/{sample}.salmon/{sample}.aligned.sorted.bam.bai"
        params:
            outDir = "results/{sample}/{sample}.salmon",
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

#if config["countMatrix"]:
os.makedirs("results/counts", exist_ok=True)
rule matrixGeneration:
    input:
        expand("results/{sample}/{sample}.salmon/quant.sf", sample=samples)
    output:
        "results/counts/" + analysisName + ".counts.txt",
        "results/counts/" + analysisName + ".tpm.txt",
        "results/counts/" + analysisName + ".txi.rds"
    params:
        configFn = configFilename,
        tx2gene = config[genomeBuild][config["countGeneSymbols"]],
        tpm = "results/counts/" + analysisName + ".tpm.txt",
        counts = "results/counts/" + analysisName + ".counts.txt",
        txi = "results/counts/" + analysisName + ".txi.rds"
    log:
        "results/counts/matrixGeneration.log"
    shell:
        """
        mkdir -p results/counts
        module load r/3.6.0
        Rscript /proj/fureylab/bin/countMatrixGeneration.R {params.configFn} \
          {params.tx2gene} --tpm {params.tpm} --counts {params.counts} \
          --txi {params.txi} --pipeline > {log}
        """

# rule deconvolution:
#     input:
#     output:
#     log:
#     shell:

rule bamCoverage:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/star/{sample}.RPKM.bw"
    params:
        binSize = config["bamCoverage"]["binSize"],
        normalizeUsing = config["bamCoverage"]["normalizeUsing"],
        ignoreForNormalization = config["bamCoverage"]["ignoreForNormalization"]
    log:
        "results/{sample}/logs/bamCoverage.log"
    shell:
        """
        module load deeptools/3.2.0
        bamCoverage --ignoreForNormalization {params.ignoreForNormalization} \
          -o {output} --binSize {params.binSize} \
          --normalizeUsing {params.normalizeUsing} --bam {input} > {log}
        """

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

os.makedirs("results/multiqc_raw", exist_ok=True)
rule multiqc_raw:
    input:
        expand("results/{sample}/QC/fastqc/{sample}_{pair}_fastqc.zip",
            sample=samples, pair=["1","2"])
    output:
        "results/multiqc_raw/multiqc_raw.html"
    log:
        "results/multiqc_raw/multiqc_raw.log"
    shell:
        """
        module load multiqc/1.7
        multiqc -n {output} results/*/QC/fastqc > {log}
        """

rule fastQscreen:
    input:
        "results/{sample}/fastq/{sample}_{pair}.fastq.gz"
    output:
        "results/{sample}/QC/fastQscreen/{sample}_{pair}_screen.txt"
    params:
        aligner = config["fastQscreen"]["aligner"],
        outDir = "results/{sample}/QC/fastQscreen",
        subset = config["fastQscreen"]["subset"],
        conf = config["fastQscreen"]["conf"]
    log:
        "results/{sample}/logs/fastQscreen_{pair}.log"
    shell:
        """
        mkdir -p results/{wildcards.sample}/QC/fastQscreen
        module load bowtie2/2.3.4
        module load samtools/1.9
        perl /proj/fureylab/bin/fastq_screen_0.14/fastq_screen \
          --conf {params.conf} --aligner {params.aligner} \
          --outdir {params.outDir} --subset {params.subset} {input} > {log}
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
        module load r/3.6.0
        module load rseqc/3.0.1
        mkdir -p results/{wildcards.sample}/QC/rseqc
        tin.py -i {input} -r {params.model} > {log}
        junction_saturation.py -i {input} -r {params.model} \
          -o {wildcards.sample} >> {log}
        mv {wildcards.sample}* results/{wildcards.sample}/QC/rseqc
        """

rule QM_bamqc:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/qualimap/{sample}.bamqc/qualimapReport.html",
        "results/{sample}/QC/qualimap/{sample}.bamqc/{sample}.genomeCoverage.txt"
    params:
        outDir = "results/{sample}/QC/qualimap/{sample}.bamqc",
        outCoverage = "results/{sample}/QC/qualimap/{sample}.bamqc/{sample}.genomeCoverage.txt",
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
          -oc {params.outCoverage} -outformat HTML \
          --sequencing-protocol {params.seqProtocol} \
          --feature-file {params.featureFile} --genome-gc-distr {params.genomeGC} \
          --java-mem-size={params.javaMemSize} {params.otherFlags} > {log}
        """

rule QM_rnaseq:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.computedCounts.txt",
        "results/{sample}/QC/qualimap/{sample}.rnaseq/qualimapReport.html"
    params:
        outDir = "results/{sample}/QC/qualimap/{sample}.rnaseq",
        outCounts = "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.computedCounts.txt",
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
        -oc {params.outCounts} --sequencing-protocol {params.seqProtocol} \
        -gtf {params.featureFile} --java-mem-size={params.otherFlags} \
        {params.otherFlags} > {log}
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

#############################
#### CLEAN UP AND MOVING ####
#############################

# For the purposes of multiqc, some directories had to be sample specific
# This rule fixes those names.
rule name_clean:
    input:
        "results/multiqc/multiqc.html",
        "results/{sample}/{sample}.salmon/{sample}.aligned.sorted.bam.bai",
        "results/{sample}/QC/qualimap/{sample}.rnaseq/qualimapReport.html",
        "results/{sample}/QC/qualimap/{sample}.bamqc/qualimapReport.html",
        "results/{sample}/star/{sample}.RPKM.bw"
    output:
        directory("results/{sample}/salmon"),
        directory("results/{sample}/QC/qualimap/rnaseq"),
        directory("results/{sample}/QC/qualimap/bamqc"),
        touch("temp/{sample}/name_clean_complete.flag")
    log:
        "results/{sample}/logs/cleanup.log"
    shell:
        """
        mv results/{wildcards.sample}/{wildcards.sample}.salmon \
          results/{wildcards.sample}/salmon
        mv results/{wildcards.sample}/QC/qualimap/{wildcards.sample}.rnaseq \
          results/{wildcards.sample}/QC/qualimap/rnaseq
        mv results/{wildcards.sample}/QC/qualimap/{wildcards.sample}.bamqc \
          results/{wildcards.sample}/QC/qualimap/bamqc
        cp project_config.yaml results/{wildcards.sample}
        """

# if config["moveOutFiles"]:
#
#     rule moveOut:
#         input:
#             "temp/{sample}/name_clean_complete.flag"
#         output:
#             expand("{path}/snakemakeRNA/QC/qualimap/bamqc/qualimapReport.html",
#                 , path=samplePaths)
#         log:
#             "{path}/snakemakeRNA/moveOut.log"
#         shell:
#             """
#             mkdir results/{wildcards.sample}/snakemakeRNA
#             mv results/{wildcards.sample}/* results/{wildcards.sample}/snakemakeRNA
#             cp -rf results/{wildcards.sample}/snakemakeRNA \
#               {wildcards.path}/snakemakeRNA
#             """
