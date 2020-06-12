# Ben Keith
# Last updated 2020.06.09
# Furey Lab Pipeline 2020
# Snakemake 1.1

########################
#### Initial set up ####
########################

import os
import glob
import re
import sys
import os.path
from os import path
from datetime import date

configfile: "project_config.yaml"
configFilename = "project_config.yaml"
today = date.today()

#genome variables
analysisName = config["analysis"]["projectName"]
organism = config["analysis"]["organism"]
genomeBuild = config["analysis"]["genomeBuild"]

#index variables
# This block of variables handles differences in the processing that will occur
# due to sequencing read length. Rather than manually changing the paths in the
# index pointer section of the config, this will change the index paths based on
# the read length in the config.
# NOTE: This assumes that the index at the chosen read length is in correct
# genomes folder (see README for more instructions)
# NOTE: Star indexing requires a specific overlap parameter which is "readLength - 1".
# Salmon only requires a "k" parameter related to read legnth, which remains consistent
# over read lengths of 75. Therefore, we may have more read-specific star index files
# compared to Salmon.
readLength = config["analysis"]["readLength"]

starOverhang = readLength - 1
trimLengths = [50,75,150]
cutadaptTrimLength = min(trimLengths, key=lambda x:abs(x-readLength))

starIndex = "star_" + str(readLength) + "bp"
salmonIndex = "salmon_" + str(readLength) + "bp" if readLength < 75 else "salmon_75bp+"
config[genomeBuild]["starIndex"] = config[genomeBuild]["starIndex"] + "/" + starIndex
config[genomeBuild]["salmonIndex"] = config[genomeBuild]["salmonIndex"] + "/" + salmonIndex

# Rather than setting rules for file handling and wrangling, I thought it
# much easier to use a bit of python before executing the pipline.
# These lines handle setting up fastq directories for either SRA or local files
# inputs.
# The config["moveOutFiles"] conditional loop handles the moving of processed
# files to the directory specificed in the configuration file. As mentioned in
# the README and config file, this flag should only be set AFTER the pipline
# has finished. Nothing will be copied if this is set as the pipeline runs, but
# for safety it is best to move the files after the pipeline has successfully
# finished.
if config["useSRA"]:
    samples = config["samples"]
    [os.makedirs("results/" + str(sample) + "/logs", exist_ok=True) for sample in samples]
else:
    samplePaths = config["samples"]
    samplePaths = [re.sub(r"\W+$", "", i) for i in samplePaths]
    samples = [i.rsplit('/', 1)[1] for i in samplePaths]

    [os.makedirs("results/" + str(sample) + "/logs", exist_ok=True) for sample in samples]
    [os.makedirs("results/" + str(sample) + "/fastq", exist_ok=True) for sample in samples]

    for path in samplePaths:
        samp = str(path.rsplit('/', 1)[1])

        if config["moveOutFiles"] and not config["useSRA"]:
            os.makedirs(str(path) + "/snakemakeRNA_" + str(genomeBuild), exist_ok=True)
            cmd = "cp -rf results/" + str(samp) + "/* " + str(path) + \
              "/snakemakeRNA_" + str(genomeBuild)
            os.system(cmd)
            print("moved results outputs for " + str(samp) + " to " + str(path))
            os.makedirs(str(path) + "/snakemakeRNA_" + str(genomeBuild), exist_ok=True)

        for file in glob.glob(str(path + "/fastq/*gz")):
            if "_R1_" in file or re.search("_1\.f*q\.gz$", file):
                cmd = "ln -s " + str(file) + " results/" + \
                  samp + "/fastq/" + samp + "_1.fastq.gz" +  " >/dev/null 2>&1"
                os.system(cmd)
            else:
                cmd = "ln -s " + str(file) + " results/" + \
                  samp + "/fastq/" + samp + "_2.fastq.gz" +  " >/dev/null 2>&1"
                os.system(cmd)

if config["moveOutFiles"]:
    projectDir = config["projectDir"] + config["analysis"]["projectName"] \
      + "_" + today.strftime("%Y_%m_%d")
    os.makedirs(projectDir, exist_ok=True)
    cmd = "cp -rf results/multiqc " + projectDir +\
      "; cp -rf results/counts " + projectDir
    os.system(cmd)
    print("Moving project files to " + projectDir)

    print("Files moved! Exiting...")
    print("The SystemExit message below this is normal!")
    sys.exit()

#setting up for the feature file for QC
if not glob.glob('temp/*gtf'):
    cmd = "mkdir temp; cp " + config[genomeBuild]["featureFile"] + " temp; gunzip temp/*gtf.gz"
    os.system(cmd)

################################################################################
################################### Pipeline ###################################
################################################################################

####################
#### all target ####
####################

# Target to run whole workflow:
rule all:
    input:
        "results/multiqc/" + config["analysis"]["projectName"] + "_multiqc.html",
        "results/counts/" + analysisName + ".txi.rds",
        expand("results/{sample}/logs/cleanup.log", sample=samples)

############################
#### "Checkpoint rules" ####
############################
# These sets of rules are usually needed prior to a job that is executed only
# once. For most rules in the pipeline, a separate job is submitted for each
# sample. For example, Jobs like multiqc requires every prior QC job to be
# finished. Using a rule like below can output a single file that can be used
# as the dependency. I've found this to be a bit safer in testing that Using
# the individual samples as the dependecy.

# StarSalmon flag
rule alignmentQuantification:
    input:
        "results/{sample}/{sample}.salmon/quant.sf",
        "results/{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        "results/{sample}/{sample}.salmon/{sample}.aligned.sorted.bam"
    output:
        touch("temp/starSalmon_run.flag")

rule QC:
    input:
        expand("results/{sample}/QC/rseqc/{sample}.Aligned.sortedByCoord.out.summary.txt",
            sample=samples),
        expand("results/{sample}/QC/rseqc/{sample}.junctionSaturation_plot.pdf",
            sample=samples),
        expand("results/{sample}/QC/qualimap/{sample}.rnaseq/qualimapReport.html",
            sample=samples),
        expand("results/{sample}/QC/qualimap/{sample}.bamqc/qualimapReport.html",
            sample=samples),
        expand("results/{sample}/{sample}.salmon/quant.sf",
            sample=samples),
        expand("results/{sample}/QC/fastq_screen/{sample}_1_screen.txt",
            sample=samples),
        "results/multiqc_raw/" + config["analysis"]["projectName"] + "_multiqc_raw.html"
    output:
        touch("temp/QC_complete.flag")

#################################################
####### useSRA: Fetch and split SRA files #######
#### else: Create links to local fastq files ####
#################################################

if config["useSRA"]:
    # fetching .sra files based on IDs provided in the config file
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

    # Splitting the .sra file into .fastq files and zipping them up.
    rule fastq_dump:
        input:
            "temp/{sample}/fastq/{sample}.sra"
        output:
            "results/{sample}/fastq/{sample}_1.fastq.gz"
        params:
            fastq_dir = "results/{sample}/fastq/",
            other_flags = "--split-files"
        log:
            "results/{sample}/logs/fastq_dump.log"
        run:
            if config["end"] == "paired":
                shell("""
                module load sratoolkit/2.10.1
                fastq-dump {params.other_flags} -L 5 -O {params.fastq_dir} {input} \
                  > {log}
                gzip {params.fastq_dir}/{wildcards.sample}_1.fastq
                gzip {params.fastq_dir}/{wildcards.sample}_2.fastq
                """)
            if config["end"] == "single":
                shell("""
                module load sratoolkit/2.10.1
                fastq-dump {params.other_flags} -L 5 -O {params.fastq_dir} {input} \
                  > {log}
                gzip {params.fastq_dir}/{wildcards.sample}_1.fastq
                """)

##################################
#### Preprocessing - cutadapt ####
##################################

# Remove adaptors.
rule cutadapt:
    input:
        "results/{sample}/fastq/{sample}_1.fastq.gz"
    output:
        "temp/{sample}/fastq/{sample}_1.fastq.trimmed.gz"
    params:
        a = config["cutadapt"]["a"],
        A = config["cutadapt"]["A"],
        qualityCutoff = config["cutadapt"]["qualityCutoff"],
        minimumLength = config["cutadapt"]["minimumLength"],
        threads = config["cutadapt"]["threads"],
        length = cutadaptTrimLength,
        basename = "{sample}/fastq/{sample}"
    log:
        "results/{sample}/logs/cutadapt.log"
    run:
        if config["adapterTrimming"]:
            if config["end"] == "paired":
                shell("""
                module load cutadapt/2.9
                cutadapt -a {params.a} -A {params.A} \
                  --quality-cutoff {params.qualityCutoff} \
                  --cores {params.threads} \
                  --length {params.length} \
                  --minimum-length {params.minimumLength} \
                  -o temp/{params.basename}_1.fastq.trimmed.gz \
                  -p temp/{params.basename}_2.fastq.trimmed.gz \
                  results/{params.basename}_1.fastq.gz \
                  results/{params.basename}_2.fastq.gz \
                  > {log}
                """)
            if config["end"] == "single":
                shell("""
                module load cutadapt/2.9
                cutadapt -a {params.a} \
                  --quality-cutoff {params.qualityCutoff} \
                  --cores {params.threads} \
                  --length {params.length} \
                  --minimum-length {params.minimumLength} \
                  -o temp/{params.basename}_1.fastq.trimmed.gz \
                  results/{params.basename}_1.fastq.gz \
                  > {log}
                """)
        if not config["adapterTrimming"]:
            if config["end"] == "paired":
                shell("""
                cp results/{params.basename}_1.fastq.gz \
                  temp/{params.basename}_1.fastq.trimmed.gz
                cp results/{params.basename}_2.fastq.gz \
                  temp/{params.basename}_2.fastq.trimmed.gz
                """)
            if config["end"] == "single":
                shell("""
                cp results/{params.basename}_1.fastq.gz \
                  temp/{params.basename}_1.fastq.trimmed.gz
                """)

##################################################
#### Alignment and quantification - STAR/rsem ####
##################################################
# NOTE: Not currently avaible, but can be added based on commands used in
# previous RNA-seq pipeline
# If this is incorporated in the future, you should be able to use the
# star indexes that have already been generated.

# if config["quantification"] == "rsem":
#     rule rsem:
#         input:
#             ""
#         output:
#             ""
#         params:
#             param1 = ""
#         log:
#             ""
#         shell:
#             """
#
#             """
#
#     rule rsem_bam_sort:
#         input:
#             ""
#         output:
#             ""
#         log:
#             ""
#         shell:
#             """
#
#             """
#
#     rule rsem_bam_index:
#         input:
#             ""
#         output:
#             ""
#         log:
#             ""
#         shell:
#             """
#
#             """

####################################################
#### Alignment and quantification - STAR/salmon ####
####################################################

if config["quantification"] == "salmon":

    # Aligning adaptor trimmed reads using star. This rule performs a number of
    # steps:
    # 1. Perform star alignment
    # 2. Indexing of output bam files
    # 3. Copying bam to temp dir (necessary for proper multiqc processing) AND
    #    reindexing this file to prevent "old index" warning messages
    # 4. Copying of feature file (.gtf) to temp directory and unzip for QC
    # 5. Moving of star log files to the proper directory.
    rule star:
        input:
            "temp/{sample}/fastq/{sample}_1.fastq.trimmed.gz"
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
            sjdbOverhang = starOverhang,
            outFileNamePrefix = "results/{sample}/star/{sample}.",
            featureFile = config[genomeBuild]["featureFile"],
            basename = "{sample}/fastq/{sample}"
        log:
            "results/{sample}/logs/star/star.log"
        run:
            if config["end"] == "paired":
                shell("""
                mkdir -p results/{wildcards.sample}/star/
                module load star/2.7.3a
                star --runThreadN {params.threads} --sjdbOverhang {params.sjdbOverhang} \
                  --genomeDir {params.genomeDir} --readFilesCommand {params.readFilesCommand} \
                  --outSAMtype {params.outSAMtype} --outSAMunmapped Within \
                  --quantMode {params.quantMode} \
                  --outFileNamePrefix {params.outFileNamePrefix} \
                  --readFilesIn temp/{params.basename}_1.fastq.trimmed.gz \
                  temp/{params.basename}_2.fastq.trimmed.gz \
                  > {log}

                module load samtools/1.9
                samtools index \
                  results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam

                cp results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam \
                  temp
                cp results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam.bai \
                  temp

                mv {params.outFileNamePrefix}Log.final.out \
                  results/{wildcards.sample}/logs/star/
                mv {params.outFileNamePrefix}Log.progress.out \
                  results/{wildcards.sample}/logs/star/
                mv {params.outFileNamePrefix}Log.out \
                  results/{wildcards.sample}/logs/star/
                """)
            if config["end"] == "single":
                shell("""
                mkdir -p results/{wildcards.sample}/star/
                module load star/2.7.3a
                star --runThreadN {params.threads} --sjdbOverhang {params.sjdbOverhang} \
                  --genomeDir {params.genomeDir} --readFilesCommand {params.readFilesCommand} \
                  --outSAMtype {params.outSAMtype} --outSAMunmapped Within \
                  --quantMode {params.quantMode} \
                  --outFileNamePrefix {params.outFileNamePrefix} \
                  --readFilesIn temp/{params.basename}_1.fastq.trimmed.gz \
                  > {log}

                module load samtools/1.9
                samtools index \
                  results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam

                cp results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam \
                  temp
                cp results/{wildcards.sample}/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam.bai \
                  temp

                mv {params.outFileNamePrefix}Log.final.out \
                  results/{wildcards.sample}/logs/star/
                mv {params.outFileNamePrefix}Log.progress.out \
                  results/{wildcards.sample}/logs/star/
                mv {params.outFileNamePrefix}Log.out \
                  results/{wildcards.sample}/logs/star/
                """)

    # Salmon quantification
    rule salmon:
        input:
            "temp/{sample}/fastq/{sample}_1.fastq.trimmed.gz"
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
            tmpDir = "temp/{sample}/salmon",
            basename = "{sample}/fastq/{sample}"
        log:
            "results/{sample}/logs/salmon.log"
        run:
            if config["end"] == "paired":
                shell("""
                mkdir -p {params.outDir}
                mkdir -p {params.tmpDir}
                module load salmon/1.2.1
                salmon quant --libType {params.libType} {params.otherFlags} \
                  --numBootstraps={params.numBootstraps} --threads {params.threads} \
                  --writeMappings={params.tmpDir}/{wildcards.sample}.aligned.sam \
                  -i {params.index} -o {params.outDir} \
                  -1 temp/{params.basename}_1.fastq.trimmed.gz \
                  -2 temp/{params.basename}_2.fastq.trimmed.gz \
                  > {log}
                mv {params.outDir}/logs/salmon_quant.log \
                  results/{wildcards.sample}/logs
                rm -r {params.outDir}/logs
                """)
            if config["end"] == "single":
                shell("""
                mkdir -p {params.outDir}
                mkdir -p {params.tmpDir}
                module load salmon/1.2.1
                salmon quant --libType {params.libType} {params.otherFlags} \
                  --numBootstraps={params.numBootstraps} --threads {params.threads} \
                  --writeMappings={params.tmpDir}/{wildcards.sample}.aligned.sam \
                  -i {params.index} -o {params.outDir} \
                  -r temp/{params.basename}_1.fastq.trimmed.gz \
                  > {log}
                mv {params.outDir}/logs/salmon_quant.log \
                  results/{wildcards.sample}/logs
                rm -r {params.outDir}/logs
                """)

    # Conversion of salmon output sam file to bam, followed by indexing.
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

# dir needed before execution of the rule
os.makedirs("results/counts", exist_ok=True)
# Runs an R script to generate count matrices for downstream processing.
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

# NOTE: deconvolution will have to be added after samples are reprocessed
# The samples that are used in decolvolution reqire a little bit of
# proprocessing before they can be used for new samples.
# rule deconvolution:
#     input:
#     output:
#     log:
#     shell:

# Deeptools bamCoverage is used to create a bigWig file for RNA-seq
# visualisation on the genome browser.
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
#### Quality control #### End-specific QC
#########################

if config["end"] == "paired":
    # Standard fastqc.
    rule fastqc:
        input:
            expand("results/{sample}/fastq/{sample}_{pair}.fastq.gz",
                  sample=samples, pair=["1","2"])
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

    # dir needed before execution of the rule
    os.makedirs("results/multiqc_raw", exist_ok=True)
    # This multiqc using just the fastqc results. This serves as a way to quickly
    # view a QC report as the pipeline is still processing.
    rule multiqc_raw:
        input:
            expand("results/{sample}/QC/fastqc/{sample}_{pair}_fastqc.zip",
                sample=samples, pair=["1","2"])
        output:
            "results/multiqc_raw/" +config["analysis"]["projectName"] +"_multiqc_raw.html"
        log:
            "results/multiqc_raw/multiqc_raw.log"
        shell:
            """
            module load multiqc/1.7
            multiqc -n {output} results/*/QC/fastqc > {log}
            """

    # fastq_screen can be used to identify potential contamination within a sample.
    rule fastq_screen:
        input:
            "results/{sample}/fastq/{sample}_{pair}.fastq.gz"
        output:
            "results/{sample}/QC/fastq_screen/{sample}_{pair}_screen.txt"
        params:
            aligner = config["fastq_screen"]["aligner"],
            outDir = "results/{sample}/QC/fastq_screen",
            subset = config["fastq_screen"]["subset"],
            conf = config["fastq_screen"]["conf"]
        log:
            "results/{sample}/logs/fastq_screen_{pair}.log"
        shell:
            """
            mkdir -p results/{wildcards.sample}/QC/fastq_screen
            module load bowtie2/2.3.4
            module load samtools/1.9
            perl /proj/fureylab/bin/fastq_screen_0.14/fastq_screen \
              --conf {params.conf} --aligner {params.aligner} \
              --outdir {params.outDir} --subset {params.subset} {input} > {log}
            """

if config["end"] == "single":
    # Standard fastqc.
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

    # dir needed before execution of the rule
    os.makedirs("results/multiqc_raw", exist_ok=True)
    # This multiqc using just the fastqc results. This serves as a way to quickly
    # view a QC report as the pipeline is still processing.
    rule multiqc_raw:
        input:
            expand("results/{sample}/QC/fastqc/{sample}_{pair}_fastqc.zip",
                sample=samples, pair=["1"])
        output:
            "results/multiqc_raw/" +config["analysis"]["projectName"] +"_multiqc_raw.html"
        log:
            "results/multiqc_raw/multiqc_raw.log"
        shell:
            """
            module load multiqc/1.7
            multiqc -n {output} results/*/QC/fastqc > {log}
            """

    # fastq_screen can be used to identify potential contamination within a sample.
    rule fastq_screen:
        input:
            "results/{sample}/fastq/{sample}_{pair}.fastq.gz"
        output:
            "results/{sample}/QC/fastq_screen/{sample}_{pair}_screen.txt"
        params:
            aligner = config["fastq_screen"]["aligner"],
            outDir = "results/{sample}/QC/fastq_screen",
            subset = config["fastq_screen"]["subset"],
            conf = config["fastq_screen"]["conf"]
        log:
            "results/{sample}/logs/fastq_screen_{pair}.log"
        shell:
            """
            mkdir -p results/{wildcards.sample}/QC/fastq_screen
            module load bowtie2/2.3.4
            module load samtools/1.9
            perl /proj/fureylab/bin/fastq_screen_0.14/fastq_screen \
              --conf {params.conf} --aligner {params.aligner} \
              --outdir {params.outDir} --subset {params.subset} {input} > {log}
            """

#########################
#### Quality control #### Post-alignment QC
#########################

# Running tin score portions of RSeQC.
rule rseqc_tin:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/rseqc/{sample}.Aligned.sortedByCoord.out.summary.txt"
    params:
        model = config[genomeBuild]["rseqcModel"]
    log:
        "results/{sample}/logs/rseqc_tin.log"
    shell:
        """
        touch temp/{wildcards.sample}.Aligned.sortedByCoord.out.bam.bai
        module load r/3.6.0
        module load rseqc/3.0.1
        mkdir -p results/{wildcards.sample}/QC/rseqc
        tin.py -i {input} -r {params.model} > {log}
        mv {wildcards.sample}.Aligned.sortedByCoord* results/{wildcards.sample}/QC/rseqc
        """

# Running junction saturation portion of RSeQC
rule rseqc_junctionSaturation:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/rseqc/{sample}.junctionSaturation_plot.pdf"
    params:
        model = config[genomeBuild]["rseqcModel"]
    log:
        "results/{sample}/logs/rseqc_junctionSaturation.log"
    shell:
        """
        touch temp/{wildcards.sample}.Aligned.sortedByCoord.out.bam.bai
        module load r/3.6.0
        module load rseqc/3.0.1
        mkdir -p results/{wildcards.sample}/QC/rseqc
        junction_saturation.py -i {input} -r {params.model} \
          -o {wildcards.sample} > {log}
        mv {wildcards.sample}.junction* results/{wildcards.sample}/QC/rseqc
        """

# Running the bamqc portion of qualimap
rule QM_bamqc:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/qualimap/{sample}.bamqc/qualimapReport.html",
        "temp/{sample}/bamqc/{sample}.genomeCoverage.txt"
    params:
        outDir = "results/{sample}/QC/qualimap/{sample}.bamqc",
        tmpDir = "temp/{sample}/bamqc",
        seqProtocol = config["qualimap"]["seqProtocol"],
        featureFile = "temp/*gtf",
        genomeGC = organism,
        javaMemSize = config["qualimap"]["javaMemSize"],
        otherFlags = config["qualimap"]["bamqcFlags"]
    log:
        "results/{sample}/logs/QM_bamqc.log"
    shell:
        """
        module load qualimap/2.2.1
        mkdir -p {params.outDir}
        qualimap bamqc -bam {input} -outdir {params.outDir} \
          -oc {params.tmpDir}/{wildcards.sample}.genomeCoverage.txt -outformat HTML \
          --sequencing-protocol {params.seqProtocol} \
          --feature-file {params.featureFile} --genome-gc-distr {params.genomeGC} \
          --java-mem-size={params.javaMemSize} {params.otherFlags} > {log}
        """

# Running the bamqc portion of qualimap
rule QM_rnaseq:
    input:
        "temp/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.computedCounts.txt",
        "results/{sample}/QC/qualimap/{sample}.rnaseq/qualimapReport.html"
    params:
        outDir = "results/{sample}/QC/qualimap/{sample}.rnaseq",
        outCounts = "results/{sample}/QC/qualimap/{sample}.rnaseq/{sample}.computedCounts.txt",
        seqProtocol = config["qualimap"]["seqProtocol"],
        featureFile = "temp/*gtf",
        javaMemSize = config["qualimap"]["javaMemSize"]
    log:
        "results/{sample}/logs/QM_rnaseq.log"
    run:
        if config["end"] == "paired":
            shell("""
            module load qualimap/2.2.1
            mkdir -p {params.outDir}
            qualimap rnaseq -bam {input} -outdir {params.outDir} -outformat HTML \
              -oc {params.outCounts} --sequencing-protocol {params.seqProtocol} \
              -gtf {params.featureFile} --java-mem-size={params.javaMemSize} \
              --paired > {log}
            """)
        if config["end"] == "single":
            shell("""
            module load qualimap/2.2.1
            mkdir -p {params.outDir}
            qualimap rnaseq -bam {input} -outdir {params.outDir} -outformat HTML \
              -oc {params.outCounts} --sequencing-protocol {params.seqProtocol} \
              -gtf {params.featureFile} --java-mem-size={params.javaMemSize} \
              > {log}
            """)

# dir needed before execution of the rule
os.makedirs("results/multiqc", exist_ok=True)
# multiqc is performed after all previous jobs are completed. This will produce
# a single report incorporating all QC reports for all samples in the run.
rule multiqc:
    input:
        "temp/QC_complete.flag"
    output:
        html = "results/multiqc/" +config["analysis"]["projectName"] +"_multiqc.html",
        data = directory("results/multiqc/" +config["analysis"]["projectName"] +"_multiqc_data")
    params:
        configFile = config[genomeBuild]["multiqcConfig"]
    log:
        "results/multiqc/multiqc.log"
    shell:
        """
        module load multiqc/1.7
        multiqc -n {output.html} -c {params.configFile} results/ > {log}
        """

#############################
#### CLEAN UP AND MOVING ####
#############################

# For the purposes of multiqc, some directories had to have sample specific names
# This rule fixes those names.
rule name_clean:
    input:
        "results/multiqc/" + config["analysis"]["projectName"] + "_multiqc.html",
        "results/{sample}/{sample}.salmon/{sample}.aligned.sorted.bam.bai",
        "results/{sample}/QC/qualimap/{sample}.rnaseq/qualimapReport.html",
        "results/{sample}/QC/qualimap/{sample}.bamqc/qualimapReport.html",
        "results/{sample}/star/{sample}.RPKM.bw"
    output:
        directory("results/{sample}/salmon"),
        directory("results/{sample}/QC/qualimap/rnaseq"),
        directory("results/{sample}/QC/qualimap/bamqc"),
        touch("temp/{sample}/name_clean_complete.flag")
    params:
        snakemakeDir = "results/{sample}/snakemakeRNA_" + str(genomeBuild)
    log:
        "results/{sample}/logs/cleanup.log"
    run:
        if config["useSRA"]:
            shell("""
            mv results/{wildcards.sample}/{wildcards.sample}.salmon \
              results/{wildcards.sample}/salmon
            mv results/{wildcards.sample}/QC/qualimap/{wildcards.sample}.rnaseq \
              results/{wildcards.sample}/QC/qualimap/rnaseq
            mv results/{wildcards.sample}/QC/qualimap/{wildcards.sample}.bamqc \
              results/{wildcards.sample}/QC/qualimap/bamqc
            cp project_config.yaml results/{wildcards.sample}
            mkdir -p {params.snakemakeDir}
            mv results/{wildcards.sample}/* {params.snakemakeDir} >/dev/null 2>&1
            """)
        else:
            shell("""
            mv results/{wildcards.sample}/{wildcards.sample}.salmon \
              results/{wildcards.sample}/salmon
            mv results/{wildcards.sample}/QC/qualimap/{wildcards.sample}.rnaseq \
              results/{wildcards.sample}/QC/qualimap/rnaseq
            mv results/{wildcards.sample}/QC/qualimap/{wildcards.sample}.bamqc \
              results/{wildcards.sample}/QC/qualimap/bamqc
            cp project_config.yaml results/{wildcards.sample}
            """)
