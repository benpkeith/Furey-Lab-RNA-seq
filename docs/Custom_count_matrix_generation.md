2020.06.19 - Ben Keith

# Custom Count Matrix Generation

Although the snakemake pipeline generates a count matrix that you can use for downstream analyses, there will be scenarios where you need to generate a count matrix from additional data or from data that is already stored in the Furey lab data directory. This serves as a simple manual to generate count matrices from samples that have already been processed by this snakemake pipeline and have been stored in permanent space.

### Setting up the config file

To generate this count matrix, you'll need to generate a config file that points to the samples you want a count matrix for. For this example, we'll call this config file _countMat.yaml_. This file will only need one section, **samples**, with paths to each sample. This config will look like the following:

```
samples:
  - /proj/fureylab/data/RNA-seq/human/ileum_tissue/nonIBD/SRR10104064/snakemakeRNA_hg38
  - /proj/fureylab/data/RNA-seq/human/ileum_tissue/nonIBD/SRR10104067/snakemakeRNA_hg38
  - /proj/fureylab/data/RNA-seq/human/ileum_tissue/CD_inflamed/SRR10104068/snakemakeRNA_hg38
  - /proj/fureylab/data/RNA-seq/human/ileum_tissue/CD_inflamed/SRR10104069/snakemakeRNA_hg38
```

For your count matrix, you just need to add the relevant paths for the samples of interest. **It is important that this path ends in the snakemake directory that you want to pull data from**. The script will look for the quantification files relative to this path.

### Running the script

Now that you have the config file, to run the script enter a command like the following:

```
module load r/3.6.0
Rscript /proj/fureylab/bin/countMatrixGeneration.R config.yaml \
/proj/fureylab/genomes/human/hg38_reference/RNA_annotation/tx2gene.GRCh38.ENSEMBL.csv \
--counts project.counts.txt \
--tpm project.tpm.txt \
--txi project.txi.rds
```

For more information about to structure this command check out the help page for this script by entering:

```
module load r/3.6.0
Rscript /proj/fureylab/bin/countMatrixGeneration.R --help
```

If you're moving from this to the RNA-seq tutorial (_RNA-seq_analysis.md_), you'll need the raw counts file ending "counts.txt".
