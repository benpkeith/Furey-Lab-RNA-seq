# Furey Lab RNA-seq

This analysis utilises a snakemake pipeline to process RNA-seq data. The only prerequisite to running the pipeline is to load Python 3 using:

```
module load python/3.6.6
```

config.yaml contain sample IDs for fetching from sra, and sample information that is used by R during downstream analysis. An updated pointer to the rsem reference directory should added before running the pipeline. 

config.json is a SLURM configuration file. This file contains general rule submission parameters that can be customised for specific rules. For example, the rsem step of the pipeline requires more that the default 15GB of memory, and so is bumped up to 100GB for the submission of these jobs.

### Running the pipeline

The pipeline is designed to run on a system using the SLURM job scheduler on UNC's longleaf HPC. The rulegraph for the pipeline can be generated using the command:

```
snakemake --rulegraph | dot -Tpng > Snakemake_graph.png
```

To run the pipeline completely, use:

```
sbatch -o snakePipe.log -t 1-0 --wrap='snakemake -p -j 100 
--cluster-config config.json --cluster "sbatch -t {cluster.time} --mem {cluster.mem}"'
```

This command will kick the pipeline off, submitting the pipeline as its own SLURM job and each individual step as a SLURM job. 'logs', 'raw', 'results', and 'temp' directories will be created, and final files of interest will be found in the 'results' directory.
