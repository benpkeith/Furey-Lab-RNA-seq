#!/usr/bin/env bash
#Ben Keith
#Usage: bash tx2gene_CCstrains_gtf.sh <gtfFile>
#
#NOTES: The gtf file is the gtf file used for build genomic indexes.
# This script will use the gtf file to create a lookup file of transcript id
# to gene id conversion that can be used be the pipeline to generate the a gene
# level count matrix,. For normal gencode genomes, such as hg38 and mm10, the transcripts
# fastq file contains all this information. For the CC strains, we have to use the 
# gtf file.

gtf_file="$1"
echo "$gtf_file"

output_file="${2:-"tx2gene"}"
echo "$output_file"

echo "Generating '$(basename "$output_file")' from \
'$(basename "$gtf_file")'."

echo "TXID,GENEID" > $output_file.ENSEMBL.csv

#In the cut line, 1,2 for ensembl tx/gene id
gunzip -c "$gtf_file" \
    | grep "transcript_id" \
    | cut -f9,10 \
    | cut -d ';' -f1,2 \
    | cut -d '"' -f4,2 \
    | sed -E 's/\"/,/g' \
    | awk -F "," '{ print $2 "," $1 }' \
    | uniq \
    >> "$output_file.ENSEMBL.csv"


# Return the number of transcripts.
count=$(wc -l $output_file.ENSEMBL.csv)
echo "${count} transcripts detected."
