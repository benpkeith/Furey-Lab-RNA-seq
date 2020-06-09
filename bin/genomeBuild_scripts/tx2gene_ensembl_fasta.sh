#!/usr/bin/env bash
#Ben Keith
#Usage: bash tx2gene_CCstrains_gtf.sh <transcriptsFasta>
#
#NOTES: The input is the transcripts fasta file used to build 
# the gentrome file in salmon indexing. the headers of this file
# contain transcript and gene id information that we can used to 
# create this lookup file for tximport.

fasta_file="$1"
echo "$fasta_file"

output_file="${2:-"tx2gene.csv"}"
echo "$output_file"

echo "Generating '$(basename "$output_file")' from \
'$(basename "$fasta_file")'."

echo "TXID,GENEID" > $output_file

#In the cut line, 1,7 for ensembl tx ID/HGNC symbol and 1,4 for ensembl tx/gene id
gunzip -c "$fasta_file" \
    | grep '>' \
    | cut -d ' ' -f1,7 \
    | awk '!a[$0]++' \
    | tr -d '>' \
    | sed -E 's/ [a-z]+:/,/g' \
    >> "$output_file"

# Return the number of transcripts.
count=$(wc -l $output_file)
echo "${count} transcripts detected."
