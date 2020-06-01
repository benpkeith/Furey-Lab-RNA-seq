#!/usr/bin/env bash

fasta_file="$1"
echo "$fasta_file"

output_file="${2:-"tx2gene"}"
echo "$output_file"

echo "Generating '$(basename "$output_file")' from \
'$(basename "$fasta_file")'."

echo "TXID,GENEID" > $output_file.HGNC.csv

#In the cut line, 1,6 for ensembl tx ID/HGNC symbol
gunzip -c "$fasta_file" \
    | grep '>' \
    | cut -d '|' -f1,6 \
    | awk '!a[$0]++' \
    | tr -d '>' \
    | sed -E 's/\|/,/g' \
    >> "$output_file.HGNC.csv"

echo "TXID,GENEID" > $output_file.ENSEMBL.csv

#In the cut line, 1,2 for ensembl tx/gene id
gunzip -c "$fasta_file" \
    | grep '>' \
    | cut -d '|' -f1,2 \
    | awk '!a[$0]++' \
    | tr -d '>' \
    | sed -E 's/\|/,/g' \
    >> "$output_file.ENSEMBL.csv"


# Return the number of transcripts.
count=$(wc -l $output_file.ENSEMBL.csv)
echo "${count} transcripts detected."
