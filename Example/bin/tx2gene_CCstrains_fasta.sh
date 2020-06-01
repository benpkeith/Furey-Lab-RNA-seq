#!/usr/bin/env bash

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
