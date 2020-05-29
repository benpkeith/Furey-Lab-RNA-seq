#!/usr/bin/env bash

gtf_file="$1"
echo "$gtf_file"

house_keeping_file="$2"
echo "$house_keeping_file"

output_file="${3:-"houseKeepingGeneModel.bed12"}"
echo "$output_file"

module load ucsctools/320
gtfToGenePred "$gtf_file" temp.genePred
genePredToBed temp.genePred temp.bed12
module load bedtools/2.29
intersectBed -a temp.bed12 -b "$house_keeping_file" | sort -R | head -7500 \
  > temp.subset.bed12
sortBed -i temp.subset.bed12 > "$output_file"
