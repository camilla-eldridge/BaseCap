#!/bin/bash

# Basecall and trim
ab1_dir="$1"
primers="$2"
threshold="$3"
gene_name="$4"


# remove empty lines in primer file..
sed -i '/^$/d' "$primers"

#res1=$(date +%s.%N)
dir=$(pwd)

mkdir -p "$dir"/"$gene_name"

#call het bases (for homo/heterozyg info)
ttuner -3730 -het -trim_threshold "$threshold" -id "$ab1_dir" -pd "$dir"/"$gene_name" -tabd "$dir"/"$gene_name" -Q

#call all bases
ttuner -3730 -trim_threshold "$threshold" -id "$ab1_dir" -pd "$dir"/"$gene_name" -recalln -Q

cd "$dir"/"$gene_name"

python3 "$dir"/BaseCap.py "$dir"/"$gene_name" "$dir"/"$primers" "$threshold" > "$gene_name".summary

awk '1' *.fasta  > "$gene_name"_trimmedseq.fa

echo "DONE"
