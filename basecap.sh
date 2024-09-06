#!/bin/bash
# set -e

# <-.(`-')  (`-')  _  (`-').->(`-')  _           (`-')  _  _  (`-')
#  __( OO)  (OO ).-/  ( OO)_  ( OO).-/ _         (OO ).-/  \-.(OO )
# '-'---.\  / ,---.  (_)--\_)(,------. \-,-----. / ,---.   _.'    \
# | .-. (/  | \ /`.\ /    _ / |  .---'  |  .--./ | \ /`.\ (_...--''
# | '-' `.) '-'|_.' |\_..`--.(|  '--.  /_) (`-') '-'|_.' ||  |_.' |
# | /`'.  |(|  .-.  |.-._)   \|  .--'  ||  |OO )(|  .-.  ||  .___.'
# | '--'  / |  | |  |\       /|  `---.(_'  '--'\ |  | |  ||  |
# `------'  `--' `--' `-----' `------'   `-----' `--' `--'`--'

# Bash script for BaseCap workflow.
# Input directory containing ab1 chromatogram files, primer sequences, phred threshold 
# and gene name.

# Notes on usage #
if [ "$1" == "-h" ]; then
  echo "Usage:  ./baecap.sh  /path/to/ab1_dir/ primer_file  threshold  gene_name protein_reference"
  exit 0
fi

#  time ./basecap.sh cd63_reads/ test_primers.txt 20 cd63 cd63_prot.faa  

# args
ab1_dir="$1"
primers="$2"
threshold="$3"
gene_name="$4"
reference="$5"

# remove empty lines in primer file..
sed -i '/^$/d' "$primers"

#res1=$(date +%s.%N)
dir=$(pwd)

mkdir -p "$dir"/"$gene_name"

#call het bases (for homo/heterozyg info)
/usr/app/BaseCap/tracetuner_3.0.6beta/rel/Linux_64/ttuner -3730 -het -trim_threshold "$threshold" -id "$ab1_dir" -pd "$dir"/"$gene_name" -tabd "$dir"/"$gene_name" -Q

#call all bases
/usr/app/BaseCap/tracetuner_3.0.6beta/rel/Linux_64/ttuner -3730 -trim_threshold "$threshold" -id "$ab1_dir" -pd "$dir"/"$gene_name" -recalln -Q

# Part 1:Exonerate and BaseCap
{ cd "$dir"/"$gene_name" || exit 1

counter=1

for f in *.phd.1; do

    counter=$((counter + 1))

    {
        echo ">$counter"
        cut -f1 "$f" -d " " | awk '/BEGIN_DNA/{flag=1;next}/END_DNA/{flag=0}flag' | tr -d '[:space:]'
    } > "$f".exfa

    exonerate --model protein2genome --query "$dir"/"$reference" \
    --target  "$f".exfa --showtargetgff yes --showalignment no --showvulgar no --showquerygff no --bestn 1 \
    --ryo "%qi(%qab - %qae)\n%qas\n >%ti:%s(%tab - %tae)\n%tas\n >%ti:%s\n%tcs\n" >  "$f".exonerate

   if  grep  '\-\-\- START OF GFF DUMP \-\-\-' "$f".exonerate ; then
	echo "$f aligns to ref"
    else
        echo "$f Contaminant/no hit to reference"
	mkdir -p "contaminants"
	mv "$f"* "contaminants"
    fi;done

   python3 "$dir"/BaseCap.py ./ "$dir"/"$primers" "$threshold" >> "$gene_name".summary

# remove exfa files
rm -rf *.exfa

# create results directory
mkdir -p "results"

mv *.phd.1 "results"
mv *.fasta "results"
mv *.exonerate "results"
mv *.qual "results"

}

# Part 2: CAP3 merging and final cap.
# move to results dir for read merging and capping of merged reads

{ cd "$dir"/"$gene_name"/"results"

for k in *.fasta;do
	r_id="$(echo $k | cut -d "_" -f1)"

	# use .fa to prevent looping over fasta files
	# qual must be name ..FR.fa.qual for cap3 to pair with ..FR.fa
	awk '1' "$r_id"*F*.qual "$r_id"*R*.qual > "$r_id"_FR.fa.qual
	awk '1' "$r_id"*F*.fasta "$r_id"*R*.fasta > "$r_id"_FR.fa;done

# merge with cap3 using qual scores
# generates contig and new qual files
for l in *_FR.fa;do
	cap3 "$l" > "$gene_name".out;done

# remove the paired files
rm -rf *_FR.fa
rm -rf *_FR.fa.qual

# cap the low scoring bases in the merged contig seq
for m in *_FR*.contigs;do
	r2_id="$(echo $m | cut -d "_" -f1)"
	python3 "$dir"/capper.py "$threshold" "$r2_id"*.contigs.qual "$r_id"*.contigs "$r2_id";done

# awk and combine all final sequences
awk '1' *capped.merged.fasta > "$gene_name".all.final.fasta

mafft --quiet --auto --preservecase "$gene_name".all.final.fasta > "$gene_name"_aligned.basecapped.fa

}
