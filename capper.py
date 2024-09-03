#!/usr/bin/env python3.7

import sys

threshold=sys.argv[1]
qual_in=sys.argv[2]
fasta_in=sys.argv[3]
seq_id=sys.argv[4]


with open(fasta_in, "r") as fasta:
	f=fasta.read()
	merged_seq=list("".join(f.split("\n")[1:]).lower())

with open(qual_in, "r") as scores:
	s=scores.read()
	merged_scores=("".join(s.split("\n")[1:])).split()

# use enumerate to get each occurance of the score.(better than dict or single index)
indexed_scores=""

for k,j in enumerate(merged_scores):
    if int(j) <= int(threshold):
        indexed_scores += str(k)+ " "

# if no low scores pass
if indexed_scores == []:
    final_seq = "".join(">" + str(seq_id) +"\n" + "".join(merged_seq))
    final=open(str(seq_id) + "_capped.merged.fasta", "w")
    final.write(final_seq)
else:
    capped_merged = merged_seq.copy()
    indexed_scores=indexed_scores.split()

    for index in indexed_scores:
        each=int(index)
        capped_merged[each] = capped_merged[each].upper()

    final_seq = "".join(">" + str(seq_id) +"\n" + "".join(capped_merged))
    final=open(str(seq_id) + "_capped.merged.fasta", "w")
    final.write(final_seq)