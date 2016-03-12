#!/bin/bash
# This script fixes the strandinfo in the IlluminaProbe gtf file, as this
# information was missing from the original gtf. The strand information is
# retrieved from a fasta file which contains the strand info for the probenames
# usage 1 is fileToRead = IlluminaMatcheRNAseqPRobes
# 2 is fasta containing strand information

# for each line in $1
while read x; do
    # retrieve the probename from the IlluminaProbe gff file.
    HJprobename=$(echo "$x" | sed -E "s/^(.)+HJ-ProbeName \"//g" | tr -d "\"")
    # retrieve the strandinfo from $2 which is located at the end of the
    # line of the fastafile. Grep the probename from the fasta file and only
    # retrieve the last plus or minus sign
    Strand=$(grep -E "HJ-ProbeName \"$HJprobename\"" $2 | grep -E -o "(\+|\-)$")
    # if strand information is available, return the gff, replacing the 7th
    # column with the strandinfo, and renaming the HJ (Harm Jan) Probename with
    # gene_id
    if [[ "$Strand" ]]; then
        echo "$x" | awk -v strand="$Strand" 'BEGIN {FS="\t" }{$7 = strand; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | sed "s/HJ-ProbeName/gene_id/"
    fi
done < $1
