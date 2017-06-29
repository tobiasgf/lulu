#!/bin/bash
#### Description: make match lists for all centroid files in directory containing the term "plantcentroids".

#### Written by: Tobias Guldberg Frøslev - tobiasgf@bio.ku.dk  on 1-2017

for f in *plantcentroids; do
 cat $f > PM.fasta
 sed 's/;size=[0-9]*;//g' -i PM.fasta # remove size annotation if any
 makeblastdb -in PM.fasta -parse_seqids -dbtype nucl
 ~/bin/blastn -db PM.fasta -num_threads 50 -outfmt '6 std qseqid sseqid pident' -out out_blast.txt -qcov_hsp_perc 80 -perc_identity 84 -query $f
 awk -F" " '{print $1 "\t" $2 "\t" $3}' out_blast.txt > $f.matchlist
 rm PM*  out_blast.txt
done
