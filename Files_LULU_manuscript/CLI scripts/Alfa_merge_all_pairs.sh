#!/bin/bash
#### Description: merging of all paired MiSeq libraries in directory
#### forward read library must contain R1 in their name and reverse read must contain R2. furthermore they need to be gzipped and end with .gz
#### no merged pairs will be outputtet in separate files
#### Written by: Tobias Guldberg Frøslev - tobiasgf@bio.ku.dk  on 12-2016

VSEARCH=$(which vsearch)
for i in *_R1*.gz; do
 echo "\n$i" >> $(basename ${i%%_L00*})_log.txt;
 "${VSEARCH}" --fastq_mergepairs $i \
         --reverse ${i/_R1/_R2} \
         --fastqout $(basename ${i%%_L00*}).merged.fastq \
         --fastq_allowmergestagger \
         --fastqout_notmerged_fwd $(basename ${i%%_L00*}).R1.not_merged.fastq \
         --fastqout_notmerged_rev $(basename ${i%%_L00*}).R2.not_merged.fastq \
         --threads 32 2>> $(basename ${i%%_L00*})_log.txt
done
