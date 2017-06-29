#!/bin/bash
#### Description: 0% (100%) clustering of reads for processing with dbotu3

#### Written by: Tobias Guldberg Frøslev - tobiasgf@bio.ku.dk  on 6-2017

#SETTING VARIABLES
VSEARCH=$(which vsearch)
THREADS=40
LOW_ABUND_SAMPLE=1
LOW_ABUND_GLOBAL=1
PROJ="VSEARCH_100"
TMP_FASTA1=$(mktemp)
TMP_FASTA2=$(mktemp)
TMP_FASTA3=$(mktemp)
TMP_FASTA4=$(mktemp)
TMP_FASTA5=$(mktemp)
TMP_FASTA7=$(mktemp)
LOG="${PROJ}.log"

#rename fasta headers in single files
cat S[0-9][0-9][0-9].fas >> "${TMP_FASTA7}"
rename.pl
cat renamed*.fas > "${TMP_FASTA1}"
rm renamed*.fas

# Dereplicate (vsearch)
"${VSEARCH}" --derep_fulllength "${TMP_FASTA7}" \
             --sizein \
             --sizeout \
             --fasta_width 0 \
	     --minuniquesize 1 \
	     --sizein -sizeout \
             --output "${TMP_FASTA2}" 2> "$LOG"

# Sorting
"${VSEARCH}" --fasta_width 0 \
             --sortbysize "${TMP_FASTA2}" \
             --output "${TMP_FASTA4}" --sizein --sizeout 2>> "$LOG"

CLUSTER=1

  NAME=$PROJ

 #Clustering
  "${VSEARCH}" --cluster_size "${TMP_FASTA4}" \
	     --id ${CLUSTER} \
	     --sizein --sizeout \
             --dbmask none --qmask none \
	     --centroids "${TMP_FASTA5}" 2>> "$LOG"

  # Sorting & filtering
  "${VSEARCH}" --fasta_width 0 \
             --sortbysize "${TMP_FASTA5}" \
             --minsize ${LOW_ABUND_GLOBAL} \
             --output $NAME.centroids --sizein --sizeout 2>> "$LOG"

  sed 's/;size=.*//' -i $NAME.centroids

  #mapping of raw reads against OTU representatives
  "${VSEARCH}" --usearch_global "${TMP_FASTA1}" \
	     --db $NAME.centroids  \
	     --strand plus \
	     --id ${CLUSTER} \
	     --maxaccepts 0 \
             --dbmask none --qmask none \
	     --uc $NAME.uc 2>> "$LOG"

  #preparing for making table
  sed -i 's/\tS0/\tbarcodelabel=S0/g' $NAME.uc
  sed -i 's/\tS1/\tbarcodelabel=S1/g' $NAME.uc
  sed -i 's/\tS2/\tbarcodelabel=S2/g' $NAME.uc
  sed -i 's/\tS3/\tbarcodelabel=S3/g' $NAME.uc
  sed -i 's/\tS4/\tbarcodelabel=S4/g' $NAME.uc
  sed -i 's/\tS5/\tbarcodelabel=S5/g' $NAME.uc
  sed -i 's/\tS6/\tbarcodelabel=S6/g' $NAME.uc
  sed -i 's/\tS7/\tbarcodelabel=S7/g' $NAME.uc
  sed -i 's/\tS8/\tbarcodelabel=S8/g' $NAME.uc
  sed -i 's/\tS9/\tbarcodelabel=S9/g' $NAME.uc

  #make table
  python ~/bin/uc2otutab.py $NAME.uc > $NAME.otutable 2>> "$LOG"

  rm -f "${TMP_FASTQ5}"


rm -f "${TMP_FASTA1}" "${TMP_FASTQ2}" "${TMP_FASTQ4}" "${TMP_FASTQ6}" "${TMP_FASTQ7}"

#cd ..
