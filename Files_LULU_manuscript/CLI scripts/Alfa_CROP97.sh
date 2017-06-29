#!/bin/bash
#### Description: Using CROP to define centroids and using VSEARCH to map the reads back to centroids.
#### Fastafiles will be rereplicated sample wise. Dereplicated fastafiles should be in PWD and be named S001.fas to S999.fas
#### script needs access to VSEARCH and CROP

#### Written by: Tobias Guldberg Frøslev - tobiasgf@bio.ku.dk  on 12-2016

#SETTING VARIABLES
CROP=$(which CROPLinux)
VSEARCH=$(which vsearch)
THREADS=40
LOW_ABUND_SAMPLE=1
LOW_ABUND_GLOBAL=1
CLUSTER=0.97
PROJ="CROP97_PIPE/CROP_97"
TMP_FASTA1=$(mktemp)
TMP_FASTA2=$(mktemp)
TMP_FASTA3=$(mktemp)
TMP_FASTA4=$(mktemp)
TMP_FASTA5=$(mktemp)
TMP_FASTA7=$(mktemp)
LOG="${PROJ}.log"

mkdir -p CROP97_PIPE

#mkdir VSEARCH_PIPE

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

# Count uniques
UNIQ=$(grep -c ">" ${TMP_FASTA2})
BP=$(($UNIQ/50))
B_PARAM=${BP%.*}

MEAN_LEN=$(awk '{/>/&&++a||b+=length()}END{print b/a}' ${TMP_FASTA2})
RML=${MEAN_LEN%.*}
ZP=$((150000/$RML))
Z_PARAM=${ZP%.*}

echo $B_PARAM
echo $Z_PARAM

# Rereplicating
"${VSEARCH}" --rereplicate "${TMP_FASTA2}" \
             --fasta_width 0 \
             --output "${TMP_FASTA3}" 2>> "$LOG"

CROPLinux  -i "${TMP_FASTA3}" \
             -o CROPout97 \
             -b ${B_PARAM} \
             -e 2000 \
             -s \
             -m 20 \
             -z ${Z_PARAM} \
             -r 1

NAME=$PROJ

  #mapping of raw reads against OTU representatives
  "${VSEARCH}" --usearch_global "${TMP_FASTA1}" \
	     --db CROPout97.cluster.fasta  \
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
  mv CROPout97.cluster.fasta $NAME.centroids

  rm -f "${TMP_FASTQ5}"

done

rm -f "${TMP_FASTA1}" "${TMP_FASTQ2}" "${TMP_FASTQ3}" "${TMP_FASTQ4}" "${TMP_FASTQ6}" "${TMP_FASTQ7}"

#cd ..
