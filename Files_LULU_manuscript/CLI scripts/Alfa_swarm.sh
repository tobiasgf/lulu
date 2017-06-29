#!/bin/bash
#### Description: Using SWARM to define OTUs

#### Written by: Tobias Guldberg Frøslev - tobiasgf@bio.ku.dk  on 1-2016

VSEARCH=$(which vsearch)
SWARM=$(which SWARM2)
TMP_FASTA=$(mktemp --tmpdir=".")
mkdir -p SWARM_PIPE
PROJ="SWARM_"
for DVALUE in {3,5,7,10,13,15}; do
  FINAL_FASTA=$PROJ"$DVALUE.fas"
  cat "Processing "$FINAL_FASTA

  # Pool sequences
  cat S[0-9][0-9][0-9]*.fas > "${TMP_FASTA}"

  # Dereplicate (vsearch)
  "${VSEARCH}" --derep_fulllength "${TMP_FASTA}" \
             --sizein \
             --sizeout \
             --fasta_width 0 \
             --output "${FINAL_FASTA}" > /dev/null
  rm -f "${TMP_FASTA}"

  # Clustering
  THREADS=16
  TMP_REPRESENTATIVES=$(mktemp --tmpdir=".")
  SWARM2 \
    -d ${DVALUE} -t ${THREADS} -z \
    -i ${FINAL_FASTA/.fas/.struct} \
    -s ${FINAL_FASTA/.fas/.stats} \
    -w ${TMP_REPRESENTATIVES} \
    -o ${FINAL_FASTA/.fas/.swarms} < ${FINAL_FASTA}

  # Sort representatives
  "${VSEARCH}" --fasta_width 0 \
             --sortbysize ${TMP_REPRESENTATIVES} \
             --output ${FINAL_FASTA/.fas/.centroids}

  rm ${TMP_REPRESENTATIVES}

  # Chimera checking
  REPRESENTATIVES=${FINAL_FASTA/.fas/.centroids}
  UCHIME=${REPRESENTATIVES/.centroids/.uchime}
  "${VSEARCH}" --uchime_denovo "${REPRESENTATIVES}" \
             --uchimeout "${UCHIME}"

  # Build OTU table
  buildOTUtable_simple.sh ${FINAL_FASTA}

done

rm "${TMP_FASTA}"

mv $PROJ* SWARM_PIPE/
