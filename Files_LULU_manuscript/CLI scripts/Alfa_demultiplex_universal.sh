#!/bin/bash
#### Description: Demultiplexing of reads for the VSEARCH, SWARM, CROP and DBOTU pipelines
#### Fastafiles will be rereplicated sample wise. Dereplicated fastafiles should be in PWD and be named S001.fas to S999.fas
#### script needs access to VSEARCH and CROP

#### Written by: Tobias Guldberg Frøslev - tobiasgf@bio.ku.dk  on 12-2016

while read INPUT TAGS PRIMER_F PRIMER_R MIN_LENGTH ; do

# Define binaries, temporary files and output files
CUTADAPT="$(which cutadapt) -e 0 --discard-untrimmed --minimum-length ${MIN_LENGTH}"
VSEARCH=$(which vsearch)
INPUT_REVCOMP=$(mktemp)
TMP_FASTQ=$(mktemp)
TMP_FASTQ2=$(mktemp)
TMP_FASTA=$(mktemp)
OUTPUT=$(mktemp)
QUALITY_FILE="${INPUT}.qual"

# Reverse complement fastq file
"${VSEARCH}" --quiet \
             --fastx_revcomp "${INPUT}" \
             --fastqout "${INPUT_REVCOMP}"

while read TAG_NAME TAG_SEQ RTAG_SEQ; do
    LOG="${TAG_NAME}.log"
    FINAL_FASTA="${TAG_NAME}.fas"
    F_CUT=$TAG_SEQ$PRIMER_F
    R_CUT=$PRIMER_R$RTAG_SEQ
    MIN_F=$((${#F_CUT}))
    MIN_R=$((${#R_CUT}))

    # Trim tags, forward & reverse primers (search normal and antisens)
    cat "${INPUT}" "${INPUT_REVCOMP}" | \
        ${CUTADAPT} -g "${F_CUT}" -O "${MIN_F}" - 2> "${LOG}" | \
        ${CUTADAPT} -a "${R_CUT}" -O "${MIN_R}" - 2>> "${LOG}" > "${TMP_FASTQ}"

    # Discard erroneous sequences and add expected error rates
    "${VSEARCH}" \
        --quiet \
        --fastq_filter "${TMP_FASTQ}" \
        --fastq_maxns 0 \
        --relabel_sha1 \
        --fastq_maxee_rate 0.002 \
        --minseqlength 10 \
        --eeout \
        --fastqout "${TMP_FASTQ2}" 2>> "${LOG}"

    # Convert fastq to fasta (discard sequences containing Ns)
    "${VSEARCH}" \
        --quiet \
        --fastq_filter "${TMP_FASTQ}" \
        --fastq_maxns 0 \
        --fastq_maxee_rate 0.002 \
        --minseqlength 10 \
        --fastaout "${TMP_FASTA}" 2>> "${LOG}"

    # Dereplicate at the study level (vsearch)
    "${VSEARCH}" \
        --quiet \
        --derep_fulllength "${TMP_FASTA}" \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --minseqlength 10 \
        --output "${FINAL_FASTA}" 2>> "${LOG}"

    # Discard quality lines, extract sha1, expected error rates and read length
    sed 'n;n;N;d' "${TMP_FASTQ2}" | \
        awk 'BEGIN {FS = "[;=]"}
             {if (/^@/) {printf "%s\t%s\t", $1, $3} else {print length($1)}}' | \
        tr -d "@" >> "${OUTPUT}"

done < "${TAGS}"

# Produce the final quality file
sort -k3,3n -k1,1d -k2,2n "${OUTPUT}" | \
    uniq --check-chars=40 > "${QUALITY_FILE}"

# Clean
rm -rf "${TMP_FASTQ}" "${TMP_FASTQ2}" "${TMP_FASTA}" "${INPUT_REVCOMP}" "${OUTPUT}"

done < batchfile.list

#pool and sort the qual files
cat *.qual > pooled.qual
sort -k3,3n -k1,1d -k2,2n pooled.qual | uniq --check-chars=40 > sorted.qual
rm pooled.qual

