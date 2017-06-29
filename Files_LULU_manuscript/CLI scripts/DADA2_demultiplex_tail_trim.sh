mkdir DADA2_AS
mkdir DADA2_SS

while read INPUT_R1 INPUT_R2 TAGS PRIMER_F PRIMER_R MIN_LENGTH ; do

CUR_PATH=$(pwd)
#mkdir tempdir
#export TMPDIR="${CUR_PATH}/tempdir/"

# Define binaries, temporary files and output files
CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH} -e 0"
CUTADAPT2="$(which cutadapt) -e 0"
VSEARCH=$(which vsearch)
C1_FASTQ=$(mktemp)
C2_FASTQ=$(mktemp)
TMP_FASTQ=$(mktemp)
MIN_F=$((${#PRIMER_F}))
MIN_R=$((${#PRIMER_R}))

REV_PRIMER_F="$(echo $PRIMER_F | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"
REV_PRIMER_R="$(echo $PRIMER_R | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"

rev="$(echo $primer | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"


while read TAG_NAME TAG_SEQ RTAG_SEQ; do
    LOG="DADA2_SS/${TAG_NAME}_R1.log"
    FINAL_FASTQ="DADA2_SS/${TAG_NAME}_R1.fastq"

    # Trim tags, forward & reverse primers (search normal and antisens)
        cat "${INPUT_R1}" | \
        ${CUTADAPT} -g "${TAG_SEQ}" -e 0 -O "${#TAG_SEQ}" - 2> "${LOG}" | \
        ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2>> "${LOG}" > "${TMP_FASTQ}"

	cat "${TMP_FASTQ}" | ${CUTADAPT2} -a "${REV_PRIMER_R}" - 2>> "${LOG}"  > "${FINAL_FASTQ}"

    LOG="DADA2_AS/${TAG_NAME}_R2.log"
    FINAL_FASTQ="DADA2_AS/${TAG_NAME}_R2.fastq"

    # Trim tags, forward & reverse primers (search normal and antisens)
        cat "${INPUT_R2}" | \
        ${CUTADAPT} -g "${TAG_SEQ}" -e 0 -O "${#TAG_SEQ}" - 2> "${LOG}" | \
        ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2>> "${LOG}" > "${TMP_FASTQ}"

        cat "${TMP_FASTQ}" | ${CUTADAPT2} -a "${REV_PRIMER_R}" - 2>> "${LOG}" > "${FINAL_FASTQ}"

    LOG="DADA2_SS/${TAG_NAME}_R2.log"
    FINAL_FASTQ="DADA2_SS/${TAG_NAME}_R2.fastq"

    # Trim tags, forward & reverse primers (search normal and antisens)
        cat "${INPUT_R2}" | \
        ${CUTADAPT} -g "${RTAG_SEQ}" -e 0 -O "${#RTAG_SEQ}" - 2>> "${LOG}" | \
        ${CUTADAPT} -g "${PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" > "${TMP_FASTQ}"

        cat "${TMP_FASTQ}" | ${CUTADAPT2} -a "${REV_PRIMER_F}" - 2>> "${LOG}" > "${FINAL_FASTQ}"

    LOG="DADA2_AS/${TAG_NAME}_R1.log"
    FINAL_FASTQ="DADA2_AS/${TAG_NAME}_R1.fastq"

    # Trim tags, forward & reverse primers (search normal and antisens)
        cat "${INPUT_R1}" | \
        ${CUTADAPT} -g "${RTAG_SEQ}" -e 0 -O "${#RTAG_SEQ}" - 2>> "${LOG}" | \
        ${CUTADAPT} -g "${PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" > "${TMP_FASTQ}"

        cat "${TMP_FASTQ}" | ${CUTADAPT2} -a "${REV_PRIMER_F}" - 2>> "${LOG}" > "${FINAL_FASTQ}"

done < "${TAGS}"

# Clean
rm -rf "${C1_FASTQ}" "${C2_FASTQ}" "${TMP_FASTQ}"

done < batchfileDADA2.list
