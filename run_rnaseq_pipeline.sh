#!/bin/bash
set -euo pipefail

# RNA-seq PIPELINE v1.4 (Symlink + sample_sheet.tsv)
# For: 2D vs 3D × 0.3/1.5 × Anoxia/Normoxia

PROJECT_DIR="/mnt/d/OxoMech"
RAW_FASTQ_DIR="${PROJECT_DIR}/raw_fastq"
REFERENCE_DIR="${PROJECT_DIR}/reference"
SAMPLE_SHEET="${PROJECT_DIR}/sample_sheet.tsv"

HISAT2_INDEX="${REFERENCE_DIR}/mouse_index"
GTF="${REFERENCE_DIR}/Mus_musculus.GRCm39.115.gtf"

FASTQC_RAW="${PROJECT_DIR}/fastqc_raw"
FASTQC_TRIM="${PROJECT_DIR}/fastqc_trimmed"
TRIMMED="${PROJECT_DIR}/trimmed"
MAPPED="${PROJECT_DIR}/mapped"
COUNTS_DIR="${PROJECT_DIR}/counts"
QC_DIR="${PROJECT_DIR}/qc"

LOG_FILE="${PROJECT_DIR}/pipeline.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log "=== Starting RNA-seq Pipeline ==="

for dir in "$RAW_FASTQ_DIR" "$REFERENCE_DIR"; do
    if [[ ! -d "$dir" ]]; then
        log "ERROR: Directory not found: $dir"
        exit 1
    fi
done

if [[ ! -f "$SAMPLE_SHEET" ]]; then
    log "ERROR: sample_sheet.tsv not found"
    exit 1
fi

if [[ ! -f "${HISAT2_INDEX}.1.ht2" ]]; then
    log "ERROR: HISAT2 index not found"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    log "ERROR: GTF file not found"
    exit 1
fi

log "All pre-flight checks passed."

source /home/sree/mambaforge/etc/profile.d/conda.sh
conda activate mouse-rnaseq

mkdir -p "$FASTQC_RAW" "$FASTQC_TRIM" "$TRIMMED" "$MAPPED" "$COUNTS_DIR" "$QC_DIR"

TOTAL_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | wc -l)
log "Total samples to process: $TOTAL_SAMPLES"

CURRENT=0

while IFS=$'\t' read -r sample dim stiff oxy rep condition r1 r2; do
    CURRENT=$((CURRENT + 1))
    log "=== Processing sample $CURRENT of $TOTAL_SAMPLES: $sample ==="

    fastqc -t 8 -o "$FASTQC_RAW" "$r1" "$r2" >> "$LOG_FILE" 2>&1

    cutadapt -m 20 -q 20,20 --trim-n -j 8 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "${TRIMMED}/${sample}_R1.fastq.gz" \
        -p "${TRIMMED}/${sample}_R2.fastq.gz" \
        "$r1" "$r2" >> "$LOG_FILE" 2>&1

    fastqc -t 8 -o "$FASTQC_TRIM" \
        "${TRIMMED}/${sample}_R1.fastq.gz" \
        "${TRIMMED}/${sample}_R2.fastq.gz" >> "$LOG_FILE" 2>&1

    hisat2 -p 8 --rna-strandness RF -x "$HISAT2_INDEX" \
        -1 "${TRIMMED}/${sample}_R1.fastq.gz" \
        -2 "${TRIMMED}/${sample}_R2.fastq.gz" \
        --dta | samtools sort -@ 8 -o "${MAPPED}/${sample}.bam" - >> "$LOG_FILE" 2>&1

    samtools index "${MAPPED}/${sample}.bam"

    log "Finished processing: $sample"

done < <(tail -n +2 "$SAMPLE_SHEET")

log "Generating MultiQC reports..."
multiqc "$FASTQC_RAW" -o "$QC_DIR" --filename multiqc_raw_report.html >> "$LOG_FILE" 2>&1
multiqc "$FASTQC_TRIM" -o "$QC_DIR" --filename multiqc_trimmed_report.html >> "$LOG_FILE" 2>&1

log "Running featureCounts..."
featureCounts -p --countReadPairs -B -C -T 8 -s 2 \
    -t exon -g gene_id \
    -a "$GTF" \
    -o "${COUNTS_DIR}/counts.txt" \
    "${MAPPED}"/*.bam >> "$LOG_FILE" 2>&1

log "=== PIPELINE COMPLETED SUCCESSFULLY ==="
log "Total samples processed : $TOTAL_SAMPLES"
log "Counts file             : ${COUNTS_DIR}/counts.txt"
log "MultiQC reports         : ${QC_DIR}/"
log "Log file                : $LOG_FILE"
log "Pipeline finished at    : $(date)"
