#!/bin/bash
set -euo pipefail

# =============================================================================
# REVISED DIAGNOSTIC SCRIPT v2
# For new naming convention: 2D_0_3_A_Rep1_ / 3D_1_5_N_Rep2_ etc.
# =============================================================================

echo "=== Mouse 2D vs 3D Transcriptome - Diagnostic Check v2 ==="
echo "Date: $(date)"
echo

# === CONFIGURATION ===
RAW_DIR="/mnt/d/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome/2D raw fastq files"
NEW_PROJECT="/mnt/d/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome"

echo "Scanning folder: $RAW_DIR"
echo

# Check folder exists
if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: Folder not found: $RAW_DIR"
    exit 1
fi

mkdir -p "$NEW_PROJECT"

# --- Find all fastq.gz files ---
mapfile -t ALL_FASTQ < <(find "$RAW_DIR" -type f -name "*.fastq.gz" | sort)
echo "Total .fastq.gz files found: ${#ALL_FASTQ[@]}"
echo

if (( ${#ALL_FASTQ[@]} == 0 )); then
    echo "ERROR: No fastq.gz files found in $RAW_DIR"
    exit 1
fi

# --- Separate R1 and R2 ---
R1_LIST=()
R2_LIST=()

for f in "${ALL_FASTQ[@]}"; do
    base=$(basename "$f")
    if [[ "$base" =~ _R1_ ]]; then
        R1_LIST+=("$f")
    elif [[ "$base" =~ _R2_ ]]; then
        R2_LIST+=("$f")
    fi
done

echo "R1 files: ${#R1_LIST[@]}"
echo "R2 files: ${#R2_LIST[@]}"
echo

# --- Pair R1 and R2 ---
declare -A PAIRED
MISSING_R2=()

for r1 in "${R1_LIST[@]}"; do
    r1_base=$(basename "$r1")
    # Replace _R1_ with _R2_ to find matching pair
    r2_base=${r1_base/_R1_/_R2_}
    
    match=""
    for r2 in "${R2_LIST[@]}"; do
        if [[ "$(basename "$r2")" == "$r2_base" ]]; then
            match="$r2"
            break
        fi
    done
    
    if [[ -z "$match" ]]; then
        MISSING_R2+=("$r1")
    else
        PAIRED["$r1"]="$match"
    fi
done

if (( ${#MISSING_R2[@]} > 0 )); then
    echo "WARNING: ${#MISSING_R2[@]} R1 files have no matching R2:"
    printf '  %s\n' "${MISSING_R2[@]}"
fi

echo "Successfully paired samples: ${#PAIRED[@]}"
echo

# --- Parse metadata and create sample sheet ---
SAMPLE_SHEET="$NEW_PROJECT/sample_sheet.tsv"
echo -e "sample\tdim\tstiffness\toxygen\trep\tcondition\tR1\tR2" > "$SAMPLE_SHEET"

declare -A CONDITION_COUNT

for r1 in "${!PAIRED[@]}"; do
    r2="${PAIRED[$r1]}"
    base=$(basename "$r1")
    
    # Extract parts from new naming: 2D_0_3_A_Rep1_...
    if [[ "$base" =~ ^(2D|3D)_([0-9]_[0-9])_([AN])_Rep([0-9]+) ]]; then
        dim="${BASH_REMATCH[1]}"
        stiff_code="${BASH_REMATCH[2]}"
        oxy_code="${BASH_REMATCH[3]}"
        rep="${BASH_REMATCH[4]}"
        
        # Convert stiffness
        if [[ "$stiff_code" == "0_3" ]]; then
            stiff="0.3"
        elif [[ "$stiff_code" == "1_5" ]]; then
            stiff="1.5"
        else
            stiff="CHECK_ME"
        fi
        
        # Convert oxygen
        if [[ "$oxy_code" == "A" ]]; then
            oxy="anoxia"
        elif [[ "$oxy_code" == "N" ]]; then
            oxy="normoxia"
        else
            oxy="CHECK_ME"
        fi
        
        condition="${dim}_${stiff}_${oxy}_Rep${rep}"
        sample="${condition}"
        
        CONDITION_COUNT["$condition"]=$(( ${CONDITION_COUNT["$condition"]:-0} + 1 ))
        
        echo -e "${sample}\t${dim}\t${stiff}\t${oxy}\tRep${rep}\t${condition}\t${r1}\t${r2}" >> "$SAMPLE_SHEET"
    else
        echo "WARNING: Could not parse filename: $base"
    fi
done

echo "=== Condition Summary ==="
for cond in $(printf '%s\n' "${!CONDITION_COUNT[@]}" | sort); do
    echo "$cond : ${CONDITION_COUNT[$cond]} sample(s)"
done
echo

echo "Sample sheet saved to: $SAMPLE_SHEET"
echo
echo "=== DIAGNOSTIC COMPLETE ==="
echo "Paired samples       : ${#PAIRED[@]}"
echo "Conditions detected  : ${#CONDITION_COUNT[@]}"
echo
echo "Next step: Review $SAMPLE_SHEET and let me know if everything looks correct."
echo "Then we can build the full preprocessing + DESeq2 pipeline."
