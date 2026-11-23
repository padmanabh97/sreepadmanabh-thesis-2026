#!/bin/bash
set -euo pipefail

# =============================================================================
# SYMLINKS AND PATH SETUP (Cleaned Version)
# Purpose: Rename raw fastq folder and update sample_sheet.tsv
# Note: Symlink creation is now handled by master_setup.sh
# =============================================================================

echo "=== Symlinks and Path Setup ==="
echo

SYMLINK_PATH="/mnt/d/OxoMech"
SAMPLE_SHEET="$SYMLINK_PATH/sample_sheet.tsv"

# --- Step 1: Rename raw fastq folder ---
echo "Step 1: Renaming '2D raw fastq files' to 'raw_fastq' ..."

OLD_FASTQ_DIR="$SYMLINK_PATH/2D raw fastq files"
NEW_FASTQ_DIR="$SYMLINK_PATH/raw_fastq"

if [[ -d "$NEW_FASTQ_DIR" ]]; then
    echo "Folder 'raw_fastq' already exists. Skipping rename."
elif [[ -d "$OLD_FASTQ_DIR" ]]; then
    mv "$OLD_FASTQ_DIR" "$NEW_FASTQ_DIR"
    echo "Renamed: $OLD_FASTQ_DIR -> $NEW_FASTQ_DIR"
else
    echo "WARNING: Could not find folder to rename: $OLD_FASTQ_DIR"
fi

# --- Step 2: Update sample_sheet.tsv ---
echo
echo "Step 2: Updating paths in sample_sheet.tsv ..."

if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "ERROR: sample_sheet.tsv not found at $SAMPLE_SHEET"
    exit 1
fi

# Replace subfolder name (removes space)
sed -i 's|2D raw fastq files|raw_fastq|g' "$SAMPLE_SHEET"

echo "sample_sheet.tsv updated successfully."

echo
echo "=== Path Setup Complete ==="
echo "You can now run the preprocessing pipeline."
