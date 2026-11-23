#!/bin/bash
set -euo pipefail

# =============================================================================
# MASTER SETUP SCRIPT
# 2D vs 3D Oxo-Mech RNA-seq Pipeline
#
# Responsibilities:
# - Create symlink /mnt/d/OxoMech (if it doesn't exist)
# - Make all shell scripts executable
# - Guide the user through the pipeline steps
# =============================================================================

ORIGINAL_PROJECT_DIR="/mnt/d/ESoftBio Lab/Omics Analysis/Oxo-Mech 2D vs 3D transcriptome"
SYMLINK_PATH="/mnt/d/OxoMech"

echo "=============================================="
echo "   MASTER SETUP - 2D vs 3D Oxo-Mech Pipeline"
echo "=============================================="
echo

# --- Step 1: Create symlink if it doesn't exist ---
echo "Step 1: Checking/Creating symlink..."

if [[ -L "$SYMLINK_PATH" ]]; then
    echo "  ✓ Symlink already exists: $SYMLINK_PATH"
else
    echo "  Creating symlink..."
    ln -s "$ORIGINAL_PROJECT_DIR" "$SYMLINK_PATH"
    echo "  ✓ Symlink created: $SYMLINK_PATH"
fi

echo

# --- Step 2: Make all scripts executable ---
SCRIPTS=(
    "$SYMLINK_PATH/check_mouse_fastq_diagnostic_v2.sh"
    "$SYMLINK_PATH/symlinks_and_path_setup.sh"
    "$SYMLINK_PATH/run_rnaseq_pipeline.sh"
)

echo "Step 2: Making scripts executable..."
for script in "${SCRIPTS[@]}"; do
    if [[ -f "$script" ]]; then
        chmod +x "$script"
        echo "  ✓ $(basename "$script")"
    else
        echo "  ⚠ File not found: $(basename "$script")"
    fi
done

echo
echo "All scripts are now executable."
echo

# --- Step 3: Ask if user wants to continue ---
read -p "Do you want to continue with the pipeline steps? (y/n): " RUN_ALL

if [[ "$RUN_ALL" != "y" && "$RUN_ALL" != "Y" ]]; then
    echo
    echo "Setup complete."
    echo "You can run the scripts manually in this order:"
    echo "  1. $SYMLINK_PATH/check_mouse_fastq_diagnostic_v2.sh"
    echo "  2. $SYMLINK_PATH/symlinks_and_path_setup.sh   (renames folder + updates sample sheet)"
    echo "  3. $SYMLINK_PATH/run_rnaseq_pipeline.sh"
    echo "  4. Rscript $SYMLINK_PATH/2D_vs_3D_transcriptome_DESeq2.R"
    exit 0
fi

echo
echo "Starting pipeline steps..."
echo

# --- Run steps with confirmation ---
read -p "Run Step 1: Diagnostic Script? (y/n): " ans
if [[ "$ans" == "y" || "$ans" == "Y" ]]; then
    echo ">>> Running diagnostic script..."
    "$SYMLINK_PATH/check_mouse_fastq_diagnostic_v2.sh"
fi

read -p "Run Step 2: Rename folder + Update sample sheet? (y/n): " ans
if [[ "$ans" == "y" || "$ans" == "Y" ]]; then
    echo ">>> Running path setup script..."
    "$SYMLINK_PATH/symlinks_and_path_setup.sh"
fi

read -p "Run Step 3: Preprocessing Pipeline? (y/n): " ans
if [[ "$ans" == "y" || "$ans" == "Y" ]]; then
    echo ">>> Running preprocessing pipeline..."
    "$SYMLINK_PATH/run_rnaseq_pipeline.sh"
fi

echo
echo "=============================================="
echo "Pipeline steps completed."
echo
echo "Next step: Run DESeq2 analysis using:"
echo "Rscript $SYMLINK_PATH/2D_vs_3D_transcriptome_DESeq2.R"
echo "=============================================="
