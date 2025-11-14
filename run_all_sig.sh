#!/usr/bin/env bash
set -euo pipefail

### === User-configurable section ===
# Input and output paths
RDS_DIR="/Users/jiajunzhan/Desktop/RA/PDAC/Result/mutationtimer/rdata"
OUT_DIR="/Users/jiajunzhan/Desktop/RA/PDAC/Code/WGS/evolution/signature/output"
SCRIPT="/Users/jiajunzhan/Desktop/RA/PDAC/Code/WGS/evolution/signature/PCAWG11-Timing_and_Signatures-master/PCAWG_signatures_fix.R"

# COSMIC references
SBS_PATH="/Users/jiajunzhan/Desktop/RA/PDAC/Code/WGS/evolution/signature/PCAWG11-Timing_and_Signatures-master/COSMIC_v3.4_SBS_GRCh38.txt"
DBS_PATH="/Users/jiajunzhan/Desktop/RA/PDAC/Code/WGS/evolution/signature/PCAWG11-Timing_and_Signatures-master/COSMIC_v3.4_DBS_GRCh38.txt"

### =====================

# Check required files/directories exist
if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: R script not found: $SCRIPT" >&2
  exit 1
fi
if [[ ! -f "$SBS_PATH" ]]; then
  echo "ERROR: SBS reference not found: $SBS_PATH" >&2
  exit 1
fi
if [[ ! -f "$DBS_PATH" ]]; then
  echo "ERROR: DBS reference not found: $DBS_PATH" >&2
  exit 1
fi
if [[ ! -d "$RDS_DIR" ]]; then
  echo "ERROR: RDS directory not found: $RDS_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

shopt -s nullglob
# Only iterate over *.objects.rds
for RDS in "$RDS_DIR"/*.objects.rds; do
  SAMPLE="$(basename "$RDS" .objects.rds)"
  OUT_PREFIX="${OUT_DIR}/${SAMPLE}"
  OUTPUT_FILE="${OUT_PREFIX}/${SAMPLE}_sig_weights.txt"

  # If output already exists, skip
  if [[ -f "$OUTPUT_FILE" ]]; then
    echo "⏭️  Skipping ${SAMPLE} (already has output: ${OUTPUT_FILE})"
    echo "---------------------------------------------"
    continue
  fi

  echo ">>> Running sample: ${SAMPLE}"
  mkdir -p "$OUT_PREFIX"

  # Arguments must match the R script order:
  #   id, obj_rds, sbs_path, dbs_path, out_prefix, n_boot
  Rscript "$SCRIPT" \
    "$SAMPLE" \
    "$RDS" \
    "$SBS_PATH" \
    "$DBS_PATH" \
    "$OUT_PREFIX"

  echo "✅ Finished ${SAMPLE}"
  echo "---------------------------------------------"
done