#!/bin/bash
DATA_DIR="$WD_SEQ/data"

SELECTION="$WD_SEQ/inputs/selected_probes.tsv"
PROBE_DIR="$DATA_DIR/probes/tsv"
QC_DIR="$DATA_DIR/qc"
REF_DIR="$DATA_DIR/reference"
HUMAN_REF_DIR="$WD_SEQ/references"

SCRIPT_QC="$(dirname "$0")/4_select_and_qc.py"
SCRIPT_BLAST="$(dirname "$0")/5_blast_human.py"
SCRIPT_REF="$(dirname "$0")/6_generate_reference.py"

echo "[Step 1] Running probe selection and QC..."
python "$SCRIPT_QC" \
    --selection "$SELECTION" \
    --probe_dir "$PROBE_DIR" \
    --output_dir "$QC_DIR"

echo "[Step 2] BLASTing selected probes against human RNA..."
python "$SCRIPT_BLAST" \
    --probes "$QC_DIR/probe_qc.tsv" \
    --output_dir "$QC_DIR" \
    --human_blastdb "$HUMAN_REF_DIR/human_refseq/human_refseq" \
    --human_blastdb "$HUMAN_REF_DIR/human_gencode/human_gencode"

echo "[Step 3] Generating reference files..."
python "$SCRIPT_REF" \
    --input "$QC_DIR/probe_qc.tsv" \
    --output_dir "$REF_DIR" \
    --output_prefix "bacteria_probes"

echo "Done"