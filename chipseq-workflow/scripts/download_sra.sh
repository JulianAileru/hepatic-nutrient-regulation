#!/usr/bin/env bash
#SBATCH --job-name=chipseq-sra
#SBATCH --array=0-11%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=chipseq-workflow/slurm-logs/sra-%A_%a.out
#SBATCH --error=chipseq-workflow/slurm-logs/sra-%A_%a.err

set -euo pipefail

# Submit from the repository root.
REPO="$SLURM_SUBMIT_DIR"
cd "$REPO"

# Load the same modules you verified interactively.
# Replace these with the exact module names if necessary.
# module load sra-tools
# module load pigz

FASTQ_DIR="$REPO/chipseq-workflow/data/fastq"
TMP_DIR="${SLURM_TMPDIR:-$REPO/chipseq-workflow/data/tmp}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$FASTQ_DIR" "$TMP_DIR"

# Read ChIP-seq accessions from the existing SRA metadata.
RUN=$(python - "$SLURM_ARRAY_TASK_ID" <<'PY'
import csv
import sys

with open("data/metadata/SraRunTable.csv", newline="") as f:
    runs = [
        row["Run"]
        for row in csv.DictReader(f)
        if row["Assay Type"] == "ChIP-Seq"
    ]

assert len(runs) == 12, f"Expected 12 ChIP-seq runs, found {len(runs)}"
print(runs[int(sys.argv[1])])
PY
)

echo "Downloading $RUN"

if [[ -s "$FASTQ_DIR/${RUN}.fastq.gz" ]]; then
    echo "FASTQ already exists: $RUN"
    exit 0
fi

fasterq-dump "$RUN" \
    --outdir "$FASTQ_DIR" \
    --temp "$TMP_DIR" \
    --threads "$SLURM_CPUS_PER_TASK"

pigz -p "$SLURM_CPUS_PER_TASK" "$FASTQ_DIR/${RUN}.fastq"

echo "Completed $RUN"
