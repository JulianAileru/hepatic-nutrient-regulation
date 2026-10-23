#!/bin/bash
#SBATCH --job-name=chipseq
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail

REPO="/nas/longleaf/home/jaileru/hepatic-nutrient-regulation"
cd "$REPO"

export NXF_SYNTAX_PARSER=v1

nextflow run nf-core/chipseq \
    -r 2.1.0 \
    -profile test \
    -work-dir "$REPO/chipseq-workflow/work" \
    --outdir "$REPO/chipseq-workflow/results"
