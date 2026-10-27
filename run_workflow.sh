#!/bin/bash

#SBATCH -J hnr-snakemake
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 3-00:00:00
#SBATCH -o logs/controller_%j.out
#SBATCH --mail-user=jaileru@unc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

REPO="/users/j/a/jaileru/hepatic-nutrient-regulation"
cd "$REPO"

pixi run --manifest-path hnr-workspace/pixi.toml snakemake \
    --executor slurm \
    --jobs 20 \
    --sdm conda \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --latency-wait 60 \
    --retries 1 \
    --default-resources "mem_mb=8000 if attempt == 1 else 64000" "runtime=120 if attempt == 1 else 480" \
    --set-threads star_index=8 star_align=8 \
    --set-resources star_index:mem_mb=64000 star_index:runtime=360 \
                    star_align:mem_mb=40000 star_align:runtime=240
