#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 10:00:00
#SBATCH --mail-user=jaileru@unc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

REPO="/nas/longleaf/home/jaileru/hepatic-nutrient-regulation"

cd "$REPO"
mkdir -p logs

pixi run --manifest-path hnr-workspace/pixi.toml snakemake --cores 4 --sdm conda
