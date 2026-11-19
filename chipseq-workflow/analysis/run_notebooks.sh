#!/bin/bash
#SBATCH -J hnr-chipseq-nb
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o logs/notebooks_%j.out

# Executes the downstream ChIP-seq notebooks in place (outputs saved into the .ipynb).
# Run from chipseq-workflow/analysis/:  sbatch run_notebooks.sh [01|02]   (default: both)
set -eo pipefail
mkdir -p logs
PIXI="pixi run --manifest-path ../../hnr-workspace/pixi.toml -e chipseq"
nbs=("01_diffbind.ipynb" "02_peak_annotation.ipynb")
[[ "$1" == "01" ]] && nbs=("01_diffbind.ipynb")
[[ "$1" == "02" ]] && nbs=("02_peak_annotation.ipynb")
for nb in "${nbs[@]}"; do
    echo "== $nb  $(date)"
    $PIXI jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 "$nb"
done
