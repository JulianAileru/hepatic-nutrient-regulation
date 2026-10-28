#!/bin/bash
#SBATCH -J hnr-nextflow
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 3-00:00:00
#SBATCH -o logs/nextflow_%j.out
#SBATCH --mail-user=jaileru@unc.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail

source env.sh

export NXF_OPTS="-Xms512m -Xmx4g"

nextflow run nf-core/chipseq \
    -c nextflow.config \
    --input samplesheet.csv \
    --outdir results \
    --fasta /users/j/a/jaileru/hepatic-nutrient-regulation/resources/genome.fasta \
    --gtf /users/j/a/jaileru/hepatic-nutrient-regulation/resources/genome.gtf \
    --aligner bwa \
    --macs_gsize 1870000000 \
    --narrow_peak \
    -resume \
    -with-report reports/execution_report.html \
    -with-timeline reports/timeline.html \
    -with-trace reports/trace.txt

echo "ChIP-seq pipeline completed successfully"

