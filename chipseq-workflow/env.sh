#!/usr/bin/env bash

# Load HPC dependencies
module load apptainer
module load nextflow/26.10.2

# Nextflow compatibility
export NXF_SYNTAX_PARSER=v1

# Project paths
export CHIPSEQ_ROOT="/users/j/a/jaileru/hepatic-nutrient-regulation/chipseq-workflow"
export CHIPSEQ_RESOURCES="/users/j/a/jaileru/hepatic-nutrient-regulation/resources"

export CHIPSEQ_FASTA="$CHIPSEQ_RESOURCES/genome.fasta"
export CHIPSEQ_GTF="$CHIPSEQ_RESOURCES/genome.gtf"

# Nextflow cache and working directories
export NXF_WORK="$CHIPSEQ_ROOT/work"
export NXF_APPTAINER_CACHEDIR="$NXF_WORK/singularity"

# Avoid Quay.io authentication issues
unset APPTAINER_DOCKER_USERNAME
unset APPTAINER_DOCKER_PASSWORD

mkdir -p "$NXF_WORK" "$NXF_APPTAINER_CACHEDIR"
