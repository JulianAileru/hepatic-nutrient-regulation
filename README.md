# Hepatic nutrient regulation

Reanalysis of mouse liver RNA-seq and ChIP-seq from male C57BL/6J mice kept on
three diets: normal chow, high-fat diet (HFD) and calorie restriction (CR). The
RNA-seq measures how each diet changes liver gene expression. The ChIP-seq maps
where the nuclear receptors PPARα and RXRα bind in CR and HFD liver.

| | |
|---|---|
| Source | NCBI BioProject [PRJNA208469](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA208469), SRA study SRP026039 |
| Organism | *Mus musculus*, C57BL/6J, male |
| Sequencer | Illumina Genome Analyzer II |
| Reference | GRCm39, Ensembl release 114 (shared by both pipelines) |
| Sample metadata | [data/metadata/SraRunTable.csv](data/metadata/SraRunTable.csv) |

The BioProject also contains 4 DNase-seq runs (SRR901061–SRR901064). They are not
analyzed here.

## Status

| Analysis | Status |
|---|---|
| RNA-seq | Complete. Results are in `results/`. |
| ChIP-seq | Incomplete. The 24 Oct 2026 run stopped at peak calling because quay.io was temporarily unavailable when the MACS3 container was downloaded. FastQC and trimming finished for all 12 libraries, and filtered BAMs exist for 7 of 12. The run needs to be resumed (see [ChIP-seq](#chip-seq)). |

## Repository layout

```
config/                 RNA-seq configuration (config.yaml, samples.tsv, units.tsv)
workflow/Snakefile      Imports snakemake-workflows/rna-seq-star-deseq2 v3.1.2
run_workflow.sh         SLURM submission script for the RNA-seq workflow
hnr-workspace/          pixi environment (Snakemake + SLURM executor plugin)
data/metadata/          SRA run table, accession list, RNA-seq sample sheet
chipseq-workflow/       nf-core/chipseq setup: run script, config, samplesheet, helper scripts
  scripts/              download_sra.sh, generate_samplesheet.py
  reports/              Nextflow execution report and timeline
```

These are generated, are not tracked in git, and are recreated by the pipelines:
`resources/` (genome, annotation, STAR index), `results/`,
`chipseq-workflow/results/`, `chipseq-workflow/work/`, FASTQs and BAMs.

---

## RNA-seq

### Samples

Paired-end, 2 × 36 bp.

| Sample | SRA run | Diet (`diet` column) |
|---|---|---|
| GSM1163164 | SRR901041 | chow |
| GSM1163165 | SRR901042 | chow |
| GSM1163166 | SRR901043 | chow |
| GSM1163167 | SRR901044 | hfd |
| GSM1163168 | SRR901045 | hfd |
| GSM1163169 | SRR901046 | calorie_restriction |
| GSM1163170 | SRR901047 | calorie_restriction |
| GSM1163171 | SRR901048 | calorie_restriction |

The workflow reads the diet labels from [config/samples.tsv](config/samples.tsv).
The `high_fat` label in `data/metadata/samples.tsv` is not used.

### Method

The workflow is [snakemake-workflows/rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
v3.1.2, configured in [config/config.yaml](config/config.yaml).

1. **Download and reference.** Reads come directly from SRA. The genome and
   annotation are Ensembl 114 GRCm39.
2. **Trimming.** fastp.
3. **Alignment and counting.** STAR 2.7.11b aligns the reads, and STAR's
   `--quantMode GeneCounts` counts reads per gene. The library is unstranded
   (RSeQC `infer_experiment` shows about a 40/40 split), so the unstranded count
   column is used.
4. **DESeq2.**
   - Design `~ diet`, with `chow` as the reference level.
   - Genes with ≤1 read in total across all samples are removed (28,478 of
     78,317 genes remain).
   - Default `DESeq()` fit, Wald test.
   - Log2 fold changes are shrunk with `lfcShrink(type = "ashr")`.
5. **PCA.** `rlog`-transformed counts (`blind = FALSE`), coloured by diet.

**Contrasts.** A positive log2 fold change means the gene is higher in the first group:

| Contrast | Comparison |
|---|---|
| `hfd_vs_chow` | HFD vs chow |
| `cr_vs_chow` | CR vs chow |
| `cr_vs_hfd` | CR vs HFD |

### Running

```bash
sbatch run_workflow.sh   # from the repo root
```

The script submits Snakemake through the SLURM executor (up to 20 concurrent jobs,
conda software environments). If a job fails, it is retried once with more memory
and runtime.

### Results guide

All paths are relative to `results/`.

#### `counts/`: raw counts

| File | Contents |
|---|---|
| `all.tsv` | Raw read counts. Rows are Ensembl gene IDs (`ENSMUSG…`), columns are samples. This is the only input to DESeq2. |
| `all.symbol.tsv` | The same matrix with gene symbols in place of IDs. Display only, not used downstream. |

Raw counts are not normalized for sequencing depth, so compare samples with
`deseq2/normcounts.tsv` instead.

#### `deseq2/`: fitted model and normalized counts

| File | Contents |
|---|---|
| `all.rds` | The fitted `DESeqDataSet`. Load it in R with `readRDS()` to run new contrasts, `vst()`/`rlog()`, or `plotCounts()`. |
| `normcounts.tsv` | Counts divided by DESeq2 size factors (sequencing-depth normalized), after the ≤1-read filter. |
| `normcounts.symbol.tsv` | The same, with gene symbols. |

Normalized counts are comparable between samples for the same gene. They are not
adjusted for gene length, so they should not be used to compare expression
between different genes.

#### `diffexp/`: differential expression, one set per contrast

| File | Contents |
|---|---|
| `<contrast>.diffexp.tsv` | DESeq2 results with Ensembl IDs, sorted by `padj`. |
| `<contrast>.diffexp.symbol.tsv` | The same, with gene symbols. |
| `<contrast>.ma-plot.svg` | MA plot (y-axis fixed to ±2). |

Columns in the results tables:

| Column | Meaning |
|---|---|
| `gene` | Ensembl ID or gene symbol |
| `baseMean` | Mean normalized count across all 8 samples |
| `log2FoldChange` | ashr-shrunken log2 fold change (first group / second group) |
| `lfcSE` | Posterior standard deviation of the shrunken log2 fold change |
| `pvalue` | Wald test p-value. It comes from the unshrunken estimate, so it isn't affected by shrinkage. |
| `padj` | Benjamini–Hochberg adjusted p-value |

Notes for reading the tables:

- **`padj` is `NA` for 10,000–14,000 genes per contrast.** DESeq2's independent
  filtering removes these low-count genes before multiple-testing correction.
  Treat them as "not tested", not "not significant".
- **The fold changes are shrunk.** Shrinkage pulls the fold changes of noisy,
  low-count genes toward 0. That makes them suitable for ranking and for GSEA,
  but they are smaller than unshrunken values.
- **Some symbols are duplicated.** 18 symbols map to more than one Ensembl ID.
  Genes with no symbol keep their Ensembl ID. Use the `.diffexp.tsv` files when
  you need unique identifiers.
- **Significant genes per contrast** (`padj < 0.05`, \|log2FC\| ≥ 1):

  | Contrast | DE genes | Up | Down |
  |---|---|---|---|
  | `hfd_vs_chow` | 1,617 | 874 | 743 |
  | `cr_vs_chow` | 503 | 314 | 189 |
  | `cr_vs_hfd` | 3,367 | 1,663 | 1,704 |

#### `qc/` and `pca.diet.svg`

| Path | Contents |
|---|---|
| `qc/multiqc_report.html` | Combined QC report: fastp, STAR, RSeQC |
| `qc/rseqc/` | Per-sample RSeQC output: strandedness, read distribution, duplication, GC, junctions, inner distance |
| `pca.diet.svg` | PCA of rlog-transformed counts, coloured by diet |
| `star/<sample>/` | STAR BAMs, `Log.final.out` (mapping stats), `ReadsPerGene.out.tab` |

### QC summary

| Sample | Diet | Read pairs | Uniquely mapped |
|---|---|---|---|
| GSM1163164 | chow | 36.6 M | 73.7% |
| GSM1163165 | chow | 31.1 M | 73.2% |
| GSM1163166 | chow | 34.2 M | 72.9% |
| GSM1163167 | hfd | 36.3 M | 73.5% |
| GSM1163168 | hfd | 24.6 M | 78.8% |
| GSM1163169 | CR | 37.4 M | 79.4% |
| GSM1163170 | CR | 19.5 M | 76.6% |
| GSM1163171 | CR | 30.1 M | 79.5% |

Points to keep in mind when interpreting the results:

- **All samples are male.** Y-chromosome genes are expressed and *Xist* is about 0,
  which matches the SRA metadata.
- **HFD has only 2 replicates, and they differ.** GSM1163167 is the most distinct
  sample in the dataset (Spearman ρ 0.78–0.87 with every other sample).
  GSM1163168 falls between chow and GSM1163167 on the PCA. Both show the same
  direction of change, but effect sizes in the HFD contrasts depend heavily on
  GSM1163167.
- **The ribosomal-protein share of reads varies tenfold** (0.3–3.2% of normalized
  counts): lowest in HFD and highest in CR. This is a likely source of technical
  difference between samples. Many *Rpl*/*Rps* genes in the `cr_vs_hfd` results
  probably reflect this rather than biology.
- **STAR's per-base mismatch rate is high** (1.9–3.2%). That is typical of older
  Genome Analyzer II data.

---

## ChIP-seq

### Samples

12 libraries, single-end 35–36 bp, two biological replicates per group. The
samplesheet ([chipseq-workflow/samplesheet.csv](chipseq-workflow/samplesheet.csv))
was built from the SRA run table by
[generate_samplesheet.py](chipseq-workflow/scripts/generate_samplesheet.py).
Each ChIP replicate uses the IgG library with the same diet and replicate number
as its control.

| Diet | PPARα | RXRα | IgG control |
|---|---|---|---|
| CR | SRR901050, SRR901054 | SRR901049, SRR901053 | SRR901057, SRR901058 |
| HFD | SRR901052, SRR901056 | SRR901051, SRR901055 | SRR901059, SRR901060 |

The study has no ChIP-seq for chow-fed mice.

### Method

The pipeline is [nf-core/chipseq](https://nf-co.re/chipseq) v2.1.0, run with
Nextflow 26.10.2 on SLURM with Apptainer containers.

- **Reads:** FastQC, then Trim Galore 0.6.7.
- **Alignment:** BWA-MEM 0.7.18 to the same GRCm39 / Ensembl 114 genome as the
  RNA-seq.
- **Library merging and filtering:** Picard MarkDuplicates, and filtering out
  duplicates, multi-mappers and blacklisted regions.
- **Peak calling:** MACS3 against the matched IgG control, with
  `--macs_gsize 1870000000`.
- **Downstream:** HOMER peak annotation, consensus peaks with featureCounts and
  DESeq2 across conditions, deepTools bigWigs and plots, and MultiQC.

### Running

```bash
# 1. Download FASTQs (array job; run from the repo root, needs sra-tools and pigz)
sbatch chipseq-workflow/scripts/download_sra.sh

# 2. Run the pipeline (from chipseq-workflow/; -resume reuses completed tasks)
cd chipseq-workflow
sbatch run_nf_chipseq.sh
```

[env.sh](chipseq-workflow/env.sh) sets the paths and loads Apptainer.
[nextflow.config](chipseq-workflow/nextflow.config) sets the SLURM executor and a
1-hour container download timeout. Containers are cached in
`chipseq-workflow/work/singularity/`.

### Outputs (`chipseq-workflow/results/`)

| Path | Contents | Status |
|---|---|---|
| `fastqc/`, `trimgalore/` | Read QC and trimming | Done |
| `bwa/merged_library/*.mLb.clN.sorted.bam` | Filtered, deduplicated alignments | 7 of 12 |
| `bwa/merged_library/{picard_metrics,preseq,samtools_stats}/` | Alignment QC | Partial |
| `bwa/merged_library/macs3/` | Peaks, consensus peaks, HOMER annotation, peak-level DESeq2 | Pending |
| `bwa/merged_library/bigwig/`, `deeptools/` | Coverage tracks and profile plots | Pending |
| `multiqc/` | Combined QC report | Pending |

---

## Software

| Component | Version |
|---|---|
| Snakemake | ≥ 9.27.1 (pixi, `hnr-workspace/`) |
| snakemake-executor-plugin-slurm | ≥ 2.8.0 |
| rna-seq-star-deseq2 | v3.1.2 |
| STAR | 2.7.11b |
| nf-core/chipseq | v2.1.0 (revision b6c771e) |
| Nextflow | 26.10.2 |

The exact versions of every tool are recorded in `hnr-workspace/pixi.lock`, the
Snakemake conda environments, and
`chipseq-workflow/results/pipeline_info/nf_core_chipseq_software_mqc_versions.yml`.

## Biological interpretation

To be added.
