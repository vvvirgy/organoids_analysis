# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

R-based bioinformatics analysis of single-cell RNA-seq (scRNA), whole-genome sequencing (WGS), and bulk proteomics data from colorectal cancer organoids. The central question is how copy number alterations (CNAs) relate to transcriptomic and proteomic expression changes — specifically, whether cells "compensate" for gene dosage effects at the RNA or protein level (the **Compensation Score**).

## Running analyses

Scripts are run directly with Rscript. There is no build system or test suite. Heavy computations run on the ORFEO HPC cluster via SLURM:

```bash
# Submit SLURM jobs (run from the cluster)
sbatch jovoniR/bash_scripts/run_rna.sh          # scRNA differential expression (200 GB, 4 h)
sbatch jovoniR/bash_scripts/run_dna.sh          # DNA differential expression (64 GB, 2 h)
sbatch jovoniR/bash_scripts/run_multimodal_analysis.sh  # array job over method x stability combos

# Run a single analysis script locally
Rscript R/3.compensation_score/1.computing_cs.R
```

SLURM scripts activate the `process` conda environment before calling Rscript.

## Data paths

Hardcoded paths in scripts point to the ORFEO cluster. The canonical roots are:

- Data: `/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/`
- Results: `/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/res/` and `/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_proj/results/`

`R/functions_utils/constants.R` defines `data_path` and `res_path` variables used throughout.

## Analysis pipeline

The numbered directories under `R/` encode execution order:

| Step | Directory | What it does |
|------|-----------|-------------|
| 1 | `R/1.preprocessing/` | Sample name mapping; CNAqc peak analysis on WGS; proteomics cleaning; scRNA QC filtering and pseudobulk generation |
| 2 | `R/2.differential_expression/` | Protein differential expression using DEP + limma (vsn normalization + `lmFit`/`eBayes`) |
| 3 | `R/3.compensation_score/` | Compute CS = `|DNA_lfc − expression_lfc|`; classify genes; enrichment analysis |
| 4 | `R/4.functional_enrichment/` | GSEA with Hallmark/Reactome gene sets; semantic similarity clustering |
| 5 | `R/5.noise_model/` | Noise-aware differential testing (diploid vs. tetraploid) using limma on a bootstrapped design matrix |

`R/functions_utils/` holds shared utilities sourced at the top of most scripts. `R/do_not_use/` and `old_scripts/` contain deprecated code — do not modify or rely on these.

## Figures

- `jovoniR/` contains sub-analysis scripts that produce intermediate `.rds` plot objects saved under `jovoniR/img/`
- `prep_figures.R` and `main1.R` at the repo root assemble publication figures from those `.rds` objects using `patchwork`, then write to `img/`
- The `img/figs/` subdirectory holds additional assembled figures

## Key packages

`tidyverse`, `DEP`, `limma`, `CNAqc`, `patchwork`, `SummarizedExperiment`, `DESeq2`, `clusterProfiler`, `ReactomePA`

## Compensation Score definition

CS for a given gene × karyotype is:

```
CS = if (DNA_lfc > 0):  DNA_lfc − expression_lfc
     else:              expression_lfc − DNA_lfc
```

A positive CS means expression is attenuated relative to DNA dosage (compensation). Scripts in `R/3.compensation_score/` operate on the merged RNA/protein LFC table joined with `DNA_lfc.rds`.
