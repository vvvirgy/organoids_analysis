rm(list=ls())

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')

library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)


# pkgs <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment", "sesameData", "sesame")
# sapply(pkgs, require, character.only = TRUE)

# Define project
# list_of_datasets <- c("TCGA-COAD", #colon
#                       "TCGA-READ" #rectal )

project <- "TCGA-COAD"

## Build queries ##

# cna
query_cnv <- GDCquery(
  project = project,
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  # sample.type = "Primary Tumor", 
  access = 'open', 
  sample.type = c('Blood Derived Normal', 'Primary Tumor'), 
)

# query_allele_specific <- GDCquery(
#   project = project,
#   data.category = "Copy Number Variation",
#   data.type = "Allele-specific Copy Number Segment",
#   sample.type = "Primary Tumor",
#   access = 'open'
# )

# rna
query_rna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = 'STAR - Counts', 
  sample.type = "Primary Tumor",
  access = 'open'
)

# proteome
query_proteins <- GDCquery(
  project = project,
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification", 
  sample.type = "Primary Tumor",
  access = 'open'
)

# mutations 
query_muts <- GDCquery(
  project = project,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation", 
  # sample.type = "Primary Tumor",
  access = 'open'
)

# clinical
query_clinical <- GDCquery(
  project = project, 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab", 
  access = 'open'
)

queries_list = list(
  'clinical' = query_clinical, 
  muts = query_muts, 
  cna = query_cnv, 
  rna = query_rna, 
  protein = query_proteins
)

# removing duplicated samples
# identify duplicates
queries_list_v2 = c(
  lapply(queries_list[-1], function(x) {
    x$results[[1]] = x$results[[1]][which(!duplicated(x$results[[1]]$cases)), ] 
    return(x)
  }), 
  list(clinical = queries_list$clinical))
names(queries_list_v2)

# saveRDS(queries_list, 'data/TCGA/queries.rds')
# queries_list = readRDS('data/TCGA/queries.rds')

download_data = lapply(queries_list_v2, GDCdownload)

# all_data = lapply(queries_list_v2[-5], GDCprepare)

cna = GDCprepare(queries_list_v2$cna, summarizedExperiment = F)
protein = GDCprepare(queries_list_v2$protein)
rna = GDCprepare(queries_list_v2$rna)
muts = GDCprepare(queries_list_v2$muts)
clinical =  GDCprepare(queries_list_v2$clinical)

saveRDS(cna, 'data/TCGA/cna_data.rds')
saveRDS(protein, 'data/TCGA/protein_data.rds')
saveRDS(rna, 'data/TCGA/rna_data.rds')
saveRDS(muts, 'data/TCGA/muts_data.rds')
saveRDS(clinical, 'data/TCGA/clinical_data.rds')

# integrate data?
query_cnv$results[[1]] = query_cnv$results[[1]][which(!duplicated(query_cnv$results[[1]]$cases)), ] 
cna = GDCprepare(query_cnv, summarizedExperiment = F)

dw = GDCdownload(query_allele_specific)
query_allele_specific$results[[1]] = query_allele_specific$results[[1]][which(!duplicated(query_allele_specific$results[[1]]$cases)), ] 
allele_specific = GDCprepare(query_allele_specific, summarizedExperiment = F)

