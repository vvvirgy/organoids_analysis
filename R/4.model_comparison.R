# comparing fits results
library(glmnet)
library(tidyverse)

# transcriptomics
rna_ploidy_mut_fit = readRDS('data/glm_fit_ridge_transcriptomics.rds')
rna_ploidy_mut_fit = Filter(Negate(is.null), rna_ploidy_mut_fit)
names(rna_ploidy_mut_fit) = lapply(rna_ploidy_mut_fit, function(x) {x$gene}) %>% unlist

rna_ploidy_fit = readRDS('data/glm_fit_ridge_transcriptomics_cna_only.rds')
rna_ploidy_fit = Filter(Negate(is.null), rna_ploidy_fit)
names(rna_ploidy_fit) = lapply(rna_ploidy_fit, function(x) {x$gene}) %>% unlist

# proteomics
prot_ploidy_fit = readRDS('data/glm_fit_ridge_proteomics_cna_only.rds')
prot_ploidy_fit = Filter(Negate(is.null), prot_ploidy_fit)
names(prot_ploidy_fit) = lapply(prot_ploidy_fit, function(x) {x$gene}) %>% unlist

prot_ploidy_mut_fit = readRDS('data/glm_fit_ridge_proteomics.rds')
prot_ploidy_mut_fit = Filter(Negate(is.null), prot_ploidy_mut_fit)
names(prot_ploidy_mut_fit) = lapply(prot_ploidy_mut_fit, function(x) {x$gene}) %>% unlist

# select only the genes that are avaliable for both models



delta_rna  <- cv_rna_ploidy$cvm[cv_rna_ploidy$lambda == cv_rna_ploidy$lambda.min] - 
  cv_rna_full$cvm[cv_rna_full$lambda == cv_rna_full$lambda.min]

delta_prot <- cv_prot_ploidy$cvm[cv_prot_ploidy$lambda == cv_prot_ploidy$lambda.min] - 
  cv_prot_full$cvm[cv_prot_full$lambda == cv_prot_full$lambda.min]
