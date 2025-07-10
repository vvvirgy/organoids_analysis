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
# RNA
common_rna = intersect(names(rna_ploidy_mut_fit), names(rna_ploidy_fit))
rna_ploidy_mut_fit = rna_ploidy_mut_fit[common_rna]
rna_ploidy_fit = rna_ploidy_fit[common_rna]

# prot
common_prot = intersect(names(prot_ploidy_mut_fit), names(prot_ploidy_fit))
prot_ploidy_mut_fit = prot_ploidy_mut_fit[common_prot]
prot_ploidy_fit = prot_ploidy_fit[common_prot]

# minimize the error!
rna = lapply(common_rna, function(x) {
  md1 = rna_ploidy_fit[[x]]$cv_glm_best_alpha
  md2 = rna_ploidy_mut_fit[[x]]$cv_glm_best_alpha
  
  compare_models(md1, md2)
})
names(rna) = common_rna

rna = tibble(
  gene = names(rna),
  model_diff_err = unlist(unname(rna))
)

rna_improvement = lapply(common_rna, function(x) {
  md1 = rna_ploidy_fit[[x]]$cv_glm_best_alpha
  md2 = rna_ploidy_mut_fit[[x]]$cv_glm_best_alpha
  
  compute_improvement(md1, md2)
})
names(rna_improvement) = common_rna
rna_improvement = tibble(
  gene = names(rna_improvement),
  improvement = unlist(unname(rna_improvement))
)

rr = rna_improvement %>% 
  ggplot(aes(improvement)) + 
  geom_histogram(binwidth = 0.1) + 
  theme_bw() + 
  ggtitle('model performance difference estimation - RNA')

rna %>% 
  ggplot(aes(model_diff_err)) + 
  geom_histogram(binwidth = 0.01) + 
  theme_bw() + 
  ggtitle('model difference error RNA')

prot = lapply(common_prot, function(x) {
  md1 = prot_ploidy_fit[[x]]$cv_glm_best_alpha
  md2 = prot_ploidy_mut_fit[[x]]$cv_glm_best_alpha
  
  compare_models(md1, md2)
})
names(prot) = common_prot

prot = tibble(
  gene = names(prot),
  model_diff_err = unlist(unname(prot))
)
prot_improvement = lapply(common_prot, function(x) {
  md1 = prot_ploidy_fit[[x]]$cv_glm_best_alpha
  md2 = prot_ploidy_mut_fit[[x]]$cv_glm_best_alpha
  
  compute_improvement(md1, md2)
})
names(prot_improvement) = common_prot
prot_improvement = tibble(
  gene = names(prot_improvement),
  improvement = unlist(unname(prot_improvement))
)

pp = prot_improvement %>% 
  ggplot(aes(improvement)) + 
  geom_histogram(binwidth = 0.1) + 
  theme_bw() + 
  ggtitle('model performance difference estimation - Protein')

prot %>% 
  ggplot(aes(model_diff_err)) + 
  geom_histogram(binwidth = 0.01) + 
  theme_bw() + 
  ggtitle('model difference error Proteins')

all_fits = list(
  'prot_ploidy_fit' = prot_ploidy_fit, 
  'prot_ploidy_mut_fit' = prot_ploidy_mut_fit, 
  'rna_ploidy_fit' = rna_ploidy_fit, 
  'rna_ploidy_mut_fit' = rna_ploidy_mut_fit
)
coeffs = lapply(all_fits, function(x) {
  lapply(x, function(g) {
    coef(g$cv_glm_best_alpha, s = 'lambda.min')
  })
})

saveRDS(coeffs, 'data/coefficients_fits_glm.rds')

og_coad = transcriptomics_data %>% 
  dplyr::filter(CGC_role_COAD == 'oncogene') %>% 
  pull(hgnc_symbol)

coeffs_og_coad = lapply(coeffs, function(x) {
  x[og_coad]
})

