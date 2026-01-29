library(tidyverse)

rna_fit = readRDS('data/glm_fit_transcriptomics.rds')
prot_fit = readRDS('data/glm_fit_proteomics.rds')

rna_fit = Filter(function(x) !all(is.null(x)), rna_fit) 
prot_fit = Filter(function(x) !all(is.null(x)), prot_fit) 

index_errors = which(lapply(rna_fit, function(x) {length(x)}) %>% unlist() == 1)
rna_fit_filtered = rna_fit[-index_errors]
names(rna_fit_filtered) = lapply(rna_fit_filtered, function(x) {x$gene}) %>% unlist

index_errors = which(lapply(prot_fit, function(x) {length(x)}) %>% unlist() == 1)
prot_fit_filtered = prot_fit[-index_errors]
names(prot_fit_filtered) = lapply(prot_fit_filtered, function(x) {x$gene}) %>% unlist

rna_res = lapply(rna_fit_filtered, function(x) {
  
  coeffs = x[['coefficients']] %>% 
    as.matrix()
  
  info = rownames(coeffs)  
  coeffs %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('info') %>% 
    dplyr::mutate(gene = x[['gene']])
}) %>% 
  bind_rows()

prot_res = lapply(prot_fit_filtered, function(x) {
  
  coeffs = x[['coefficients']] %>% 
    as.matrix()
  
  info = rownames(coeffs)  
  coeffs %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('info') %>% 
    dplyr::mutate(gene = x[['gene']])
}) %>% 
  bind_rows()


# visualize the results
rna_res %>% 
  dplyr::filter(info != '(Intercept)') %>% 
  # tidyr::pivot_wider(names_from = info, values_from = s0) %>% 
  # filter(tot_cna > 0) %>% 
  filter(s0 > 0) %>% 
  ggplot(aes(s0)) +
  geom_histogram(binwidth = 0.01) + 
  theme_bw() + 
  facet_wrap(vars(info), scales = 'free') +
  ggtitle('RNA')

prot_res %>% 
  dplyr::filter(info != '(Intercept)') %>% 
  tidyr::pivot_wider(names_from = info, values_from = s0) %>%
  filter(tot_cna > 0) %>%
  filter(s0 > 0) %>% 
  ggplot(aes(s0)) +
  geom_histogram(binwidth = 0.01) + 
  theme_bw() + 
  facet_wrap(vars(info), scales = 'free') +
  ggtitle('Protein')
