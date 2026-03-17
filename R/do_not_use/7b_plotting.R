.libPaths('/u/area/vgazziero/R/rstudio_4_4')
library(tidyverse)
rm(list=ls())
setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')

# loading
testing_sublinear_all_genes_rna = readRDS('data/testing_sublinear_all_genes_rna.rds')
testing_sublinear_all_genes_prot = readRDS('data/testing_sublinear_all_genes_prot.rds')

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes_v5.rds')

transcriptomics_data = transcriptomics_data %>% 
  dplyr::mutate(is_mutated = ifelse(mut_consequence == 'synonymous_variant', FALSE, is_mutated)) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(expression = value) %>% 
  dplyr::filter(!is.na(expression))

proteomics_data = readRDS('data/proteogenomics_data_all_genes_new_norm_v4.rds')

proteomics_data = proteomics_data %>% 
  dplyr::mutate(is_mutated = ifelse(mut_consequence == 'synonymous_variant', FALSE, is_mutated)) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(expression = mean_intensity) %>% 
  dplyr::filter(!is.na(expression))


all_anova_prot = lapply(testing_sublinear_all_genes_prot %>% names, function(g) {
  testing_sublinear_all_genes_prot[[g]]$anova %>%
    as.data.frame() %>%
    mutate(gene = g) %>% 
    rownames_to_column('test')
}) %>%
  bind_rows()
#
classes_genes = all_anova_prot %>% 
  mutate(significance = ifelse(`Pr(>F)` <= 0.05, 'significant', 'ns')) %>% 
  filter(test == 'H1: h != 1') 

res_test = lapply(testing_sublinear_all_genes_prot %>% names, function(g) {
  
  testing_sublinear_all_genes_prot[[g]]$fit_alt %>% 
    broom::tidy() %>% 
    mutate(gene = g)
  
}) %>% bind_rows() %>% 
  full_join(., classes_genes, by = 'gene') %>% 
  filter(term == 'scaled_cna')

res_test = res_test %>% 
  mutate(class = case_when(
    significance == 'ns' ~ 'linear', 
    (significance != 'ns' & estimate < 1) ~ 'sublinear',
    (significance != 'ns' & estimate > 1) ~ 'superlinear'
  ))

res_test %>% 
  ggplot(aes(class, fill = class)) + 
  geom_bar(stat = 'count') + 
  theme_bw() + 
  ggsci::scale_fill_bmj() + 
  labs(y = 'Number of genes', x = 'Classes')
ggsave('res/results_protein.png', width = 6, height = 6)
  

plot_classes_genes = function(cls, colors) {
  cls %>% 
    ggplot(aes(significance, fill = significance)) + 
    geom_bar(stat = 'count') + 
    theme_bw()
}




nonlinear_genes = all_anova %>%
  filter(`Pr(>F)` <= 0.05) %>%
  pull(gene)

pik3ca = plot_expression_with_regression_lines(proteomics_data, gene = 'GTF2B', fit = testing_sublinear_all_genes_prot, K_ref = 2)
ggsave(plot = pik3ca, 'res/pik3ca_sublinear.png', width = 6, height = 6)



#
non_linear_coefficients = lapply(nonlinear_genes, function(g) {
  broom::tidy(testing_sublinear_all_genes_rna[[g]]$fit_alt) %>%
    mutate(gene = g) %>%
    bind_cols(., confint(testing_sublinear_all_genes_rna[[g]]$fit_alt)) %>%
    rename(CI_low = '2.5 %') %>%
    rename(CI_high = '97.5 %')
}) %>% bind_rows() %>%
  relocate(gene, .before = term)

non_linear_coefficients_v2 = lapply(nonlinear_genes, function(g) {
  broom::tidy(testing_sublinear_all_genes_rna[[g]]$fit_null) %>%
    mutate(gene = g) %>%
    select(gene, estimate) %>%
    rename(intercept_null = estimate)
}) %>% bind_rows() #

non_linear_coefficients = non_linear_coefficients %>%
  select(gene, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`)

df = proteomics_data %>%
  filter(hgnc_symbol == 'PARP16') %>%
  filter(mutation_status == 'Wild-type') %>%
  group_by(hgnc_symbol) %>%
  full_join(., non_linear_coefficients, by = join_by('hgnc_symbol' == 'gene')) %>%
  full_join(., non_linear_coefficients_v2, by = join_by('hgnc_symbol' == 'gene')) %>%
  mutate(cna_log = log(tot_cna / K_ref))





df %>%
  # filter(term != '(Intercept)') %>%
  ggplot(aes(x = cna_log, y = expression)) +
  geom_point() +
  geom_abline(aes(intercept = intercept_null,
                  slope = 1),
              linetype = "dashed",
              size = 1.2,
              alpha = 0.8,
              color = "red") +
  geom_abline(aes(intercept = intercept,
                  slope = scaled_cna),
              # linetype = "dashed",
              size = 1.2,
              alpha = 0.8,
              color = "blue") +
  facet_wrap(~hgnc_symbol) +
  theme_bw() +
  labs(title = 'Sublinearity test (results for non linear genes)',
       subtitle = "Red dashed: linear expectation (slope = 1)\nBlue: fitted slope",
       x = "log(CNA/2)", y = "vst expression (log scaled)")




# transcriptomics_data %>%
#   filter(hgnc_symbol == 'ATM') %>%
#   filter(mutation_status == 'Wild-type') %>%
#   mutate(scaled_cna = log(tot_cna / K_ref)) %>%
#   ggplot(aes(x = scaled_cna, y = expression)) +
#   geom_point() +
#   ylim(c(0, 14)) +
#   geom_abline(intercept = coef(testing_sublinear_all_genes_rna[['ATM']]$fit_null)[1], slope = 1,
#               linetype = "dashed", size = 1.2, alpha = 0.8,
#               color = "red") +
#   geom_abline(intercept = coef(testing_sublinear_all_genes_rna[['ATM']]$fit_alt)[1], slope = coef(testing_sublinear_all_genes_rna[['ATM']]$fit_alt)[2],
#               linetype = "dashed", size = 1.2, alpha = 0.8,
#               color = "blue")


g = 'SMAD2'
transcriptomics_data %>%
  filter(hgnc_symbol == 'SMAD2') %>%
  filter(mutation_status == 'Wild-type') %>%
  mutate(scaled_cna = log(tot_cna / K_ref)) %>%
  ggplot(aes(x = scaled_cna, y = expression)) +
  geom_point() +
  geom_abline(intercept = coef(testing_sublinear_all_genes_rna[[g]]$fit_null)[1], slope = 1,
              linetype = "dashed", size = 1.2, alpha = 0.8,
              color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Sublinearity test for", g),
       subtitle = "Red dashed: linear expectation (slope = 1)\nBlue: fitted slope",
       x = "log(CNA/K_ref)", y = "vst expression")
