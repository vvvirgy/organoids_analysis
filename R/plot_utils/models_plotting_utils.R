# plots utils

plot_estimates_pvals = function(x, 
                                pval_th = 0.05, 
                                what, 
                                cols = c('#86B0BD', 'grey30')
) {
  alpha_vals = setNames(nm = c(paste0('Significant (P <= ', pval_th, ')'),
                               paste0('Not significant (P > ', pval_th, ')')), 
                        object = c(1, 0.3))
  
  colors_vals = setNames(nm = c(paste0('Significant (P <= ', pval_th, ')'),
                                paste0('Not significant (P > ', pval_th, ')')), 
                         object = cols)
  
  x %>% 
    dplyr::filter(term != '(Intercept)') %>% 
    dplyr::filter(!is.na(p.value)) %>% 
    dplyr::mutate(sign = ifelse(p.value <= pval_th, 
                                paste0('Significant (P <= ', pval_th, ')'),
                                paste0('Not significant (P > ', pval_th, ')'))) %>% 
    ggplot(aes(y = -log10(p.value), x = estimate, alpha = sign, colour = sign)) + 
    geom_point() + 
    theme_bw() + 
    scale_alpha_manual(values = alpha_vals) + 
    scale_color_manual(values = colors_vals) + 
    geom_hline(yintercept = -log10(pval_th), 
               linetype = 'dashed', 
               colour = 'Firebrick', 
               show.legend = T, 
               linewidth = 1
    ) + 
    facet_wrap(~term) + 
    xlab(paste0(what, ' coefficients'))
                    
}

plot_predictor_higher_impact = function(x, palette = wes_palette("GrandBudapest2")) {
  x %>% 
    dplyr::filter(term != '(Intercept)') %>% 
    dplyr::select(term, gene, estimate, p.value, method) %>% 
    tidyr::pivot_wider(names_from = term, values_from = c(estimate, p.value)) %>% 
    mutate(higher_predictor = case_when(
      abs(estimate_ploidy_diff) > abs(estimate_mutation_multiplicity) ~ 'ploidy_diff', 
      abs(estimate_ploidy_diff) < abs(estimate_mutation_multiplicity) ~ 'mutation_multiplicity', 
      (is.na(estimate_ploidy_diff) & !is.na(estimate_mutation_multiplicity)) ~ 'mutation_multiplicity', 
      (!is.na(estimate_ploidy_diff) & is.na(estimate_mutation_multiplicity)) ~ 'ploidy_diff', 
      .default = NA
    )) %>% 
    group_by(method) %>% 
    mutate(tot_method = length(unique(gene))) %>% 
    group_by(higher_predictor, method) %>%
    reframe(frac = n()/tot_method) %>% 
    distinct() %>% 
    ggplot(aes(y = frac, x = higher_predictor, fill = method)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    theme_bw() +
    ylab('Relative fraction') + 
    xlab('Predictor with the higher impact') +
    labs(fill = 'Omics') + 
    scale_fill_manual(values = palette)
}

betas_association = function(x, color, pth = 0.1) {
  x %>% 
    dplyr::filter(term != '(Intercept)') %>% 
    dplyr::select(term, gene, estimate, p.value, method) %>% 
    # dplyr::filter(p.value <= pth) %>%
    tidyr::pivot_wider(names_from = term, values_from = c(estimate, p.value), values_fill = 0) %>% 
    filter(p.value_ploidy_diff <= pth | p.value_mutation_multiplicity <= pth) %>% 
    ggplot(aes(estimate_mutation_multiplicity, estimate_ploidy_diff)) + 
    geom_point(color = color, alpha = 0.7) + 
    theme_bw() + 
    xlab('Mutation multiplicity') +
    ylab('Ploidy diff')
}


# gettting datas 
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/classify_muts.R')
proteogenomics_data = readRDS('data/proteogenomics_data_all_genes_new_norm.rds') %>%
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>%
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>%
  dplyr::rename(protein_expression = mean_intensity) %>%
  dplyr::filter(!is.na(protein_expression))

proteogenomics_data = classify_mutations(proteogenomics_data)
proteogenomics_data = proteogenomics_data %>% 
  dplyr::filter(category != 'truncating')

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes_v2.rds') %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(rna_expression = value) %>% 
  dplyr::filter(!is.na(rna_expression))

transcriptomics_data = classify_mutations(transcriptomics_data) %>% 
  dplyr::filter(category != 'truncating')
