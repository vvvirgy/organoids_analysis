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
    xlab(paste0(what, ' coefficients')) + 
    theme(legend.position = 'bottom')
                    
}

plot_predictor_higher_impact = function(x, palette = wes_palette("GrandBudapest2"), pth = 0.05) {
  x %>% 
    dplyr::filter(term != '(Intercept)') %>% 
    dplyr::select(term, gene, estimate, p.value, method) %>% 
    filter(p.value <= pth) %>% 
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

classify_betas = function(x, pth = 0.05) {
  x %>% 
    dplyr::filter(term != '(Intercept)') %>% 
    dplyr::select(term, gene, estimate, p.value, method) %>% 
    dplyr::filter(p.value <= pth) %>%
    tidyr::pivot_wider(names_from = term, values_from = c(estimate, p.value), values_fill = 0) %>% 
    mutate(class = 
             case_when(
               (estimate_ploidy_diff > 0 & estimate_mutation_multiplicity > 0) ~ 'Both estimators positive', 
               (estimate_ploidy_diff < 0 & estimate_mutation_multiplicity < 0) ~ 'Both estimators negative',
               (estimate_ploidy_diff > 0 & estimate_mutation_multiplicity < 0) ~ 'Ploidy positive, multiplicity negative',
               (estimate_ploidy_diff < 0 & estimate_mutation_multiplicity > 0) ~ 'Ploidy negative, multiplicity positive',
               (estimate_ploidy_diff > 0 & estimate_mutation_multiplicity == 0) ~ 'Ploidy positive, multiplicity not significant',
               (estimate_ploidy_diff < 0 & estimate_mutation_multiplicity == 0) ~ 'Ploidy negative, multiplicity not significant',
               (estimate_ploidy_diff == 0 & estimate_mutation_multiplicity > 0) ~ 'Ploidy not significant, multiplicity positive',
               (estimate_ploidy_diff == 0 & estimate_mutation_multiplicity < 0) ~ 'Ploidy not significant, multiplicity negative',
             )) 
}

cols = setNames(nm = c("Ploidy positive, multiplicity not significant", "Ploidy not significant, multiplicity positive", 
                       "Ploidy not significant, multiplicity negative","Ploidy positive, multiplicity negative",        
                       "Ploidy negative, multiplicity not significant", "Both estimators positive",                     
                       "Ploidy negative, multiplicity positive", "Both estimators negative"),
                object = rcartocolor::carto_pal(n = 9, name = 'Prism')[1:8])

betas_association = function(x, cols = cols, pth = 0.05) {
  
  classify_betas(x, pth = pth) %>%
    # filter(p.value_ploidy_diff <= pth | p.value_mutation_multiplicity <= pth) %>%
    ggplot(aes(estimate_mutation_multiplicity, estimate_ploidy_diff, color = class)) + 
    geom_point(alpha = 0.7) + 
    theme_bw() + 
    xlab('Mutation multiplicity') +
    ylab('Ploidy diff') + 
    theme(legend.position = 'bottom') + 
    scale_color_manual(values = cols) + 
    guides(colour = guide_legend(nrow = 4))
}

# kmeans clustering

kmeans_estimates_by_group = function(x, seed = 12345, pth = 0.05) {
  # cluster univariate the genes per each coefficient
  x = x %>% 
    filter(term != '(Intercept)') 
  
  significance_class = x %>% 
    dplyr::select(gene, estimate, p.value, term) %>%
    pivot_wider(names_from = term, values_from = c(estimate, p.value)) %>% 
    # filter(p.value_ploidy_diff <= pth & p.value_mutation_multiplicity <= pth)
    mutate(Significance = case_when(
      (p.value_ploidy_diff <= 0.05 & p.value_mutation_multiplicity <= 0.05) ~ 'both', #) %>% 
      (p.value_ploidy_diff <= 0.05 & p.value_mutation_multiplicity > 0.05 | p.value_ploidy_diff <= 0.05 & is.na(p.value_mutation_multiplicity) ) ~ 'ploidy_diff',
      (p.value_ploidy_diff > 0.05 & p.value_mutation_multiplicity < 0.05 | is.na(p.value_ploidy_diff)  & p.value_mutation_multiplicity <= 0.05) ~ 'multiplicity', 
      (p.value_ploidy_diff > 0.05 & p.value_mutation_multiplicity > 0.05) ~ 'None', 
      .default = 'None'
    )) %>% 
    filter(Significance != 'None') %>% 
    group_by(Significance) %>% 
    group_split()
    
  lapply(significance_class, function(s) {
    
    cls = s$Significance %>% unique

    if(cls == 'both') {
      df = s %>% 
        dplyr::select(gene, starts_with('estimate')) %>% 
        tibble::column_to_rownames('gene')
    } else {
      
      df = s %>% 
        dplyr::select(gene,(ends_with(cls) & starts_with('estimate'))) %>% 
        tibble::column_to_rownames('gene')
      
    }
    
    set.seed(seed)
    
    sil_plot = fviz_nbclust(df, FUNcluster = kmeans, method = "silhouette")
    k = sil_plot@layers$geom_vline$data$xintercept
    
    clustering = kmeans(df, k, nstart = 10)
    
    clusters = tibble(
      gene = names(clustering$cluster), 
      cluster = unname(clustering$cluster)
    )
    
    return(list(kmeans = clustering, sil_plot = sil_plot, clusters = clusters, data = s, class = cls))
  })
}  

kmeans_estimates = function(x, seed = 12345, pth = 0.05) {
  # cluster univariate the genes per each coefficient
  x = x %>%
    filter(term != '(Intercept)') %>%
    dplyr::filter(p.value <= pth) %>%
    select(gene, term, estimate) %>%
    tidyr::pivot_wider(names_from = term, values_from = estimate, values_fill = 0) %>%
    tibble::column_to_rownames('gene')
  
  set.seed(seed)
  
  sil_plot = fviz_nbclust(x, FUNcluster = kmeans, method = "silhouette")
  k = sil_plot@layers$geom_vline$data$xintercept
  
  clustering = kmeans(x, k, nstart = 10)
  
  clusters = tibble(
    gene = names(clustering$cluster), 
    cluster = unname(clustering$cluster)
  )
  
  return(list(kmeans = clustering, sil_plot = sil_plot, clusters = clusters, data = x))
}

plot_clusters_1D = function(x) {
  
  lapply(x, function(s){
    s$data %>% 
      full_join(., s$clusters, by = 'gene') %>%
      mutate(cluster = paste0('C', cluster)) %>% 
      ggplot(aes(estimate, fill = cluster)) + 
      geom_histogram(binwidth = 0.1) + 
      theme_bw() + 
      ggsci::scale_fill_npg()
  })
}

plot_clusters_2D = function(x) {
    x$data %>% 
      tibble::rownames_to_column('gene') %>% 
      full_join(., x$clusters, by = 'gene') %>%
      mutate(cluster = paste0('C', cluster)) %>% 
      ggplot(aes(ploidy_diff, mutation_multiplicity, color = cluster)) + 
      geom_point() + 
      theme_bw() + 
      ggsci::scale_color_npg()
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
