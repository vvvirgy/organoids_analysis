### utils

# Computing single allele expression and statistics --------------------------------------------------------------------------------

compute_expr_single_allele = function(x, column, which = 'Wild-type', use_2n = TRUE) {
  if(use_2n == TRUE) {
    x = x %>% 
      filter(mutation_status == which) %>% 
      # filter(karyotype == '1:1')
      filter(tot_cna == 2)
  } else {
    x = x %>% 
      filter(mutation_status == which)
  }
  x %>% 
    mutate(observed_expr = get(column)/tot_cna)%>% 
    group_by(hgnc_symbol) %>% 
    mutate(mean_wt_expr_observed = mean(observed_expr), sd_expr = sd(observed_expr))
    # summarise(mean_wt_expr_observed = mean(single_allele_expr), sd_expr = sd(single_allele_expr)) 
}

compute_expr_mut_allele = function(x, use_2n = TRUE, column) {
  
  if(use_2n == TRUE) {
    x = x %>% 
      filter(tot_cna != 2 | mutation_status == 'Mutated')
      # filter(karyotype != '1:1') 
  } else {
    x = x %>% 
      filter(mutation_status == 'Mutated')
  }
  x %>% 
    mutate(observed_expr = get(column)/tot_cna)
}

compute_expectation_expression = function(x, ploidy_seq = seq(2:6)) {
  x %>% 
    bind_cols(., 
              map_dfc(ploidy_seq, ~ x$mean_wt_expr_observed * .x) %>%
                set_names(paste0("n", ploidy_seq))
    ) %>% 
    pivot_longer(cols = starts_with('n'), names_to = 'ploidy', values_to = 'estimate_expression') %>% 
    mutate(ploidy = as.numeric(gsub('n', '', ploidy)))
}

join_allele_expr = function(x, column, which = 'Wild-type', use_2n = TRUE) {
  wt = compute_expr_single_allele(x, column = column, which = which, use_2n = use_2n)
  mut = compute_expr_mut_allele(x, column = column, use_2n = use_2n)
  
  p_seq = range(mut$tot_cna %>% unique)
  
  expected_expression = compute_expectation_expression(wt, ploidy_seq = seq(p_seq[1], p_seq[2])) %>% 
    dplyr::select(chr, hgnc_symbol, ploidy, estimate_expression) %>% 
    distinct()
  
  x = wt %>%
    dplyr::select(hgnc_symbol,mean_wt_expr_observed, sd_expr) %>%
    distinct() %>%
    full_join(., mut) %>%
    bind_rows(., wt) %>% 
    left_join(., expected_expression, by = join_by('hgnc_symbol' == 'hgnc_symbol', 
                                                   'tot_cna' == 'ploidy', 
                                                   'chr' == 'chr'
                                                   )) %>% 
    mutate(alteration = case_when(
      (tot_cna == 2 & mutation_status != 'Wild-type') ~ 'Mutated, no CNA',
      (tot_cna != 2 & mutation_status == 'Wild-type') ~ 'Wild-type, CNA',
      (tot_cna != 2 & mutation_status != 'Wild-type') ~ 'Mutated, CNA', 
      (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'Wild-type, no CNA' 
    )) %>% 
    mutate(alteration_classes = case_when(
      # (tot_cna == 2 & mutation_status != 'Wild-type') ~ 'Mutated, no CNA',
      # (tot_cna != 2 & mutation_status == 'Wild-type') ~ 'Wild-type, CNA',
      # (tot_cna != 2 & mutation_status != 'Wild-type') ~ 'Mutated, CNA', 
      (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'No alteration', 
      .default = 'Alteration'
    )) 

  return(x)
}



# plotting function
plot_single_allele_expr = function(df, gene) {
  df %>%  
    dplyr::filter(hgnc_symbol == gene) %>%
    ggplot() + 
    geom_segment(aes(
      y = sample,
      x = observed_expr, 
      xend = mean_wt_expr_observed, 
      colour = alteration
    )) +
    geom_rect(aes(xmin = mean_wt_expr_observed-sd_expr, 
                  xmax = mean_wt_expr_observed+sd_expr, 
                  ymax = Inf, 
                  ymin = 0), 
              alpha = 0.1, 
              fill = 'paleturquoise') +
    geom_point(aes(
      y = sample, 
      x = observed_expr, 
      colour = alteration
    )) +
    ggsci::scale_color_futurama() +
    geom_vline(aes(xintercept = mean_wt_expr_observed), linetype = 'dashed', show.legend = F, colour = 'firebrick') + 
    theme_bw() + 
    facet_wrap(~alteration, scales = 'free_y') + 
    ggtitle(gene) + 
    xlab('Observed expression/ploidy') + 
    ylab('')
}



# Expression testing ------------------------------------------------------

# t test for the expression


# test_expression_groups = function(x) {
#   x %>% 
#     summarise(
#       n_WT2N = sum(alteration_classes == "No alteration"),
#       n_Altered = sum(alteration_classes == "Alteration"),
#       test = if (n_WT2N >= 2 & n_Altered >= 2) {
#         list(t.test(observed_expr ~ alteration_classes, var.equal = FALSE))
#       } else {
#         list(NA_real_)
#       }, 
#       mean_WT2N = mean(observed_expr[alteration_classes == "WT_2N"], na.rm = TRUE),
#       mean_Altered = mean(observed_expr[alteration_classes == "Altered"], na.rm = TRUE)
#     ) %>% 
#     mutate(tidy = map(test, broom::tidy)) %>% 
#     unnest(tidy) %>%
#     mutate(FDR = p.adjust(p.value, method = "fdr")) %>% 
#     select(-c(test, x))
# }

# trying to include also bootstrap - testing
# test_belonging_alterations = function(x, 
#                                       pth = .05, 
#                                       n_bootstrap = 10000) {
#   x %>% 
#     group_by(hgnc_symbol) %>% 
#     group_modify(~ {
#       
#       ref <- .x %>% filter(alteration_classes == "No alteration") %>% pull(observed_expr)
#       alt <- .x %>% filter(alteration_classes == 'Alteration')
#       
#       print(ref)
#       
#       if (length(ref) == 1 | nrow(alt) == 0) {
#         return(tibble())
#       } else {
#         
#         if(n_bootstrap > 0) {
#           
#           bootstrap_data = (mosaic::do(n_bootstrap) * mosaic::resample(ref))
#           sigma = apply(bootstrap_data,1,sd) %>% mean
#           mu = apply(bootstrap_data,1,mean) %>% mean
#           ci_int = abs(pth/2+ c(0, -1))
#           ci = quantile(apply(bootstrap_data,1,mean), probs = ci_int)
#           
#         } else {
#           mu <- mean(ref, na.rm = TRUE)
#           
#           if (length(ref) == 2) {
#             sigma <- abs(ref[1]-ref[2])/2
#           }
#           
#           sigma <- sd(ref, na.rm = TRUE)
#           ci = NULL
#         }
#         
#         alt %>%
#           mutate(
#             sd_not_alt = sigma, 
#             mu_not_alt = mu,
#             z = (observed_expr - mu) / sigma,
#             p_value = 2 * pnorm(-abs(z)),
#             significant = p_value <= pth, 
#             CI_lower = min(ci), 
#             CI_upper = max(ci), 
#             CI = pth
#           ) 
#       }
#       
#       
#     })
# }

# classification
classify_elements = function(x) {
  x %>%
    mutate(
      side = case_when(
        (significant == TRUE & z < 0) ~ 'lower',
        (significant == TRUE & z > 0) ~ 'higher',
        .default = NA
      )
    ) %>% 
    mutate(cls_dosage = 
             case_when(
               (tot_cna != 2 & significant == FALSE) ~ 'Dose_sensitive',
               
               (tot_cna == 1 & side == 'higher') ~ 'Dosage_compensating',
               (tot_cna == 1 & side == 'lower') ~ 'Haploinsufficient',
               
               (tot_cna == 2 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
               (tot_cna == 2 & side == 'higher' & mutation_status == 'Mutated') ~ 'Mutation_sensitive',
               (tot_cna == 2 & significant == FALSE & mutation_status == 'Mutated') ~ 'Mutation_insensitive',
               # (tot_cna == 2 & significant == FALSE & mutation_status == 'Mutated') ~ 'Mutation_insensitive',
               
               (tot_cna == 3 & side == 'lower') ~ 'Dosage_insensitive',
               (tot_cna == 3 & side == 'higher') ~ 'Super_dosage_sensitive',
               (tot_cna == 4 & side == 'lower') ~ 'Dosage_insensitive',
               (tot_cna == 4 & side == 'higher') ~ 'Super_dosage_sensitive'

               # (tot_cna == 2 & significant == TRUE & mutation_status == 'Mutated') ~ 'Mutation_sensitive',
               # significant == FALSE ~ 'Dosage_sensitive', 
               # (tot_cna == 1 & side == 'higher') ~ 'Dosage_compensating',
               # (tot_cna == 1 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
               # (tot_cna == 2 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
               # (tot_cna == 2 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
               
             )
           
           ) 
}

# Summary statistics getters -------------------------------------------------------------------------------------------


# get propotions
get_props = function(x, group, filter = TRUE) {
  
  if(filter == TRUE) {
    x = x %>% 
      filter(!is.na(cls_dosage))
  }
  
  x %>%
    mutate(cls_dosage = ifelse(is.na(cls_dosage), 'Not_assigned', cls_dosage)) %>% 
    # mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
    group_by(pick(group)) %>% 
    count(cls_dosage) %>% 
    mutate(tot = sum(n)) %>% 
    mutate(prop = n/tot)
}

# bind together the rna and prot classification and check if the classification is concordant among the two 
merging_omics = function(rna_cls, prot_cls) {
  
  rna_cls = rna_cls %>% 
    dplyr::select(hgnc_symbol, mut_consequence, sample, tot_cna, mutation_status, cls_dosage, alteration)
  prot_cls = prot_cls %>% 
    dplyr::select(hgnc_symbol, mut_consequence, sample, tot_cna, mutation_status, cls_dosage, alteration)
  
  x = full_join(rna_cls, prot_cls, by = join_by(
    'hgnc_symbol' == 'hgnc_symbol', 
    'mut_consequence' == 'mut_consequence', 
    'sample' == 'sample', 
    'tot_cna' == 'tot_cna', 
    'mutation_status' == 'mutation_status', 
    'alteration' == 'alteration'
  ), 
  suffix = c('_RNA', '_Prot'))
  
  x = x %>% 
    mutate(cls_comparison = case_when(
      cls_dosage_RNA == cls_dosage_Prot ~ 'Same class', 
      cls_dosage_RNA != cls_dosage_Prot ~ 'Different classes',
      (is.na(cls_dosage_RNA) & !is.na(cls_dosage_Prot)) ~ 'Only protein', 
      (!is.na(cls_dosage_RNA) & is.na(cls_dosage_Prot)) ~ 'Only RNA'
    ))
}

# define the prevalent class for each gene

get_prevalent_class = function(x, filter = TRUE) {
  if(filter == TRUE) {
    x = x %>% 
      filter(!is.na(cls_dosage))}
  x %>% 
    filter(cls_dosage != 'Control') %>% 
    group_by(hgnc_symbol) %>% 
    count(cls_dosage) %>% 
    mutate(prop = n/sum(n)) 
}

# plotting -----------------
# plot expression 
dosage_colors = setNames(
  nm = c('Dosage_compensating', 'Dosage_sensitive', 'Dosage_insensitive', 'Haploinsufficient', 'Mutation_sensitive', 'Mutation_insensitive', 'Control', 'Not_assigned'), 
  object = c('#BF092F', '#001BB7', '#31694E', '#3B9797', '#6B3F69',  '#E195AB', '#84994F', '#CBCBCB')
)

plot_expr_gene_by_cls = function(x, gene, which, ploidy_seq = seq(2:6)) {
  df = x %>% 
    filter(hgnc_symbol == gene) %>% 
    # mutate(tot_cna = factor(tot_cna)) %>% 
    filter(!is.na(cls_dosage)) %>% 
    filter(!is.na(mutation_status)) 
  max_expr = max(df$expression) + 3
  min_expr = min(df$expression) - 2
  
  conf_int = df %>% 
    select(mu_not_alt, starts_with('CI')) %>% 
    filter(!is.na(mu_not_alt)) %>% 
    distinct() %>% 
    slice(rep(1,length(ploidy_seq))) %>% 
    mutate(ploidy = ploidy_seq, 
           CI = 1-CI) %>% 
    mutate(theorical_expression = mu_not_alt*ploidy, 
           CI_lower = CI_lower*ploidy, 
           CI_upper = CI_upper*ploidy, 
           theoric = 'Theoretical expression', 
           CI = as.character(CI))
   
  ci_col = setNames(nm = conf_int$CI %>% unique, object = 'lightskyblue')
  
  p = df %>% 
    # mutate(theoric = 'Theoretical expression') %>% 
    ggplot(aes(x = tot_cna, 
               y = expression, 
               color = cls_dosage, 
               shape = mutation_status)) + 
    geom_ribbon(data = conf_int, 
                aes(
                  x = ploidy,
                  ymin = CI_lower, 
                  ymax = CI_upper, 
                  fill = CI
                ), 
                inherit.aes = F, 
                alpha = 0.3, 
                # fill = 'lightskyblue', 
                show.legend = T) + 
    geom_point(size = 3) + 
    geom_line(data = conf_int,
                aes(y = theorical_expression, 
                  x = ploidy,
                  linetype = theoric
                  ), 
              inherit.aes = F, 
              # linetype = 'dashed',
              alpha = .6) + 
    theme_bw() + 
    # geom_hline(data = conf_int, 
    #            aes(yintercept = CI_lower))
    labs(y = paste0(which, ' expression'), 
         x = 'Ploidy') + 
    scale_color_manual(values = dosage_colors) + 
    scale_shape_manual(values = c('Wild-type' = 16, 'Mutated' = 17)) +
    scale_linetype_manual(values = c('Theoretical expression' = 'dashed')) +
    scale_fill_manual(values = ci_col) + 
    ggtitle(gene) + 
    coord_cartesian(ylim = c(min_expr, #min(df$estimate_expression), 
                             max_expr)) + 
    guides(color = guide_legend(title = 'Dosage classes', 
                                ncol = 3, order = 1, 
                                title.position = 'top', 
                                override.aes = list(fill = "white")), 
           shape = guide_legend(title = 'Mutation status', 
                                nrow = 2, order = 2,
                                title.position = 'top', 
                                override.aes = list(fill = "white")), 
           linetype = guide_legend(title = '', order = 3, 
                                   override.aes = list(fill = "white")), 
           fill = guide_legend(title = 'Confidence interval', 
                               ncol = 2)) + 
    theme(legend.position = 'bottom')
  return(p)
}

plot_cls_by_sample = function(x, dosage_colors) {
  x %>% 
    filter(!is.na(sample)) %>% 
    pivot_longer(cols = c(prop_prot, prop_rna), names_to = 'Assay', values_to = 'prop') %>% 
    mutate(Assay = gsub('prop_', '', Assay)) %>% 
    mutate(Assay = factor(Assay, levels = c('rna', 'prot'))) %>% 
    ggplot(aes(
      y = sample, x = prop, fill = cls_dosage
    )) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    theme_bw() + 
    facet_wrap(~Assay, 
               labeller = labeller(
                 Assay = setNames(
                   c('RNA', 'Protein'), 
                   c('rna', 'prot')))) + 
    scale_fill_manual(values = dosage_colors) + 
    labs(y = 'Sample', 
         x = 'Genes proportions') +
    guides(fill = guide_legend(title = 'Dosage classes', 
                                ncol = 4, order = 1, 
                                title.position = 'top')) +
    theme(legend.position = 'bottom', legend.title = element_text(hjust = 0.5))
    
}

# heatmap and oncoprint
create_annotation = function(x, ann_colors, position) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position,
                                    annotation_legend_param = list(nrow = 3, width = 12, by_row = T))
  
}

plot_heatmap = function(rna, 
                        prot, 
                        ann, 
                        genes, 
                        cols) {
  
  class_merge = merging_omics(prot_cls = prot, rna_cls = rna)
  x = class_merge %>%
    dplyr::select(hgnc_symbol, sample, cls_dosage_RNA, cls_dosage_Prot)
  
  levels_cls = c('Dosage_compensating', 'Dosage_sensitive', 'Dosage_insensitive', 
                 'Haploinsufficient', 'Mutation_sensitive', 'Mutation_insensitive', 
                 'Control', 'Not_assigned')
  
  # wide the data independently and match them after
  
  x_wide = x %>% 
    mutate(cls_dosage_RNA = ifelse(is.na(cls_dosage_RNA), 'Not_assigned', cls_dosage_RNA)) %>%
    mutate(cls_dosage_Prot = ifelse(is.na(cls_dosage_Prot), 'Not_assigned', cls_dosage_Prot)) %>%
    mutate(cls_dosage_RNA = factor(cls_dosage_RNA, levels = levels_cls))  %>%
    mutate(cls_dosage_Prot = factor(cls_dosage_Prot, levels = levels_cls))  %>%
    mutate(Protein = as.numeric(cls_dosage_Prot),
           RNA = as.numeric(cls_dosage_RNA)) %>%
    dplyr::select(-c(starts_with('cls_dosage_'))) %>% 
    dplyr::filter(hgnc_symbol %in% genes) %>% 
    pivot_wider(names_from = sample, values_from = c(Protein, RNA), values_fill = 0) %>% 
    tibble::column_to_rownames('hgnc_symbol')
  
  Heatmap(x_wide)
  
  
}



long2wide = function(x, lvs, suffix) {
  x %>% 
    dplyr::select(hgnc_symbol, sample, cls_dosage) %>%  
    mutate(cls_dosage = ifelse(is.na(cls_dosage), 'Not_assigned', cls_dosage)) %>%
    mutate(cls_dosage = factor(cls_dosage, levels = levels_cls))  %>%
    mutate(class_num = as.numeric(cls_dosage)) %>%
    mutate(sample = paste(sample, suffix, sep = '_')) %>% 
    dplyr::select(-cls_dosage) %>% 
    pivot_wider(names_from = sample, values_from = class_num) %>% 
    tibble::column_to_rownames('hgnc_symbol')
}


# plot_expression_classes = function(x, gene, which) {
#   
#   x = test_expr_rna %>% 
#     filter(hgnc_symbol == gene)
#   
#   ref = x %>% 
#     filter(cls_dosage == 'Control') 
#   alt = x %>% 
#     filter(cls_dosage != 'Control') 
#   
#   ggplot() + 
#     geom_density(data = ref, aes(observed_expr, fill = cls_dosage), 
#                  alpha = .8, 
#                  color = NA) + 
#     geom_histogram(data = alt, aes(observed_expr, fill = cls_dosage)) + 
#     scale_fill_manual(values = dosage_colors) + 
#     xlab(paste0('Observed per allele expression ', which)) + 
#     ggtitle(gene) + 
#     theme_bw() + 
#     theme(legend.position = 'bottom')
# }

# clustering ------------------------
# K means clustering

kmeans_estimates = function(x, seed = 12345, pth = 0.05) {
  # cluster univariate the genes per each coefficient
  x = x %>%
    select(hgnc_symbol, observed_expr, sample) %>%
    group_by(hgnc_symbol) %>% 
    # tibble::column_to_rownames('hgnc_symbol') %>% 
    group_split()
  names(x) = lapply(x, function(s) s$hgnc_symbol %>% unique)
  
  # tidyr::pivot_wider(names_from = term, values_from = estimate, values_fill = 0) %>%
  
  
  
  k_clus = lapply(x, function(s) {
    
    set.seed(seed)
    
    tryCatch({
      
      data = s %>% 
        tibble::column_to_rownames('sample') %>% 
        select(observed_expr)
      sil_plot = fviz_nbclust(data, FUNcluster = kmeans, method = "silhouette")
      k = sil_plot@layers$geom_vline$data$xintercept
      
      clustering = kmeans(data, k, nstart = 10)
      
      clusters = tibble(
        sample = names(clustering$cluster), 
        cluster = unname(clustering$cluster)
      )
      
      list(kmeans = clustering, sil_plot = sil_plot, clusters = clusters, data = s)
    }, 
    error = function(e) {NA})
  })
  
  return(k_clus)
  
}


# test_expression_groups = function(x, column = 'rna_expression') {
#   df = x %>% 
#     mutate(alteration = case_when(
#       (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'Not_altered',
#       .default = 'Altered'
#     )) 
#     # pivot_wider(names_from = alteration, values_from = any_of(column)) %>% 
#     # group_by(hgnc_symbol, alteration)
#     # filter(hgnc_symbol == 'TP53') %>% 
#   
#   # remove genes that do not have one of the two classes
#   genes2remove = df %>%   
#     mutate(alteration = factor(alteration, levels = c('Altered', 'Not_altered'))) %>% 
#     group_by(hgnc_symbol) %>% 
#     # filter(n() <= 4) %>%
#     group_by(hgnc_symbol, alteration) %>% 
#     count() %>% 
#     ungroup() %>%
#     group_by(hgnc_symbol) %>%
#     filter(n <= 1) %>% 
#     pull(hgnc_symbol) %>% 
#     unique
#   
#   df %>% 
#     ungroup() %>% 
#     group_by(hgnc_symbol) %>% 
#     filter(!hgnc_symbol %in% genes2remove) %>%
#     # rstatix::t_test(rna_expression ~ alteration)
#     summarise(tidy_t = list(t.test(rna_expression ~ alteration, data = cur_data()))) %>%
#     mutate(result = broom::tidy(tidy_t[[1]])) %>%
#     unnest(result)
#     # group_nest(hgnc_symbol) %>%
#     
#   
# }
# 
# test = x %>% 
#   mutate(alteration = case_when(
#     (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'Not_altered',
#     .default = 'Altered'
#   )) %>% 
#   # pivot_wider(names_from = alteration, values_from = any_of(column)) %>% 
#   group_by(hgnc_symbol) %>%
#   filter(hgnc_symbol == 'SAMD11')
# 
# test_alt = test %>% filter(alteration == 'Altered') %>% pull(rna_expression)
# test_wt = test %>% filter(alteration == 'Not_altered') %>% pull(rna_expression)
# 
# rstatix::t_test(data = test, )
