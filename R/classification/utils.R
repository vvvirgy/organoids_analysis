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
    mutate(observed_expr = get(column)/tot_cna) %>% 
    mutate(tot_cna_h1 = tot_cna - multiplicity) %>% 
    mutate(observed_expr_h1 = get(column)/tot_cna_h1)
    # mutate(tot_cna_h1 = #case_when(
    #   # IMPACT %in% c('HIGH', 'MODERATE') ~ 
    #     tot_cna-multiplicity
    # ) %>% 
    # mutate(
    #   observed_expr_h1 = #case_when(
    #     # IMPACT %in% c('HIGH', 'MODERATE') ~ 
    #     get(column)/tot_cna_h1, .default = NA
      # ))
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
# classify_elements = function(x) {
#   x %>%
#     mutate(
#       side = case_when(
#         (significant == TRUE & z < 0) ~ 'lower',
#         (significant == TRUE & z > 0) ~ 'higher',
#         .default = NA
#       )
#     ) %>% 
#     mutate(cls_dosage = 
#              case_when(
#                (tot_cna != 2 & significant == FALSE) ~ 'Dose_sensitive',
#                
#                (tot_cna == 1 & side == 'higher') ~ 'Dosage_compensating',
#                (tot_cna == 1 & side == 'lower') ~ 'Haploinsufficient',
#                
#                (tot_cna == 2 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
#                (tot_cna == 2 & side == 'higher' & mutation_status == 'Mutated') ~ 'Mutation_sensitive',
#                (tot_cna == 2 & significant == FALSE & mutation_status == 'Mutated') ~ 'Mutation_insensitive',
#                # (tot_cna == 2 & significant == FALSE & mutation_status == 'Mutated') ~ 'Mutation_insensitive',
#                
#                (tot_cna == 3 & side == 'lower') ~ 'Dosage_insensitive',
#                (tot_cna == 3 & side == 'higher') ~ 'Super_dosage_sensitive',
#                (tot_cna == 4 & side == 'lower') ~ 'Dosage_insensitive',
#                (tot_cna == 4 & side == 'higher') ~ 'Super_dosage_sensitive'
# 
#                # (tot_cna == 2 & significant == TRUE & mutation_status == 'Mutated') ~ 'Mutation_sensitive',
#                # significant == FALSE ~ 'Dosage_sensitive', 
#                # (tot_cna == 1 & side == 'higher') ~ 'Dosage_compensating',
#                # (tot_cna == 1 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
#                # (tot_cna == 2 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
#                # (tot_cna == 2 & side == 'lower' & mutation_status == 'Mutated') ~ 'Haploinsufficient',
#                
#              )
#            
#            ) 
# }

classify_elements = function(x, pth = 0.05) {
  
  # first checking mutations w/out high impact and wild type --> classified wrt the theoretical rect only
  group1 = x %>% 
    filter(is.na(p_value_h1))
  
  group1 = group1 %>% 
    mutate(class = case_when(
      (tot_cna != 2 & p_value > pth) ~ 'Dose sensitive (no mutations/low impact mutation)', 
      (tot_cna == 2 & mutation_status == 'Mutated' & p_value > pth) ~ 'Dose sensitive (no mutations/low impact mutation)', 
      (tot_cna != 2 & p_value <= pth & z > 0) ~ 'Enhanced expression (no mutations/low impact mutations)',
      (tot_cna != 2 & p_value <= pth & z < 0) ~ 'Reduced expression (no mutations/low impact mutations)',
      (mutation_status == 'Mutated' & p_value <= pth & z > 0) ~ 'Enhanced expression (no mutations/low impact mutations)',
      (mutation_status == 'Mutated' & p_value <= pth & z < 0) ~ 'Reduced expression (no mutations/low impact mutations)',
      (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'Control', 
      .default = 'Not classified'
    ))
  
  group2 = x %>% 
    filter(!is.na(p_value_h1))
  
  group2 = group2 %>% 
    mutate(class = case_when(
      (p_value > pth & p_value_h1 <= pth) ~ 'Enhanced expression (with impactful mutation)', 
      (p_value <= pth & p_value_h1 > pth) ~ 'Dose sensitive (with impactful mutation)', 
      (p_value <= pth & p_value_h1 <= pth & z_h1 > 0 & z > 0) ~ 'Super-enhanced expression (with impactful mutation)', 
      (p_value <= pth & p_value_h1 <= pth & z_h1 < 0 & z < 0) ~ 'Higly reduced expression (with impactful mutation)', 
      (p_value <= pth & p_value_h1 <= pth & z_h1 > 0 & z < 0) ~ 'Buffering expression (with impactful mutation)', 
      .default = 'Not classified'
    ))
  
  bind_rows(group1, group2) %>% 
    arrange(hgnc_symbol)
  
}

# Summary statistics getters -------------------------------------------------------------------------------------------


# get propotions
get_props = function(x, group, filter = TRUE) {
  
  if(filter == TRUE) {
    x = x %>% 
      filter(!is.na(class))
  }
  
  x %>%
    mutate(cls_dosage = ifelse(is.na(class), 'Not_assigned', class)) %>% 
    # mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
    group_by(pick(group)) %>% 
    count(class) %>% 
    mutate(tot = sum(n)) %>% 
    mutate(prop = n/tot)
}

# bind together the rna and prot classification and check if the classification is concordant among the two 
merging_omics = function(rna_cls, prot_cls) {
  
  rna_cls = rna_cls %>% 
    dplyr::select(hgnc_symbol, mut_consequence, sample, tot_cna, mutation_status, class, alteration)
  prot_cls = prot_cls %>% 
    dplyr::select(hgnc_symbol, mut_consequence, sample, tot_cna, mutation_status, class, alteration)
  
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
      class_RNA == class_Prot ~ 'Same class', 
      class_RNA != class_Prot ~ 'Different classes',
      (is.na(class_RNA) & !is.na(class_Prot)) ~ 'Only protein', 
      (!is.na(class_RNA) & is.na(class_Prot)) ~ 'Only RNA'
    ))
}

# define the prevalent class for each gene

get_prevalent_class = function(x, filter = TRUE) {
  if(filter == TRUE) {
    x = x %>% 
      filter(!is.na(class))}
  x %>% 
    filter(class != 'Control') %>% 
    group_by(hgnc_symbol) %>% 
    count(class) %>% 
    mutate(prop = n/sum(n)) 
}

get_prevalent_class_for_gene_type = function(x, filter = TRUE) {
  if(filter == TRUE) {
    x = x %>% 
      filter(!is.na(class))}
  x %>% 
    filter(class != 'Control') %>% 
    group_by(hgnc_symbol) %>% 
    count(class) %>% 
    mutate(prop = n/sum(n)) 
}

# plotting -----------------
# plot expression 

## plot utilities ----

# ggplot theme
light_theme = function() {
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = 'bottom', 
        legend.title = element_text(hjust = 0.5, 
                                    vjust = 1))
}

# alteration colors
alt_colors = setNames(
  # object = c('#4A70A9', '#658C58', '#B95E82', '#8D5F8C'),
  c('#476EAE', '#48B3AF', '#A7E399', '#F6FF99'), 
  nm = c("Wild-type, CNA", 'Mutated, CNA', 'Mutated, no CNA', 'Wild-type, no CNA')
)

# assay colors
omics_cols = setNames(object = c('steelblue', 'goldenrod'), c('RNA', 'Protein'))

# classes 
# levels_cls = c('Not_assigned', 'Dosage_compensating', 'Dosage_sensitive', 'Dosage_insensitive', 
#                'Haploinsufficient', 'Mutation_sensitive', 'Mutation_insensitive', 
#                'Control')

dosage_colors = setNames(
  nm = c(
    'Control',
    'Enhanced expression (with impactful mutation)', 
    'Super-enhanced expression (with impactful mutation)',
    'Enhanced expression (no mutations/low impact mutations)', 
    'Dose sensitive (no mutations/low impact mutation)', 
    'Dose sensitive (with impactful mutation)', 
    'Reduced expression (no mutations/low impact mutations)', 
    'Higly reduced expression (with impactful mutation)', 
    'Buffering expression (with impactful mutation)', 
    'Not classified'
    ), 
  # object = c(
  #   '#BF092F', 
  #   '#5459AC', 
  #   '#31694E', 
  #   '#3B9797', 
  #   '#6B3F69', 
  #   '#E195AB', 
  #   # '#84994F',
  #   '#809D3C',
  #   '#CBCBCB')
  object = c(
    '#5D8736', 
    '#FF2929',
    '#8C1007', 
    '#E82561', 
    '#FAB12F', 
    '#EA7300',
    '#3B9797', 
    '#28518A', 
    '#5D2F77', 
    # '#809D3C',
    '#CBCBCB'
    )
)

levels_cls = names(dosage_colors)
# 
# cols = c(
#   '#BF092F', 
#   '#28518A', 
#   '#1E5128', 
#   '#3B9797', 
#   '#5D2F77', 
#   '#C95792', 
#   # '#809D3C',
#   '#5D8736',
#   '#CBCBCB')
# scales::show_col(cols)

## plot functions ----

### by gene ----

plot_expr_gene_by_cls = function(x, gene,ploidy_seq = seq(2:6), color_by_classes = TRUE, filter_wt = TRUE) {
  df = x %>% 
    filter(hgnc_symbol == gene) %>% 
    # mutate(tot_cna = factor(tot_cna)) %>% 
    # filter(!is.na(cls_dosage)) %>% 
    filter(!is.na(mutation_status)) 
  
  if(nrow(df) > 2) {
    
    max_expr = max(df$expression) + 3
    min_expr = min(df$expression) - 2
    
    conf_int = df %>% 
      group_by(Assay) %>% 
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
    
    if(filter_wt == TRUE) {
      df = df %>% 
        filter(mutation_status == 'Wild-type')
    }
    
    if (color_by_classes == TRUE) {
      p = df %>% 
        filter(!is.na(class)) %>% 
        # mutate(theoric = 'Theoretical expression') %>% 
        ggplot(aes(x = tot_cna, 
                   y = expression, 
                   color = class, 
                   shape = mutation_status)) + 
        geom_point(size = 3) + 
        scale_color_manual(values = dosage_colors) + 
        guides(color = guide_legend(title = 'Dosage classes', 
                                    order = 1, 
                                    ncol = 1, 
                                    title.position = 'top', 
                                    override.aes = list(fill = "white")))
      
    } else {
      p = df %>% 
        # mutate(theoric = 'Theoretical expression') %>% 
        ggplot(aes(x = tot_cna, 
                   y = expression, 
                   shape = mutation_status)) + 
        geom_point(size = 3)
    }
    
    title_p = paste0(gene, ' (', df$CGC_role_PANCANCER, ')')
    
    p = p + 
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
      
      geom_line(data = conf_int,
                aes(y = theorical_expression, 
                    x = ploidy,
                    linetype = theoric
                ), 
                inherit.aes = F, 
                # linetype = 'dashed',
                alpha = .6, 
                color = 'deepskyblue4') + 
      theme_bw() + 
      # geom_hline(data = conf_int, 
      #            aes(yintercept = CI_lower))
      labs(y = 'Expression', 
           x = 'Ploidy') + 
      scale_shape_manual(values = c('Wild-type' = 16, 'Mutated' = 17)) +
      scale_linetype_manual(values = c('Theoretical expression' = 'dashed')) +
      scale_fill_manual(values = ci_col) + 
      ggtitle(title_p) + 
      coord_cartesian(ylim = c(min_expr, #min(df$estimate_expression), 
                               max_expr)) + 
      guides(
        shape = guide_legend(title = 'Mutation status', 
                             nrow = 2, 
                             order = 2,
                             title.position = 'top', 
                             override.aes = list(fill = "white")), 
        linetype = guide_legend(title = '', 
                                order = 3, 
                                override.aes = list(fill = "white")), 
        fill = guide_legend(title = 'Confidence interval', 
                            ncol = 2)) + 
      theme(legend.position = 'right') + 
      facet_wrap(~Assay, nrow = 2)
  } else {
    p = ggplot() +
      ggtitle(title_p)
  }
  return(p)
}

# plot proportion of samples with a specific class
plot_classes_by_gene = function(x, gene, classes_cols) {
  x = x %>%
    filter(hgnc_symbol == gene) 
  if(nrow(x) > 0) {
    x %>% 
      pivot_longer(cols = c(prop_prot, prop_rna), names_to = 'assay', values_to = 'prop') %>% 
      mutate(assay = gsub('prop_', '', assay)) %>% 
      mutate(assay = factor(assay, levels = c('rna', 'prot'))) %>% 
      ggplot(aes(x = class, y = prop, fill = class)) + 
      geom_bar(stat = 'identity') + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            legend.position = 'bottom', 
            legend.title = element_text(hjust = 0.5, 
                                        vjust = 0)
      ) + 
      xlab('Dosage class') + 
      ylab('Proportion of genes') +
      scale_fill_manual(values = classes_cols) + 
      # scale_fill_manual(values = c('rna' = 'steelblue', 'prot' = 'goldenrod'), 
      #                   labels = c("rna" = "RNA", "prot" = "Protein")) +
      guides(fill = guide_legend(
        title = 'Dosage classes', 
        ncol = 3, 
        title.position = 'top'
      )) + 
      facet_wrap(~assay, labeller = labeller(assay = c('rna' = 'RNA', 'prot' = 'Protein'))) + 
      ggtitle(gene)
  } else {
      ggplot()
    }
}


### by cohort ----

# prevalent class distribution
plot_prevalent_class = function(x, classes_cols, omics_cols) {
  x %>% 
    group_by(hgnc_symbol, Assay) %>% 
    dplyr::slice_max(prop) %>% 
    dplyr::group_by(class, Assay) %>% 
    count() %>% 
    ggplot(aes(x = reorder(class, +n), y = n, fill = class)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(values = classes_cols) + 
    theme_bw() + 
    light_theme() + 
    labs(x = 'Dosage classes', y = 'Number of genes', 
         title = 'Prevalent class') +
    coord_flip() + 
    facet_wrap(~Assay) + 
    # facet_wrap2(~Assay, 
    #             strip = strip_themed(
    #               background_x = elem_list_rect(fill = omics_cols, alpha = 0.1)
    #             )) +
    guides(fill = guide_legend(title = 'Dosage classes', 
                               ncol = 3, 
                               title.position = 'top'))
}

# samples barplot 
plot_cls_by_sample = function(x, dosage_colors) {
  x %>% 
    filter(!is.na(sample)) %>% 
    pivot_longer(cols = c(prop_prot, prop_rna), names_to = 'Assay', values_to = 'prop') %>% 
    mutate(Assay = gsub('prop_', '', Assay)) %>% 
    mutate(Assay = factor(Assay, levels = c('rna', 'prot'))) %>% 
    ggplot(aes(
      y = sample, x = prop, fill = class
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
                                ncol = 3, order = 1, 
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


levels_cls_numeric = factor(levels_cls, levels = levels_cls) %>% as.numeric()

# ht_colors = setNames(
#   nm = levels_cls_numeric,
#   object = c(
#     '#CBCBCB',
#     '#BF092F', 
#     '#1E5128', 
#     '#3B9797', 
#     '#28518A', 
#     '#5D2F77', 
#     '#C95792', 
#     # '#809D3C',
#     '#5D8736')
#   # object = c('#BF092F', '#001BB7', '#31694E', '#3B9797', '#6B3F69',  '#E195AB', '#84994F', '#CBCBCB')
# )

ht_colors = setNames(
  nm = levels_cls_numeric, 
  object = dosage_colors %>% unname()
)
annotation_colors = list('Assay' = setNames(object = c('steelblue', 'goldenrod'), c('RNA', 'Protein')))

plot_heatmap = function(rna, 
                        prot, 
                        genes, 
                        cols = ht_colors,
                        annotation_colors, 
                        levels_cls, 
                        title) {
  
  class_merge = merging_omics(prot_cls = prot, rna_cls = rna)
  x = class_merge %>%
    dplyr::select(hgnc_symbol, sample, class_RNA, class_Prot)
  
  # levels_cls = c('Not_assigned', 'Dosage_compensating', 'Dosage_sensitive', 'Dosage_insensitive', 
  #                'Haploinsufficient', 'Mutation_sensitive', 'Mutation_insensitive', 
  #                'Control')
  
  # wide the data independently and match them after
  
  x_wide = x %>% 
    mutate(class_RNA = ifelse(is.na(class_RNA), 'Not_assigned', class_RNA)) %>%
    mutate(class_Prot = ifelse(is.na(class_Prot), 'Not_assigned', class_Prot)) %>%
    mutate(class_RNA = factor(class_RNA, levels = levels_cls))  %>%
    mutate(class_Prot = factor(class_Prot, levels = levels_cls))  %>%
    mutate(Protein = as.numeric(class_Prot),
           RNA = as.numeric(class_RNA)) %>%
    dplyr::select(-c(starts_with('class_'))) %>% 
    dplyr::filter(hgnc_symbol %in% genes) %>% 
    pivot_wider(names_from = sample, values_from = c(Protein, RNA), values_fill = 1) %>% 
    tibble::column_to_rownames('hgnc_symbol')
  
  # create annotation 
  ann_data = tibble(
    sample = colnames(x_wide)
  ) %>% 
    mutate(Assay = str_extract(sample, 'Protein|RNA')) %>% 
    select(Assay)
  
  ann_ht = create_annotation(x = ann_data, ann_colors = annotation_colors, position = 'column')
  ht = Heatmap(x_wide, 
          col = cols, 
          name = 'Dosage classes', 
          bottom_annotation = ann_ht,
          column_title = title,
          heatmap_legend_param =
            list(at = names(cols), labels = levels_cls, direction = "horizontal", nrow = 3))
          
  draw(ht, heatmap_legend_side = "bottom", 
       annotation_legend_side = "bottom")
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

# further analysis -----

## existing relationships with the mutation type and the classification 

check_corr_mut_type = function(x) {
  
  x %>% 
    pull(cls_dosage) %>% unique
  
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


# K means clustering

# kmeans_estimates = function(x, seed = 12345, pth = 0.05) {
#   # cluster univariate the genes per each coefficient
#   x = x %>%
#     select(hgnc_symbol, observed_expr, sample) %>%
#     group_by(hgnc_symbol) %>% 
#     # tibble::column_to_rownames('hgnc_symbol') %>% 
#     group_split()
#   names(x) = lapply(x, function(s) s$hgnc_symbol %>% unique)
#   
#   # tidyr::pivot_wider(names_from = term, values_from = estimate, values_fill = 0) %>%
#   
#   
#   
#   k_clus = lapply(x, function(s) {
#     
#     set.seed(seed)
#     
#     tryCatch({
#       
#       data = s %>% 
#         tibble::column_to_rownames('sample') %>% 
#         select(observed_expr)
#       sil_plot = fviz_nbclust(data, FUNcluster = kmeans, method = "silhouette")
#       k = sil_plot@layers$geom_vline$data$xintercept
#       
#       clustering = kmeans(data, k, nstart = 10)
#       
#       clusters = tibble(
#         sample = names(clustering$cluster), 
#         cluster = unname(clustering$cluster)
#       )
#       
#       list(kmeans = clustering, sil_plot = sil_plot, clusters = clusters, data = s)
#     }, 
#     error = function(e) {NA})
#   })
#   
#   return(k_clus)
#   
# }


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
