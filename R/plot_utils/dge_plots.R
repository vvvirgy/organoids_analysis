# plot utils 
expr_cols = setNames(
  nm = c('up', 'down', 'not differential'), 
  object = c('firebrick', 'darkblue', 'grey30')
)

classify_genes = function(x, pth = .05, fth = .75) {
  
  x %>% 
    ungroup %>% 
    mutate(significance = ifelse(adj_pval <= pth, 'significant', 'ns')) %>% 
    mutate(fc_cls = case_when(
      (significance == 'significant' & lfc >= -fth & lfc <= fth) ~ 'not differential',
      (lfc <= -fth & significance == 'significant') ~ 'down', 
      (lfc >= fth & significance == 'significant') ~ 'up', 
      .default = 'ns'
    )) %>% 
    mutate(fc_cls = factor(fc_cls, levels = c('up', 'down', 'not differential', 'ns')))
  
}

# volcano plot 
plot_volcano = function(x, omic_list, cols, pth = .05, fth = .75) {
  x %>% 
    filter(omic %in% omic_list) %>% 
    # mutate(sign = ifelse(adj_pval <= pth, TRUE, FALSE)) %>% 
    # mutate(cls = case_when(
    #   (sign == TRUE & lfc >= fth ) ~ 'up', 
    #   (sign == TRUE & lfc <= fth ) ~ 'down', 
    #   .default = 'ns'
    # )) %>% 
    filter(!is.na(lfc)) %>% 
    ggplot(aes(x = lfc, 
               y = -log10(adj_pval), 
               color = fc_cls, 
               alpha = significance)) + 
    geom_point() + 
    theme_bw() + 
    facet_grid(omic~karyotype, scales = 'free_y') + 
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = setNames(c(1, .3), c('significant', 'ns'))) + 
    guides(alpha = guide_legend(theme = theme(legend.position = 'None'))) + 
    geom_hline(yintercept = -log10(pth), linetype = 'dashed', alpha = .7, colour = '#151515') + 
    geom_vline(xintercept = c(-fth, fth), linetype = 'dashed', alpha = .7, colour = '#151515')
}

multi_omics_fc_plt = function(x, gene, colors, fth = .75) {
  x %>% 
    filter(name == gene) %>% 
    ggplot(aes(y = lfc, 
               x = karyotype, 
               fill = omic, 
               alpha = significance, 
               color = significance)) + 
    geom_bar(stat = 'identity', 
             position = 'dodge') +
    scale_alpha_manual(values = setNames(nm = c('significant', 'ns'), c(1, .5))) +
    scale_color_manual(values = setNames(nm = c('significant', 'ns'), c('black', 'transparent'))) + 
    scale_fill_manual(values = colors) + 
    theme_bw() + 
    geom_hline(yintercept = c(-fth, fth), linetype = 'dashed', alpha = .7, colour = '#151515') +
    ggtitle(gene)
}


# plot multivariate among the two omics
bg_cols = c('Same sign' = '#A5C89E', 'Opposite sign' = '#E5BA41')

plot_omics_comparison = function(x, 
                                 filter = T, 
                                 bg_colors
) {
  
  if(filter == TRUE) {
    x = x %>% 
      filter(sign_RNA =='significant', sign_prot=='significant')
    
    title = 'genes significant in both tests'
  } else {
    title = 'all genes'
  }
  
  x %>%
    ggplot(aes(
      x = lfc_RNA, 
      y = lfc_protein
    )) + 
    geom_rect(data = 
                data.frame(xmin = 0, 
                           xmax = Inf, 
                           ymin = 0, 
                           ymax = Inf, 
                           cls = 'Same sign'), 
              aes(xmin=xmin,
                  xmax=xmax,
                  ymin=ymin,
                  ymax=ymax, 
                  fill = cls), 
              inherit.aes = F, 
              alpha = 0.3) + 
    geom_rect(data = 
                data.frame(xmin = 0, 
                           xmax = -Inf, 
                           ymin = 0, 
                           ymax = -Inf, 
                           cls = 'Same sign'), 
              aes(xmin=xmin,
                  xmax=xmax,
                  ymin=ymin,
                  ymax=ymax, 
                  fill = cls), 
              inherit.aes = F, 
              alpha = 0.3) + 
    geom_rect(data = 
                data.frame(xmin =  0,
                           xmax =  -Inf,
                           ymin =  0,
                           ymax =  Inf,
                           cls = 'Opposite sign'), 
              aes(xmin=xmin,
                  xmax=xmax,
                  ymin=ymin,
                  ymax=ymax, 
                  fill = cls), 
              inherit.aes = F, 
              alpha = 0.3)  + 
    geom_rect(data = 
                data.frame(xmin =  0,
                           xmax =  Inf,
                           ymin =  0,
                           ymax =  -Inf,
                           cls = 'Opposite sign'), 
              aes(xmin=xmin,
                  xmax=xmax,
                  ymin=ymin,
                  ymax=ymax, 
                  fill = cls), 
              inherit.aes = F, 
              alpha = 0.3)  + 
    geom_point() + 
    facet_wrap(~karyotype) + 
    labs(y = 'Protein FC', 
         x = 'RNA FC') + 
    theme_light() + 
    scale_fill_manual(values = bg_colors) + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    guides(fill = guide_legend(title = 'FC signs')) + 
    theme(legend.position = 'bottom') + 
    ggtitle(title)
  
}


# plot densities
plot_densities = function(x, filter = T) {

  if(filter == TRUE) {
    x = x %>% 
      filter(significance == 'significant') 
  } 
  
  x %>% 
    filter(!is.na(lfc)) %>% 
    filter(!is.na(fc_cls)) %>% 
    ggplot(aes(x = karyotype, fill = karyotype, y = lfc)) + 
    geom_violin() + 
    theme_bw() + 
    facet_wrap(~omic) + 
    scale_fill_brewer(palette = 'Pastel2') + 
    ggpubr::stat_compare_means(comparisons = list(
      c('1:0', '2:0'), 
      c('1:0', '2:1'),
      c('1:0', '2:2'),
      c('2:0', '2:1'), 
      c('2:0', '2:2'), 
      c('2:2', '2:1')
    ))

}

# barplots to compare the number of up and down genes
plot_barplots = function(x, filter_significance = TRUE, cols) {

  if(filter_significance == TRUE) {
    x = x %>%   
      filter(significance == 'significant') 
  }
  
  x %>% 
    filter(!is.na(lfc)) %>% 
    filter(!is.na(fc_cls)) %>% 
    ggplot(aes(karyotype, fill = fc_cls)) + 
    geom_bar(stat = 'count', position = 'dodge') + 
    theme_light() + 
    facet_wrap(~omic) + 
    scale_fill_manual(values = cols) + 
    labs(x = 'Karyotype', 
         y = 'Number of genes') + 
    guides(fill = guide_legend(title = 'FC class'))
  
}

# plotting sankeys among the results
sankey_cols = setNames(
  nm = c('not estimated', 'ns', 'significant, up', 'significant, down', "significant, not differential"), 
  object = c('#FFF0CE', '#BFC9D1', 'firebrick', 'darkblue', '#5D0E41')
)



plot_sankey = function(x, 
                       strat = 'karyotype', 
                       karyos = c('1:0', '2:0', '2:1', '2:2'), 
                       omics = c('RNA', 'protein'), 
                       facet = 'omic', 
                       cols) {
  
  x = x %>% 
    mutate(cls = paste(significance, fc_cls, sep = ', ')) %>% 
    mutate(cls = ifelse(
      significance == 'ns', 'ns', cls
    )) %>% 
    mutate(cls = ifelse(
      (is.na(significance) & is.na(fc_cls)), 'not estimated', cls)
    ) %>% 
    mutate(stratum = factor(.data[[strat]])) %>% 
    group_by(karyotype, omic, cls) %>% 
    mutate(freq = row_number()/n())
    
  
  x %>% 
    # filter(omic == 'RNA') %>% 
    filter(karyotype %in% karyos) %>%
    filter(omic %in% omics) %>% 
    ggplot(
      aes(x = stratum, 
          stratum = cls, 
          alluvium = name,
          y = freq,
          fill = cls) 
    ) + 
    scale_fill_manual(values = cols)+
    # geom_flow(stat = "alluvium", lode.guidance = "frontback") +
    geom_flow(alpha = .7) +
    geom_stratum() +
    theme_light() + 
    facet_wrap(as.formula(paste0('~', facet)), scales = 'free')
    # facet_wrap(~ {{ facet }}, scales = 'free')
}


# plot the results of the enrichment

plot_enrichment_results = function(df, highlight, colors, order_by = NULL, sources) {
  
  if(highlight == TRUE) {
    df = df %>% 
      filter(highlighted == TRUE)
  }
  
  df %>% 
    mutate(ratio = intersection_size/term_size) %>% 
    ggplot(aes(
      x = ratio, 
      y = term_name, 
      color = p_value, 
      size = intersection_size
    )) + 
    geom_point() + 
    theme_light() + 
    scale_color_continuous(palette = 'Blues')
    
    
  
}

enrichment_heatmap = function(x, source_list, pth = .05, genes_number = 4, highlight = TRUE) {
  
  df = x$result
  
  if(highlight == TRUE) {
    df = df %>% 
      filter(highlighted == TRUE)
  } else {
    
    df = df %>% 
      filter(source %in% source_list) 
  }
  
  df = df %>% 
    separate_rows(intersection, sep = ',') %>% 
    select(-parents) %>% 
    mutate(adj_pval = p.adjust(p_value, method = 'BH')) %>% 
    filter(adj_pval <= pth) %>% 
    filter(intersection_size >= genes_number)
  
  df %>% 
    ggplot(aes(
      y = reorder(term_name, +intersection_size), 
      x = intersection, 
      fill = -log10(adj_pval)
    )) + 
    geom_tile() + 
    theme_light() + 
    labs(y = 'Biological term', 
         x = 'Gene') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
    facet_wrap(~source, scales = 'free', ncol = 3)

}

# create a heatmap using complexheatmap
# 
# prepare_heatmap_data = function(res, 
#                                 deg, 
#                                 source_list, 
#                                 pth = .05, 
#                                 genes_number = 4, 
#                                 highlight = TRUE, 
#                                 kayo, 
#                                 omics = c('RNA', 'protein'), 
#                                 direction =) {
#   
#   genes = res$meta$query_metadata$queries$query_1 %>% unique
#   # filter the deg table 
#   
#   deg = deg %>% 
#     filter(karyotype == karyo) %>% 
#     filter(omic %in% omics) %>% 
#     filter(name %in% genes)
#     
#   
#   
#   
#   
# }


enrichment_heatmap = function(x, source_list, pth = .05, genes_number = 4, highlight = TRUE) {
  
  df = x$result
  
  if(highlight == TRUE) {
    df = df %>% 
      filter(highlighted == TRUE)
  } else {
    
    df = df %>% 
      filter(source %in% source_list) 
  }
  
  df = df %>% 
    separate_rows(intersection, sep = ',') %>% 
    select(-parents) %>% 
    mutate(adj_pval = p.adjust(p_value, method = 'BH')) %>% 
    filter(adj_pval <= pth) %>% 
    filter(intersection_size >= genes_number)
  
  colnames(df)
  
  
  df %>% 
    ggplot(aes(
      y = reorder(term_name, +intersection_size), 
      x = intersection, 
      fill = -log10(adj_pval)
    )) + 
    geom_tile() + 
    theme_light() + 
    labs(y = 'Biological term', 
         x = 'Gene') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
    facet_wrap(~source, scales = 'free', ncol = 3)
  
}


  
  
  
  
  
  
  



   