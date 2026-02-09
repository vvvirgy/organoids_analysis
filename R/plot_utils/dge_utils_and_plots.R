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
    guides(alpha = guide_legend(theme = theme(legend.position = 'None'), 
                                title = 'Significance class (padj <= 0.05)'), 
           color = guide_legend(title = 'FC class')) + 
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
    filter(!is.na(lfc), !is.na(adj_pval)) %>% 
    mutate(cls = paste(significance, fc_cls, sep = ', ')) %>% 
    mutate(cls = ifelse(
      significance == 'ns', 'ns', cls
    )) %>% 
    mutate(cls = ifelse(
      (is.na(significance) & is.na(lfc)), 'not estimated', cls)
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
fc_cols = setNames(
  nm = c('up', 'down', 'not differential', 'not present'), 
  object = c('firebrick', 'darkblue', '#5D0E41', '#FFF0CE')
)

sankey_annotations_by_omic = function(x, 
                                      filter_terms = TRUE,
                                      karyo, 
                                      omic_list = c('RNA', 'protein'),
                                      pth = .05, 
                                      source_list = NULL, 
                                      cols
                                      
) {
  
  # rna = paste('RNA', karyo, c('up', 'down', 'not differential'), sep = '_')
  
  classes = cross(list(omic_list, karyo, c('up', 'down', 'not differential'))) %>% map_chr(paste, sep = "_", collapse = "_")
  
  x = x[c(classes)]
  
  df = lapply(x %>% names, function(data) {
    x[[data]]$result %>% 
      mutate(category = data)
  }) %>% bind_rows() 
  
  if(filter_terms == TRUE) {
    df = df %>% 
      filter(p_value <= pth)
  }
  
  if(!is.null(source_list)) {
    df = df %>% 
      filter(source %in% source_list)
  }
  
  # merge the results in a better way
  
  df = df %>% 
    separate(category, sep = '_', into = c('omic', 'karyotype', 'fc_cls')) %>% 
    dplyr::select(term_name, fc_cls, omic, karyotype) %>% 
    distinct() %>% 
    group_by(term_name, karyotype) %>% 
    filter(n() == 1) %>% 
    ungroup() %>% 
    complete(
      term_name,
      karyotype,
      omic,
      fill = list(
        fc_cls = 'not present'
      )
    ) %>% 
    # mutate(omic = factor(omic, levels = c('RNA', 'protein'))) %>%
    mutate(karyotype = factor(karyotype)) %>% 
    mutate(fc_cls = factor(fc_cls, levels = c('not present', 'down', 'not differential', 'up') %>% rev)) %>% 
    group_by(karyotype) %>% 
    mutate(freq = row_number()/n())
    
  df %>% 
    # filter(omic == 'RNA') %>% 
    ggplot(
      aes(x = karyotype, 
          stratum = fc_cls, 
          alluvium = term_name,
          y = freq,
          fill = fc_cls
          ) 
    ) + 
    scale_fill_manual(values = cols)+
    # geom_flow(stat = "alluvium", lode.guidance = "frontback") +
    geom_flow(alpha = .7) +
    geom_stratum() +
    theme_light() + 
    # ggtitle(paste('karyotype', karyo)) + 
    labs(y = '',
         x = 'Omic Assay') + 
    guides(fill = guide_legend(title = 'FC class', position = 'bottom')) + 
    facet_wrap(~omic, scales = 'free') + 
    ggtitle('Annotations among karyotypes')
  
}


sankey_annotations = function(x, 
                              # by = 'omic', 
                              filter_terms = TRUE,
                              karyo, 
                              pth = .05, 
                              source_list = NULL, 
                              cols
                              
) {
  
  rna = paste('RNA', karyo, c('up', 'down', 'not differential'), sep = '_')
  protein = paste('protein', karyo, c('up', 'down', 'not differential'), sep = '_')
  
  x = x[c(rna, protein)]
  
  df = lapply(x %>% names, function(data) {
    x[[data]]$result %>% 
      mutate(category = data)
  }) %>% bind_rows() 
  
  if(filter_terms == TRUE) {
    df = df %>% 
      filter(p_value <= pth)
  }
  
  if(!is.null(source_list)) {
    df = df %>% 
      filter(source %in% source_list)
  }
  
  # merge the results in a better way
  
  df = df %>% 
    separate(category, sep = '_', into = c('omic', 'karyotype', 'fc_cls')) %>% 
    select(term_name, fc_cls, omic) %>% 
    distinct() %>% 
    group_by(term_name, omic) %>% 
    filter(n() == 1) %>% 
    ungroup() %>% 
    complete(
      term_name,
      omic,
      fill = list(
        fc_cls = 'not present'
      )
    ) %>% 
    mutate(omic = factor(omic, levels = c('RNA', 'protein'))) %>%
    mutate(fc_cls = factor(fc_cls, levels = c('not present', 'down', 'not differential', 'up') %>% rev)) %>% 
    group_by(omic) %>% 
    mutate(freq = row_number()/n())
  
  df %>% 
    # filter(omic == 'RNA') %>% 
    ggplot(
      aes(x = omic, 
          stratum = fc_cls, 
          alluvium = term_name,
          y = freq,
          fill = fc_cls
      ) 
    ) + 
    scale_fill_manual(values = cols)+
    # geom_flow(stat = "alluvium", lode.guidance = "frontback") +
    geom_flow(alpha = .7) +
    geom_stratum() +
    theme_light() + 
    ggtitle(paste('karyotype', karyo)) + 
    labs(y = '',
         x = 'Omic Assay') + 
    guides(fill = guide_legend(title = 'FC class', position = 'bottom'))
  
}


plot_enrichment_results = function(df, highlight, colors, order_by = NULL, sources) {
  
  if(highlight == TRUE) {
    df = df %>% 
      filter(highlighted == TRUE)
  } else if(!is.null(sources)) {
    df = df %>% 
      filter(source %in% sources)
  }
  
  df %>% 
    # mutate(ratio = intersection_size/term_size) %>% 
    ggplot(aes(
      x = precision, 
      y = reorder(term_name, +precision), 
      color = p_value, 
      size = intersection_size
    )) + 
    geom_point() + 
    theme_light() + 
    scale_fill_distiller(palette = 'Blues') + 
    labs(
      y = 'Biological terms', 
      x = 'Precision'
    )

    # scale_color_brewer(palette = 'Blues')

}

# enrichment_barplot

enrichment_barplot = function(x, 
                              pth = .05, 
                              source_list = NULL, 
                              cols) {
  
  df = lapply(names(x), function(data){
    x[[data]]$result %>% 
      mutate(category = data)
  }) %>% bind_rows() %>% 
    tidyr::separate(category, into = c('omic' ,'karyotype', 'cls'), sep = '_', convert = T)
  
  if(!is.null(source_list)) {
    df = df %>% 
      filter(source %in% source_list)
  }
  
  df = df %>% 
    filter(p_value <= pth)
  
  df %>% 
    ggplot(aes(
      cls, 
      fill = cls
    )) + 
    geom_bar(stat = 'count') + 
    theme_light() + 
    ggh4x::facet_nested_wrap(vars(omic, karyotype), ncol = 1, strip.position = 'right') +   
    coord_flip() + 
    scale_fill_manual(values = cols) +
    labs(
      y = 'Number of terms', 
      x = 'FC class'
    )
  
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
    # mutate(adj_pval = p.adjust(p_value, method = 'BH')) %>% 
    filter(p_value <= pth) %>% 
    filter(intersection_size >= genes_number)
  
  df %>% 
    ggplot(aes(
      y = reorder(term_name, +intersection_size), 
      x = intersection, 
      fill = -log10(p_value)
    )) + 
    geom_tile() + 
    theme_light() + 
    labs(y = 'Biological term', 
         x = 'Gene') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
    facet_wrap(~source, scales = 'free', ncol = 3)

}


plot_terms_fc = function(res, 
                         deg, 
                         source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP'), 
                         pth = .05, 
                         highlight = FALSE,
                         genes_number = 0, 
                         karyo, 
                         omics) {
  
  genes = res$meta$query_metadata$queries$query_1 %>% unique
  df = res$result %>% 
    separate_rows(intersection, sep = ',') %>% 
    filter(intersection_size >= genes_number) 
  
  if(highlight == TRUE) {
    df = df %>% 
      filter(highlighted == TRUE)
  } else {
    df = df %>% 
      filter(source %in% source_list)
  }
  
  # filter the deg table
  
  deg = deg %>%
    filter(karyotype == karyo) %>%
    filter(omic %in% omics) %>%
    filter(name %in% genes) %>% 
    mutate(omic = as.character(omic))
  
  df = df %>% 
    left_join(., deg, by = join_by('intersection' == 'name')) %>%
    dplyr::select(term_name, intersection, lfc, omic, karyotype, source, intersection_size) 
  
  # ann_data = df %>% 
  #   select(source, term_name, intersection_size) %>% 
  #   distinct() %>% 
  #   select(source, intersection_size)
  # 
  # ann = create_annotation(ann_data %>% select(source), 
  #                         ann_colors = ann_colors, 
  #                         position = 'row')
  # bar_ann = rowAnnotation('number of genes' = anno_barplot(ann_data$intersection_size, baseline = 0,  
  #                                                          axis_param = list(direction = "reverse")))
  
  df = df %>% 
    mutate(omic = as.character(omic)) %>% 
    group_by(term_name, omic, source) %>% 
    summarise(mean_fc = mean(lfc,na.rm = T), 
              sd_fc = sd(lfc, na.rm = T)) %>% 
    mutate(term_name = 
             ifelse(duplicated(.data[['term_name']]), paste0(term_name, ' (', source, ')'), term_name)) 
  
  
  df %>% 
    ggplot(aes(y = reorder(term_name, +mean_fc), 
               x = mean_fc, 
               color = mean_fc)) + 
    geom_point() + 
    theme_light() +
    geom_errorbar(aes(xmin = mean_fc-sd_fc, xmax = mean_fc+sd_fc)) + 
    labs(
      x = 'Mean FC per term', 
      y = 'Biological term'
    ) + 
    facet_wrap(~source, scales = 'free_y', ncol = 1)
  
}

# create a heatmap using complexheatmap


# preprocess_data_for_heatmap = function(res,
#                                        deg,
#                                        source_list,
#                                        pth = .05,
#                                        genes_number = 4,
#                                        highlight = FALSE,
#                                        karyo,
#                                        omics = c('RNA', 'protein'), 
#                                        direction) {
#   
#   if(names(res) )
#   
# }

create_annotation = function(x, ann_colors, position) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position,
                                    annotation_legend_param = list(nrow = 3, width = 12, by_row = T))
  
}

source_colors = list(
  source = setNames(nm = c('GO:BP', 'GO:MF', 'GO:CC', 'REAC', 'KEGG', 'WP', 'TF', 'MIRNA', 'CORUM'), 
             object = rcartocolor::carto_pal(n = 10, name = 'Prism')[-10]))


plot_fc_heatmap = function(res,
                           deg,
                           source_list,
                           pth = .05,
                           genes_number = 4,
                           highlight = FALSE,
                           karyo,
                           omics = c('RNA', 'protein'), 
                           direction, 
                           ann_colors) {
  
  genes = res$meta$query_metadata$queries$query_1 %>% unique
  df = res$result %>% 
    separate_rows(intersection, sep = ',') %>% 
    filter(intersection_size >= genes_number) 
  
  if(highlight == TRUE) {
    df = df %>% 
      filter(highlighted == TRUE)
  } else {
    df = df %>% 
      filter(source %in% source_list)
  }
  
  # filter the deg table
  
  deg = deg %>%
    filter(karyotype == karyo) %>%
    filter(omic %in% omics) %>%
    filter(name %in% genes) %>% 
    mutate(omic = as.character(omic))

  df = df %>% 
    left_join(., deg, by = join_by('intersection' == 'name')) %>%
    select(term_name, intersection, lfc, omic, karyotype, source, intersection_size) 
  
  ann_data = df %>% 
    select(source, term_name, intersection_size) %>% 
    distinct() %>% 
    select(source, intersection_size)
  
  ann = create_annotation(ann_data %>% select(source), 
                          ann_colors = ann_colors, 
                          position = 'row')
  bar_ann = rowAnnotation('number of genes' = anno_barplot(ann_data$intersection_size, baseline = 0,  
                              axis_param = list(direction = "reverse")))
  
  df = df %>% 
    mutate(omic = as.character(omic)) %>% 
    group_by(term_name, omic, source) %>% 
    summarise(mean_fc = mean(lfc,na.rm = T)) %>% 
    mutate(term_name = 
             ifelse(duplicated(.data[['term_name']]), paste0(term_name, ' (', source, ')'), term_name)) %>% 
    select(-source) %>%
    pivot_wider(names_from = omic, values_from = mean_fc)
  
  df = df %>% 
    tibble::column_to_rownames('term_name') %>% 
    as.matrix()
  
  if(direction == 'up') {
    col_scale = circlize::colorRamp2(breaks = c(0, mean(df), max(df)), 
                                     colors = c('white', '#CF0F0F', '#740A03'))
    ht = ComplexHeatmap::Heatmap(df, 
                            col = col_scale,
                            right_annotation = ann, 
                            left_annotation = bar_ann, 
                            name = 'mean(FC)'
                            )
    draw(ht, 
         merge_legend = TRUE, 
         heatmap_legend_side = "bottom", 
         annotation_legend_side = "bottom")
    
  } else if(direction == 'down') {
    col_scale = circlize::colorRamp2(breaks = c(min(df), mean(df), 0) %>% rev, 
                                     colors = c('white', '#9ECAE1', '#3182BD'))
    ht = ComplexHeatmap::Heatmap(df, 
                            col = col_scale, 
                            right_annotation = ann, 
                            left_annotation = bar_ann, 
                            name = 'mean(FC)')
    draw(ht, 
         merge_legend = TRUE, 
         heatmap_legend_side = "bottom", 
         annotation_legend_side = "bottom")
    
  } else {
    ht = ComplexHeatmap::Heatmap(df, 
                                 right_annotation = ann, 
                                 left_annotation = bar_ann, 
                                 name = 'mean(FC)')
    
    draw(ht, 
         merge_legend = TRUE, 
         heatmap_legend_side = "bottom", 
         annotation_legend_side = "bottom")
    
  }
  
}



# RUNNING OF GSEA ON DIFFERENT DATABASES ------

run_gsea = function(genes, 
                    databases = c('GO', 'KEGG', 'WP', 'REAC')) {
  
  res = list()
  
  if('GO' %in% databases) {
    
    res$GO <- gseGO(
      geneList = genes,
      OrgDb = org.Hs.eg.db, 
      pvalueCutoff = 0.05, 
      ont = "ALL",
      minGSSize = 10,
      maxGSSize = 250,
      keyType = "ENTREZID",
      pAdjustMethod = "BH",
      by = "fgsea"
    )
  }
  
  if('REAC' %in% databases) {
    res$REAC <- gsePathway(
      geneList = genes,
      organism = "human",
      pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, 
      minGSSize = 10,
      maxGSSize = 250,
      by = "fgsea"
    )
    
  }
  
  if('WP' %in% databases) {
    res$WP = gseWP(genes, 
                   organism = "Homo sapiens", 
                   pvalueCutoff = 0.05, 
                   minGSSize = 10,
                   maxGSSize = 250,
                   pAdjustMethod = "BH")
  }
  
  if('KEGG' %in% databases) {
    res$KEGG <- gseKEGG(
      geneList = genes,
      pvalueCutoff = 0.05,
      minGSSize = 10,
      maxGSSize = 250,
      organism = "hsa",
      pAdjustMethod = "BH"
    )
  }
  
  return(res)
}

plot_nes = function(x) {
  x %>% 
    ggplot(aes(y = Description, 
               x = NES, 
               color = p.adjust)) + 
    geom_point() + 
    scale_color_viridis_c(option = 'viridis') + 
    theme_light()
}
  
# enrichment_heatmap = function(x, source_list, pth = .05, genes_number = 4, highlight = F) {
#   
#   df = x$result
#   
#   if(highlight == TRUE) {
#     df = df %>% 
#       filter(highlighted == TRUE)
#   } else {
#     
#     df = df %>% 
#       filter(source %in% source_list) 
#   }
#   
#   df = df %>% 
#     separate_rows(intersection, sep = ',') %>% 
#     select(-parents) %>% 
#     mutate(adj_pval = p.adjust(p_value, method = 'BH')) %>% 
#     filter(adj_pval <= pth) %>% 
#     filter(intersection_size >= genes_number)
#   
#   colnames(df)
#   
#   
#   df %>% 
#     ggplot(aes(
#       y = reorder(term_name, +intersection_size), 
#       x = intersection, 
#       fill = -log10(adj_pval)
#     )) + 
#     geom_tile() + 
#     theme_light() + 
#     labs(y = 'Biological term', 
#          x = 'Gene') + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
#     facet_wrap(~source, scales = 'free', ncol = 3)
#   
# }


  
  
  
  
  
  
  



   