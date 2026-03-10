# compensation score plots 

# colors 
category_colors <- c(
  "(RNA-prot heavy)" = "#AD002AB2",
  "(Prot-heavy)" = "#E18727B2",
  "(RNA-heavy)" = "#20854Eb2",
  "(RNA-Prot light)" = "#00468BB2",
  "Intermediate/Other" = "gainsboro"
)

# proportion of genes

plot_props_by_karyo = function(x, colors, gene_list) {
  
  x %>% 
    filter(name %in% gene_list) %>% 
    ungroup() %>% 
    group_by(karyotype) %>% 
    mutate(tot = n()) %>% 
    group_by(karyotype, reg_group) %>% 
    mutate(prop = n()/tot) %>% 
    dplyr::select(karyotype, reg_group, prop) %>% 
    distinct() %>% 
    ggplot(aes(x = karyotype, 
               y = prop,
               fill = reg_group)) + 
    geom_bar(stat = 'identity', alpha = 1) + 
    theme_bw() + 
    scale_fill_manual(values = colors)
    
}


# plot props

plot_props = function(x, colors, gene_list) {
  
  x %>% 
    filter(name %in% gene_list) %>% 
    ungroup() %>% 
    ggplot(aes(x = reg_group, 
               fill = reg_group)) + 
    geom_bar(stat = 'count', alpha = 1) + 
    theme_bw() + 
    scale_fill_manual(values = colors)
  
}


