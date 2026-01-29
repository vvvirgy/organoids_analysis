#  create the annotations ------

create_design = function(x, cna_state, gene, cond = 'karyotype') {
  
  cna_state = cna_state %>% 
    filter(hgnc_symbol == gene) %>% 
    group_by(sample, hgnc_symbol) %>% 
    filter(n()==1)
  
  ann = x %>% 
    full_join(., cna_state, by = join_by('fixed_name' == 'sample')) %>% 
    filter(!is.na(.data[[cond]])) %>% 
    filter(!is.na(PDO))
  
  ann_new = ann %>% 
    group_by_at(cond) %>%
    rename(replicate_name = replicate) %>% 
    filter(!is.na(replicate_name)) %>%
    mutate(replicate = seq_len(n())) %>% 
    dplyr::rename(label = replicate_name) %>%
    rename(condition := !!sym(cond)) 
  
  return(ann_new)
  
}

dep_preprocesing = function(data, design) {
  
  data = data %>% 
    dplyr::select(PG.ProteinGroups, PG.Genes, any_of(design$label))
  # create the dep object
  data = make_unique(data, names = 'PG.Genes', ids = 'PG.ProteinGroups', delim = '_')
  samples_index = which(!colnames(data) %in% c('PG.ProteinGroups', 'PG.Genes', 'name', 'ID'))
  
  # from now values are log transformed
  dep = DEP::make_se(proteins_unique = data, columns = samples_index, expdesign = design)
  
  # filtering out proteins with missing values in all replicates
  dep_filt = filter_missval(dep, thr = 0)
  
  dep_norm = normalize_vsn(dep_filt)
  
  return(dep_norm)
}

differential_expression = function(dep_norm, ctrl) {
  # trying to test for dge over kras 
  data_diff <- test_diff(dep_norm, type = "control", control = ctrl)
  dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
  return(dep)
  
}
