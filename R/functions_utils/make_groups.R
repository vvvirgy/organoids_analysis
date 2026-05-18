make_groups = function(cna_data, 
                       genes, 
                       karyo,
                       n = 100 # number of times to sample the groups 
                       ) {
  samples_classification = lapply(genes, function(g) {
    sample_list = cna_data %>% 
      filter(hgnc_symbol == g) %>% 
      filter(karyotype == karyo) %>% 
      pull(sample)
    
    
    lapply(1:n, function(i) {
      
      sample_groups = tibble(
        sample = sample(sample_list, size = round(length(sample_list)/2)), 
        group = 'A'
      ) 
      
      sample_groups = sample_groups %>% 
        bind_rows(., 
                  tibble(
                    sample = setdiff(sample_list, sample_groups$sample), 
                    group = 'B'
                  )) %>% 
        mutate(iteration = i) %>% 
        mutate(gene = g)
      }) %>% 
      bind_rows()
  })
  names(samples_classification) = genes
  return(samples_classification)
}



