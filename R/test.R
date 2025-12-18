library(clusterProfiler)

data_path = 'data/scRNA'
data = lapply(list.files(data_path, full.names = T), function(x) {
  readRDS(x)
})

matching_samples = readxl::read_excel('data/scRNA_samples.xlsx', sheet = 1) %>% 
  as.data.frame() %>% 
  dplyr::mutate(genomics_code = ifelse(genomics_code == 'NA', NA, genomics_code)) %>% 
  dplyr::mutate(fixed_name = ifelse(fixed_name == 'NA', NA, fixed_name))

samples_check = matching_samples %>% 
  dplyr::select(fixed_name, scRNA_sample) %>%
  dplyr::distinct()

new_mapping_experiment = readRDS("data/mapping_samples.rds")

all_samples_dict = full_join(matching_samples, new_mapping_experiment, 
                             by = join_by('genomics_code' == 'genomics_code', 
                                          'fixed_name' == 'fixed_name'), na_matches = 'na')

genes_karyo = readRDS('/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_all_genes_qc_v3.rds')
genes_karyo = genes_karyo %>% 
  mutate(karyotype = gsub(' ', '', karyotype))

data = lapply(list.files(data_path, full.names = T), function(p) {
  
  x = readRDS(p)
  sample = x$sample %>% unique
  
  x = x@assays$RNA$counts
  x = as.matrix(x)
  
  gene_expression = rowSums(x)
  gene_expression = gene_expression / sum(gene_expression)
  
  genomic_sample = all_samples_dict %>% 
    filter(scRNA_sample == sample)
  
  if(nrow(genomic_sample > 0)) {
    genomic_sample = genomic_sample$fixed_name %>% unique
    
    karyo = genes_karyo %>% filter(sample == genomic_sample)
    
    dplyr::tibble(hgnc_symbol = names(gene_expression), expr = gene_expression) %>% 
      dplyr::left_join(karyo, by = "hgnc_symbol") %>% 
      na.omit()
  }
}) %>% 
  bind_rows()


housekeeping = clusterProfiler::read.gmt('HSIAO_HOUSEKEEPING_GENES.v2025.1.Hs.gmt')
housekeeping = housekeeping$gene %>% unique

data_by_sample = data %>% 
  group_by(sample) %>% 
  group_split()

s1 = data_by_sample[[1]]

housekeeping_expr = s1 %>% 
  # (housekeeping = ifelse(hgnc_symbol %in% housekeeping, 'yes', 'no')) %>% 
  # group_by(housekeeping) %>% 
  filter(hgnc_symbol %in% housekeeping) %>% 
  summarise(housekeeping_expression_mean = mean(expr), housekeeping_exp_sd = sd(expr)) %>% 
  mutate(lower_bound = housekeeping_expression_mean - housekeeping_exp_sd, upper_bound = housekeeping_expression_mean + housekeeping_exp_sd) %>% 
  pivot_longer(cols = c(lower_bound, upper_bound), values_to = 'values', names_to = 'bound')


s1 %>% 
  ggplot(aes(expr)) + 
  geom_histogram() + 
  theme_bw() + 
  geom_vline(data = housekeeping_expr, 
             aes(xintercept = values), inherit_aes = F) + 
  scale_x_log10()

# names(data) = c(45, 55)
# 
# 
# 
# x = lapply(data %>% names, function(n){
#   
#   
#   
#   karyo = genes_karyo %>% filter(sample == as.integer(n))
#   
#   dplyr::tibble(hgnc_symbol = names(gene_expression), expr = gene_expression) %>% 
#     dplyr::left_join(karyo, by = "hgnc_symbol") %>% 
#     na.omit()
# }) %>% 
#   bind_rows()

# data = data %>% 
#   separate(karyotype, sep = ':', into = c('Major', 'minor'), convert = T) %>% 
#   mutate(tot_cna = mutate(Major + minor))

data %>% 
  filter(hgnc_symbol == 'AKT1') %>% 
  ggplot(mapping = aes(x = karyotype, y = expr)) +
  # geom_point() + 
  geom_boxplot() +
  # ggpubr::stat_compare_means(comparisons = list(c("1:0", "1:1"), c("2:0", "2:1"), c("1:1", "2:0"), c('1:1', '2:1'), c('2:1', '2:2'))) +
  ggpubr::stat_compare_means(comparisons = list(c('1:1', "2:2"))) + 
  # facet_wrap(~sample, scales = 'free_x') + 
  scale_y_log10()


data %>% 
  filter(hgnc_symbol == 'PCM1') %>% 
  ggplot(aes(karyotype)) + 
  geom_bar(stat = 'count')
# 


# libsize = colSums(x)
# gene_expression = (x/libsize) %>% rowSums()




dplyr::tibble(hgnc_symbol = names(gene_expression), expr = gene_expression) %>% 
  dplyr::left_join(genes_karyo, by = "hgnc_symbol") %>% 
  na.omit() %>% 
  ggplot(mapping = aes(x = karyotype, y = log(expr))) +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = list(c("1:0", "1:1"), c("2:0", "2:1"), c("1:1", "2:0")))


DESeq2:: 