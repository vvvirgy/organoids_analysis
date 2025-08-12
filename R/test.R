
rna_objects_paths = list.files("data/scRNA", full.names = T)

p = rna_objects_paths[1]

df = lapply(rna_objects_paths, function(p) {
  sc.obj = readRDS(p)
  sample = unique(sc.obj@meta.data$sample)
  counts = rowSums(sc.obj@assays$RNA$counts)
  dplyr::tibble(counts = counts, gene = names(counts), sample = sample)
}) %>% bind_rows()

genes_cna_status = readRDS('data/karyotypes_mutations_all_genes_qc_ccf.rds')

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

cna_status_df = genes_cna_status %>% 
  # dplyr::group_by(sample) %>% 
  tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':', convert = T) %>% 
  dplyr::mutate(tot_cna = Major + minor) %>% 
  dplyr::rename(gene = hgnc_symbol) %>% 
  dplyr::select(gene, sample, tot_cna)


df %>% 
  dplyr::left_join(all_samples_dict %>% 
                     dplyr::select(fixed_name, scRNA_sample), by = join_by("sample" == "scRNA_sample")) %>% 
  na.omit() %>% 
  dplyr::left_join(cna_status_df, by = join_by("fixed_name" == "sample", "gene" == "gene")) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(fc = counts / sum(counts)) %>% 
  dplyr::filter(gene == "TP53") %>% 
  ggplot(mapping = aes(x = tot_cna, y = fc)) +
  geom_point()
