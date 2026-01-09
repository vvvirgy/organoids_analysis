rm(list=ls())
library(tidyverse)

x = tibble(sample = 
             list.files('data/scRNA/') %>% gsub('_filtered.rds', '', .))

genomic = tibble(sample = 
                   list.files('data/cnaqc/') %>% gsub('.rds', '', .))

# write.table(x, 'data/scRNA_samples.csv', sep = ',', quote = F, col.names = T)
write.table(genomic, 'data/genomics_samples.csv', sep = ',', quote = F, col.names = T)

x = x %>% 
  mutate(sample = gsub('Sample_', '', sample))

tt = read.csv('data/scRNA_samples.csv', header = T, sep = ',')

tt = tt %>% 
  select(RNA, DNA) %>% 
  mutate(RNA = gsub('Sample_', '', RNA)) %>% 
  mutate(RNA = ifelse(RNA == '', NA, RNA))

write.csv(tt, file = 'data/rna_dna_dictionary.csv', sep = ',', quote = F, col.names = T)

# adding proteomics
# proteomics_data = readxl::read_excel("~/Documents/Università/PhD/projects/organoids/DataLog&TransformedAfterNormalization_ProMeFa_NARemoved.xlsx")
proteomics_data = readxl::read_excel("data/DataLog&TransformedAfterNormalization_ProMeFa_NARemoved.xlsx")

proteomics_data = proteomics_data %>% 
  dplyr::rename("gene" = ...1)

# map sample names among the experiments (hopefully)
# mapping = readxl::read_excel("~/Documents/Università/PhD/projects/organoids/mapping_ICR_names.xlsx")
mapping = readxl::read_excel('data/mapping_ICR_names.xlsx')

mapping = mapping %>% 
  select(`Patient code (HSR)`, `Patient code (ICR)`, fixed_name) %>% 
  dplyr::rename(proteomics_code = `Patient code (HSR)`) %>% 
  dplyr::rename(genomics_code = `Patient code (ICR)`) %>% 
  dplyr::mutate(proteomics_code = gsub("#", "", proteomics_code))

mapping_bis = tibble(
  sample_proteomics_code = colnames(proteomics_data)[-1], 
  genomics_code = gsub("_a|_b", "", colnames(proteomics_data)[-1])
  # fixed_name = colnames(proteomics_data)[-1]
)

proteomics_names = full_join(mapping, mapping_bis, by = join_by("proteomics_code" == "genomics_code"))

new_mapping = proteomics_names %>% 
  mutate(fixed_name = ifelse(is.na(fixed_name), gsub("HSR", "", proteomics_code), fixed_name)) %>% 
  mutate(fixed_name = gsub("-", "_", fixed_name)) 

new_mapping_experiment = new_mapping %>% 
  # dplyr::filter(!is.na(sample_proteomics_code)) %>% 
  dplyr::mutate(fixed_name = gsub("_seg", "", fixed_name)) %>% 
  dplyr::mutate(fixed_name = case_when(fixed_name == "65_VI" ~ "65_IV", 
                                       .default = fixed_name))

saveRDS(new_mapping_experiment, "mapping_samples.rds")

tt = tt %>% 
  mutate(DNA = gsub(' ', '', DNA)) %>% 
  mutate(DNA_v2 = gsub('_PDO$', '', DNA))

full_dict = new_mapping_experiment %>% 
  full_join(., tt, by = join_by('fixed_name' == 'DNA_v2'))

saveRDS(full_dict, 'data/full_dict_dna_rna_prot.rds')


