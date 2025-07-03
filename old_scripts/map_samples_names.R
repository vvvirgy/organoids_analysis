library(tidyverse)

drivers = readRDS("drivers.rds")
proteomics_data = readxl::read_excel("DataLog&TransformedAfterNormalization_ProMeFa_NAImputed.xlsx")

proteomics_data = proteomics_data %>% 
  dplyr::rename("gene" = ...1)

# map sample names among the experiments (hopefully)
mapping = readxl::read_excel("mapping_ICR_names.xlsx")

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

tt = full_join(mapping, mapping_bis, by = join_by("proteomics_code" == "genomics_code"))

new_mapping = tt %>% 
  mutate(fixed_name = ifelse(is.na(fixed_name), gsub("HSR", "", proteomics_code), fixed_name)) %>% 
  mutate(fixed_name = gsub("-", "_", fixed_name)) 

new_mapping_experiment = new_mapping %>% 
  dplyr::filter(!is.na(sample_proteomics_code)) %>% 
  dplyr::mutate(fixed_name = gsub("_seg", "", fixed_name)) %>% 
  dplyr::mutate(fixed_name = case_when(fixed_name == "65_VI" ~ "65_IV", 
                                       .default = fixed_name))

saveRDS(new_mapping_experiment, "mapping_samples.rds")



