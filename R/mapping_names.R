library(tidyverse)

drivers = readRDS("drivers.rds")
proteomics_data = readxl::read_excel("DataLog&TransformedAfterNormalization_ProMeFa_NARemoved.xlsx")

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

# select the top genes to check first
genes_to_check = c("APC", "KRAS", "TP53", "PIK3CA")

drivers_to_check = drivers %>% 
  dplyr::filter(gene %in% genes_to_check)

drivers_to_check_correct_samples = drivers_to_check %>% 
  filter(PDO %in% new_mapping_experiment$fixed_name %>% unique)

# add the mapping to the proteomics ids
# new_names_data = dplyr::full_join(drivers_to_check, mapping, by = join_by("PDO" == "fixed_name"))

# CNA info 
# drivers_to_check_correct_samples %>% 
#   select(gene, karyotype, PDO) %>%
#   full_join(., new_mapping_experiment, by = join_by("PDO" == "fixed_name"))
  
# Filter the data to keep only selected genes
proteomic_to_check = proteomics_data %>%
  tidyr::separate(gene, into = c("uniprot", "gene_name"), sep = "_") %>% 
  dplyr::filter(gene_name %in% genes_to_check)

proteomic_genes = proteomic_to_check %>% 
  select(everything(), -uniprot) %>% 
  reshape2::melt()

proteomic_genes = proteomic_genes %>% 
  dplyr::rename(sample_id = variable) %>% 
  dplyr::mutate(PDO = gsub("_a|_b", "", sample_id)) %>% 
  dplyr::group_by(PDO, gene_name) %>% 
  dplyr::summarise(mean_expr = mean(value), sd_expr = sd(value)) 

samples_check = new_mapping_experiment %>% 
  select(proteomics_code, fixed_name) %>%
  distinct

x = drivers_to_check_correct_samples %>%
  select(gene, karyotype, PDO) %>% 
  distinct(.keep_all = F) %>% 
  left_join(., samples_check, by = join_by("PDO" == "fixed_name")) %>% 
  # filter(!is.na(gene)) %>% 
  full_join(., proteomic_genes, by = join_by("proteomics_code" == "PDO", "gene" == "gene_name")) %>% 
  mutate(karyotype = ifelse(is.na(karyotype), "no alteration", karyotype))


# first visualization

x = x %>% 
  mutate(karyotype = factor(karyotype, 
                            levels = c("no alteration", "1:1", "1:0", "2:0", "2:1", "2:2", "3:0", "3:1", "3:2", "3:3", "4:0", "5:0")))

cna_colors = setNames(nm = x$karyotype %>% levels, object = c("gainsboro", "#BBE2EC", yarrr::piratepal(palette = "basel")))

x %>% 
  ggplot(aes(x = karyotype, y = mean_expr, fill = karyotype)) + 
  geom_boxplot(outliers = F) + 
  facet_wrap(vars(gene), scales = "free") + 
  geom_jitter() + 
  theme_bw() + 
  scale_fill_manual(values = cna_colors)

