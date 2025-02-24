library(tidyverse)

drivers = readRDS("drivers.rds")
proteomics_data = readxl::read_excel("DataLog&TransformedAfterNormalization_ProMeFa_NARemoved.xlsx")

proteomics_data = proteomics_data %>% 
  dplyr::rename("gene" = ...1)

# map sample names among the experiments (hopefully)
mapping = readxl::read_excel("mapping_ICR_names.xlsx")

mapping = mapping %>% 
  select(`Patient code (HSR)`, `Patient code (ICR)`, fixed_name)

# select the top genes to check first
genes_to_check = c("APC", "KRAS", "TP53", "PIK3CA")

drivers_to_check = drivers %>% 
  dplyr::filter(gene %in% genes_to_check)

# add the mapping to the proteomics ids
new_names_data = dplyr::full_join(drivers_to_check, mapping, by = join_by("PDO" == "fixed_name"))

# CNA info 
new_names_data %>% 
  select(gene, karyotype, PDO) %>% 
  full_join(., proteomic_genes, by = "PDO")

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
  
# first visualization
proteomic_genes %>% 
  ggplot(aes(x = reorder(gene_name, +mean_expr), y = mean_expr, fill = reorder(PDO, +mean_expr))) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~gene_name, scales = "free") 

# add the CNA info
drivers_to_check$PDO %>% unique() %>% length

samples_proteomics = lapply(strsplit(x = proteomic_genes$variable %>% unique %>% as.character, split = "_"), dplyr::first) %>% 
  unlist() %>% unique %>% 
  gsub("-", "_", .)



