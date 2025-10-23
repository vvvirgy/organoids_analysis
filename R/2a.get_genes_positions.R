library(tidyverse)
library(biomaRt)
source('organoids_analysis/R/functions_utils/get_genes_genomics_positions.R')

# first, extract the positions of all the genes
# transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')
# proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds') 

rna = readRDS('data/normalized_res_pseudobulk_v2.rds')
protein = readRDS('data/proteomics_normalized.rds')

genes = unique(c(rownames(rna), protein$protein))

all_genes = get_grch38_genomics_positions(genes, 'data/cnaqc/11_PDO.rds')
saveRDS(all_genes, 'data/all_genes_positions_v2.rds')

# add some metadata to include the information on cgs and intogen annotation (will be used later)

# cancer_genes_somatic_colon = readRDS('data/genes_to_check.rds')

# CGC genes
cancer_genes = read.table('data/Cancer_gene_census_tier1.csv', sep = ",", header = T) %>% 
  dplyr::as_tibble()

# genes associated with COAD
cancer_genes_somatic_colon = cancer_genes %>% 
  dplyr::filter(Somatic == 'yes') %>% 
  dplyr::filter(grepl('*colorectal*', Tumour.Types.Somatic.)) %>% 
  dplyr::select(Gene.Symbol, Role.in.Cancer) %>% 
  dplyr::distinct()

all_genes = full_join(all_genes, cancer_genes_somatic_colon, by = join_by('hgnc_symbol' == 'Gene.Symbol')) %>% 
  dplyr::rename(CGC_role_COAD = Role.in.Cancer)
# cancer_genes_somatic_colon = get_grch38_genomics_positions(cancer_genes_somatic_colon, 'data/cnaqc/11_PDO.rds')
# saveRDS(cancer_genes_somatic_colon, 'data/cancer_genes_somatic_colon_positions.rds')

# all genes from CGC
cancer_genes_somatic = cancer_genes %>% 
  dplyr::filter(Somatic == 'yes') %>% 
  # dplyr::filter(grepl('*colorectal*', Tumour.Types.Somatic.)) %>% 
  # dplyr::pull(Gene.Symbol) %>% 
  # unique
  dplyr::select(Gene.Symbol, Role.in.Cancer) %>% 
  dplyr::distinct()

all_genes = full_join(all_genes, cancer_genes_somatic, by = join_by('hgnc_symbol' == 'Gene.Symbol')) %>% 
  dplyr::rename(CGC_role_PANCANCER = Role.in.Cancer)

# cancer_genes_somatic = get_grch38_genomics_positions(cancer_genes_somatic, 'data/cnaqc/11_PDO.rds')
# saveRDS(cancer_genes_somatic, 'data/cancer_genes_somatic_positions.rds')

# COAD known drivers
intogen_drivers = read.table('../2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep = '\t', header = T)
coad_drivers = intogen_drivers %>% 
  dplyr::filter(CANCER_TYPE == 'COAD') %>% 
  dplyr::select(SYMBOL, IS_DRIVER) %>% 
  dplyr::distinct()

all_genes = full_join(all_genes, coad_drivers, by = join_by('hgnc_symbol' == 'SYMBOL')) %>% 
  dplyr::rename(is_driver_intogen = IS_DRIVER)

missing_genes = all_genes %>% 
  dplyr::filter(is.na(chr)) %>% 
  dplyr::pull(hgnc_symbol)

missing_genes = get_grch38_genomics_positions(missing_genes, 'data/cnaqc/11_PDO.rds')
missing_genes = all_genes %>% 
  dplyr::select(hgnc_symbol, CGC_role_COAD, CGC_role_PANCANCER, is_driver_intogen) %>% 
  dplyr::filter(hgnc_symbol %in% missing_genes$hgnc_symbol) %>% 
  dplyr::full_join(., missing_genes, by = 'hgnc_symbol')
all_genes = all_genes %>% 
  dplyr::filter(!hgnc_symbol %in% missing_genes$hgnc_symbol) %>% 
  bind_rows(., missing_genes)

all_genes = all_genes %>% 
  dplyr::mutate(CGC_role_COAD = ifelse(is.na(CGC_role_COAD), 'None', CGC_role_COAD), 
                CGC_role_PANCANCER = ifelse(is.na(CGC_role_PANCANCER), 'None', CGC_role_COAD), 
                is_driver_intogen = ifelse(is.na(is_driver_intogen), FALSE, is_driver_intogen))
# coad_drivers = get_grch38_genomics_positions(coad_drivers, 'data/cnaqc/11_PDO.rds')
saveRDS(all_genes, 'data/all_genes_positions_info.rds')


