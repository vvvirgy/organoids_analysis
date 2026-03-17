rm(list=ls())

# .libPaths()
library(tidyverse)
library(DEP)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
source('organoids_analysis/R/functions_utils/fit_plots.R')
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')
source('organoids_analysis/R/proteomics_de_utils.R')

# setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
# source('organoids_analysis/R/functions_utils/proteomics_utils.R')
# source('organoids_analysis/R/functions_utils/fit_plots.R')
# source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

# get correct sample names
dict = readRDS('data/full_dict_dna_rna_prot.rds')

# get the pair proteomics-genomics correct sample name
samples_check = dict %>% 
  dplyr::select(proteomics_code, DNA, fixed_name) %>%
  dplyr::distinct()

data = readxl::read_excel('data/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
  as.data.frame()

data[which(is.na(data[,'PG.Genes'])), 'PG.Genes'] = 'Unknown'

# processing data
data = data %>% 
  # dplyr::select(-PG.ProteinGroups) %>% 
  tidyr::pivot_longer(., cols = colnames(data)[-(1:2)], 
                      names_to = 'sample', 
                      values_to = 'Intensity') %>% 
  dplyr::mutate(Intensity = as.numeric(Intensity)) %>% 
  dplyr::mutate(Intensity = ifelse(is.nan(Intensity), NA, Intensity)) 

# creating annotation
ann = tibble(
  PDO = data$sample,
  replicate = data$sample, 
) %>% 
  distinct() %>% 
  mutate(PDO = gsub('_a$|_b$', '', PDO)) %>% 
  mutate(PDO = gsub('_HSR', 'HSR', PDO)) %>% 
  full_join(., samples_check, by = join_by('PDO' == 'proteomics_code')) 

ann = ann %>% 
  filter(!is.na(PDO)) %>% 
  filter(!is.na(replicate))

data = data %>% 
  full_join(., ann, by = join_by('sample' == 'replicate')) %>%
  filter(!is.na(Intensity))

genes_cna_status = readRDS('data/processed_data/genes_filtered_karyo_mut_status_filt_ccf_08.rds')
cnas = readRDS('data/cnaqc_v2/cnas_list_v2.rds')

ploidies = lapply(cnas, function(x) {
  tibble(
    sample = x$sample, 
    ploidy = x$ploidy
  )
}) %>% bind_rows()

saveRDS(ploidies, 'data/ploidies_pdos.rds')

ploidies = readRDS('data/ploidies_pdos.rds')

data = data %>% 
  left_join(., ploidies, by = join_by('fixed_name' == 'sample'))

data = data %>% 
  rename(sample_ploidy = ploidy) %>% 
  left_join(., genes_cna_status, by = join_by('fixed_name' == 'sample', 
                                              'PG.Genes' == 'hgnc_symbol')) %>% 
  filter(!is.na(karyotype))

homeostatis_genes = clusterProfiler::read.gmt('data/GOBP_HOMEOSTATIC_PROCESS.v2026.1.Hs.gmt')
homeostatis_genes = homeostatis_genes$gene %>% unique

data %>% 
  group_by(PDO, PG.Genes) %>% 
  mutate(mean_int = mean(Intensity)) %>% 
  # summarise(sf = sum(Intensity)) %>% 
  # filter(PG.Genes == 'KRAS') %>%
  filter(PG.Genes %in% homeostatis_genes) %>% 
  # mutate(logIntensity = log(Intensity, base = 2)) %>% 
  ggplot(aes(
    x = karyotype, 
    y = log1p(mean_int)
  )) + 
  geom_violin()

data %>% 
  group_by(PDO, PG.Genes) %>% 
  mutate(mean_int = mean(Intensity)) %>% 
  # summarise(sf = sum(Intensity)) %>% 
  # filter(PG.Genes == 'KRAS') %>%
  filter(!PG.Genes %in% homeostatis_genes) %>% 
  # mutate(logIntensity = log(Intensity, base = 2)) %>% 
  ggplot(aes(
    x = karyotype, 
    y = log1p(mean_int)
  )) + 
  geom_violin()


