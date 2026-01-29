rm(list=ls())

.libPaths()
library(tidyverse)
library(DEP)

# setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
source('organoids_analysis/R/functions_utils/proteomics_utils.R')
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

data = data %>% 
  pivot_wider(names_from = sample, values_from = Intensity)

# # add annotations to the data
# data = data %>% 
#   dplyr::full_join(., ann, by = join_by('sample' == 'replicate')) %>%
#   filter(!is.na(sample))

# create the dep object
data = make_unique(data, names = 'PG.Genes', ids = 'PG.ProteinGroups', delim = '_')
samples_index = which(!colnames(data) %in% c('PG.ProteinGroups', 'PG.Genes', 'name', 'ID'))

ann = ann %>% 
  group_by(PDO) %>%
  rename(replicate_name = replicate) %>% 
  mutate(replicate = str_extract(replicate_name, 'a$|b$')) %>% 
  mutate(replicate = factor(replicate)) %>% 
  mutate(replicate = as.numeric(replicate)) 
ann = ann %>% 
  mutate(condition = 'HSR') %>% 
  rename(label = replicate_name)

# data_se <- make_se(data, 2:73, ann)
dep = DEP::make_se_parse(proteins_unique = data, columns = samples_index)





