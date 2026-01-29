rm(list=ls())

# .libPaths()
library(tidyverse)
library(DEP)

source('organoids_analysis/R/functions_utils/fit_plots.R')
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

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

# removing missing values by default
data %>% 
  group_by(PG.Genes) %>% 
  filter(!is.na(Intensity)) %>% 
  mutate(n_samples = n()) %>% 
  filter(n_samples == 72) 
  # pull(PG.Genes) %>% unique %>% length
  # filter(PG.Genes == 'TP53')
  
  
  ggplot(aes(n_samples)) + 
  geom_bar(stat = 'count')


data %>% 
  filter(!is.na(Intensity)) %>% 
  ggplot(aes(Intensity, fill = sample)) + 
  geom_histogram() 

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

# trying normalising the values according to the 2n for each gene
genes_cna_status = readRDS('data/karyotypes_mutations_all_genes_qc_ccf_v4.rds')

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
  # mutate(condition = 'HSR') %>% 
  rename(label = replicate_name) %>%
  filter(!is.na(PDO)) %>%
  filter(!is.na(label)) %>% 
  rename(condition = PDO)

# from now values are log transformed
dep = DEP::make_se(proteins_unique = data, columns = samples_index, expdesign = ann)

# filtering out proteins with missing values in both replicates
dep_filt = filter_missval(dep, thr = 0)

plot_numbers(dep_filt)

# visualise the expression

p1 = data %>% 
  select(-c(ID, name)) %>% 
  pivot_longer(cols = colnames(data)[-c(1,2,75, 76)], names_to = 'sample', values_to = 'expression') %>% 
  ggplot(aes(y = sample, x = log(expression, 2), fill = sample)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = 'None')
  
(plot_normalization(dep) + theme(legend.position = 'None')) / p1
ggsave('res/test_norm_proteomics.pdf', width = 10, height = 20)

plot_normalization(dep_filt)

dep_norm = normalize_vsn(dep_filt)

plot_normalization(dep_norm, dep_filt)


prot_norm = dep_norm@assays@data@listData[[1]]
prot_norm = prot_norm %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('gene') %>% 
  pivot_longer(cols = colnames(prot_norm), names_to = 'sample', values_to = 'expression')

metadata = dep_norm@colData %>% 
  as_tibble()

prot_norm = prot_norm %>% 
  full_join(., metadata, by = join_by('sample' == 'ID'))


# adding further metadata
genes_to_check = genes_cna_status$hgnc_symbol %>% unique

# filter out genes that are highly fragmented
coad_genes = readRDS('data/all_genes_positions_info.rds')
genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol)) %>% 
  dplyr::select(hgnc_symbol, chr, from, to)

genes_cna_status = genes_cna_status %>% 
  dplyr::mutate(mut_consequence = ifelse(is.na(mut_consequence), 'wild-type', mut_consequence))

# removing mutations falling on the same segment + multihit mutations + selecting the segments with the maximum overlap over the gene 
drivers_to_check_correct_samples = filter_fragmented_cnas(genes_cna_status, 
                                                          samples_list = (samples_check$fixed_name %>% unique), 
                                                          genes_to_check = genes_to_check, 
                                                          min_length = 10^6, 
                                                          genes_position = genes, 
                                                          strategy = 'MOv') %>% 
  mutate(multiplicity = ifelse(mut_consequence == 'wild-type', 0, multiplicity))

saveRDS(drivers_to_check_correct_samples, 'data/genes_cna_mut_status_filtered.rds')



