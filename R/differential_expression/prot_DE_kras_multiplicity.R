library(CNAqc)
library(tidyverse)
# library(Seurat)
library(DEP)
source('organoids_analysis/R/proteomics_de_utils.R')

# cnas = readRDS('data/cnaqc_v2/cnas_list_v2.rds')

# kras_status = lapply(cnas, function(x) {
#   x$phasing %>% 
#     filter(VEP.SYMBOL == 'KRAS') %>% 
#     dplyr::select(VEP.SYMBOL, driver_label, VAF, CCF, karyotype, multiplicity, sample)
# }) %>% bind_rows()

# using only KRAS mutated samples with associated CCF
# kras_status = lapply(cnas %>% names, function(x) {
#   cnas[[x]]$phasing %>% 
#     dplyr::filter(VEP.SYMBOL == 'KRAS') %>% 
#     dplyr::filter(is_driver) %>% 
#     dplyr::select(VEP.SYMBOL, driver_label, VAF, CCF, karyotype, multiplicity) %>% 
#     mutate(sample = x)
# }) %>% 
#   bind_rows()
# kras_status$sample %>% unique %>% length()
# 
# mut_samples = kras_status$sample %>% unique

# samples_to_use = samples_to_use[which(sapply(samples_to_use, nrow) > 0)] %>% names

# get correct sample names
# dict = readRDS('data/full_dict_dna_rna_prot.rds')
# wt_samples = dict %>% 
#   filter(!fixed_name %in% mut_samples) %>% 
#   pull(fixed_name) %>% 
#   unique


dict = readRDS('data/full_dict_dna_rna_prot.rds')
genes_cna_status = readRDS('data/processed_data/karyotypes_mutations_all_genes_qc_ccf.rds')

samples_check = dict %>% 
  dplyr::select(proteomics_code, fixed_name) %>%
  dplyr::distinct()

kras_status = genes_cna_status %>% 
  filter(hgnc_symbol == 'KRAS') %>% 
  dplyr::mutate(mut_consequence = ifelse(is.na(mut_consequence), 'wild-type', mut_consequence)) %>% 
  tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':', convert = T) %>% 
  mutate(k = Major + minor) %>% 
  mutate(multiplicity = ifelse(is.na(multiplicity), 0, multiplicity)) %>% 
  mutate(mut_ratio = paste0(multiplicity, '/', k)) 

kras_status_tb = kras_status %>% 
  mutate(IMPACT = ifelse(is.na(IMPACT), 'wild-type', IMPACT)) %>% 
  mutate(IMPACT = factor(IMPACT, levels = c('wild-type','LOW', 'MODIFIER', 'MODERATE', 'HIGH'))) %>% 
  group_by(sample) %>% 
  slice_max(IMPACT) %>% 
  dplyr::select(hgnc_symbol, sample, multiplicity, mut_ratio) %>% 
  # mutate(m_class = case_when(
  #   is.na(multiplicity) ~ 'wild-type', 
  #   multiplicity == 1 ~ 'single', 
  #   multiplicity > 1 ~ 'multiple'
  # )) %>% 
  filter(sample != '11') %>% 
  group_by(sample) %>% 
  # dplyr::select(-mut_consequence) %>% 
  distinct()
saveRDS(kras_status_tb, 'data/kras_multiplicity/design_matrix.rds')

kras_status_tb = readRDS('data/kras_multiplicity/design_matrix.rds')

# get the pair proteomics-genomics correct sample name
samples_check = dict %>% 
  dplyr::select(proteomics_code, DNA, fixed_name) %>%
  dplyr::distinct()

data = readxl::read_excel('data/utilities/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
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
  pivot_wider(names_from = sample, values_from = Intensity)



# trying normalising the values according to the 2n for each gene

# using metadata with removed genes that are highly fragmented
# genes_cna_status = readRDS('data/genes_cna_mut_status_filtered.rds')

# genes_cna_status = genes_cna_status %>% 
#   tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':', convert = T, remove =F) %>% 
#   dplyr::mutate(ploidy = sum(Major + minor)) 

# selecting only some genes to run the dge --> removing genes without or with low number of diploids


data_filt = data 
# diploid_genes = intersect(diploid_genes, rownames(data_filt))
g = 'KRAS'

# kras_status_test = kras_status_tb %>% 
#   filter(sample %in% tt)
# create the design matrix
design_mat = create_design(ann, kras_status_tb, gene = g, cond = 'mut_ratio') %>% 
  dplyr::select(label, condition, replicate)

if(all(unique(design_mat$label) %in% colnames(data_filt)) == TRUE & length(unique(design_mat$PDO)) > 2){
  
  # create the dep object and normalise intensities
  dep = dep_preprocesing(data = data_filt, design = design_mat)
  
  if(length(unique(dep$condition)) > 1 & any(unique(dep$condition) == 'X0.2')) {    
    # run dep 
    dep_res = differential_expression(dep, ctrl = 'X0.2')
    
    list('gene' = g, 'dep' = dep_res, 'design' = design_mat)
  }
}

get_dge_res = function(dep, gene) {
  
  dep %>% 
    get_df_long() %>% 
    filter(PG.Genes == gene) 
  
}

kras_res = get_dge_res(dep_res, 'KRAS')
kras_res_red = kras_res %>% 
  select(PG.Genes, condition, ends_with('CI.L'), ends_with('CI.R'), ends_with('diff'), 
         ends_with('p.adj'), ends_with('p.val'), ends_with('significant'))
kras_res_red = kras_res_red %>% 
  distinct()

get_fc_tb_clean = function(tb) {
  
  lapply(tb$condition %>% unique, function(x) {
    df = tb %>% 
      select(PG.Genes, starts_with(x)) %>%
      distinct()
    colnames(df) = gsub(paste0(x, '_vs_X0.2_'), '', colnames(df))
    
    df = df %>% 
      mutate(condition = x) %>% 
      relocate(condition, .after = PG.Genes) %>% 
      filter(condition != 'X0.2')
    
    return(df)
    
  }) %>% 
    bind_rows()
  
}

fc_tb_clean = get_fc_tb_clean(kras_res_red)
saveRDS(fc_tb_clean, 'data/kras_multiplicity/fc_tb_clean_prot.rds')

fc_tb_clean %>% 
  filter(condition != 'X0') %>% 
  mutate(mut_ratio = )
  # mutate(m_status = ifelse(condition == 'single', 'm = 1', 'm > 1')) %>% 
  mutate(mut_ratio = gsub('X', '', condition),
         mut_ratio = gsub('\\.', '/', mut_ratio)) %>% 
  mutate(adj.pval = p.adjust(p.val, method = 'BH')) %>% 
  mutate(significance = ifelse(adj.pval <= .05, 'significant (adjP <= 0.05', 'ns')) %>% 
  ggplot(aes(
    y = mut_ratio, 
    x = diff, 
    fill = significance
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  labs(
    x = 'log2FC', 
    y = 'm/k KRAS'
  )

