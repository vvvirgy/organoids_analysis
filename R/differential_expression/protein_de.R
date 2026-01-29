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
  pivot_wider(names_from = sample, values_from = Intensity)

# trying normalising the values according to the 2n for each gene

# using metadata with removed genes that are highly fragmented
# genes_cna_status = readRDS('data/genes_cna_mut_status_filtered.rds')
genes_cna_status = readRDS('data/processed_data/genes_filtered_karyo_mut_status_filt_ccf_08.rds')
# genes_cna_status = genes_cna_status %>% 
#   tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':', convert = T, remove =F) %>% 
#   dplyr::mutate(ploidy = sum(Major + minor)) 

# selecting only some genes to run the dge --> removing genes without or with low number of diploids
diploid_genes = genes_cna_status %>% 
  filter(karyotype == '1:1') %>% 
  group_by(hgnc_symbol) %>% 
  dplyr::count() %>% 
  filter(n>2) %>% 
  pull(hgnc_symbol)

data_filt = data %>% 
  filter(PG.Genes %in% diploid_genes)

# diploid_genes = intersect(diploid_genes, rownames(data_filt))
genes = data_filt$PG.Genes

dep_by_karyotype = parallel::mclapply(genes, function(g) {
  
  message(g) 
  # create the design matrix
  design_mat = create_design(ann, genes_cna_status, gene = g, cond = 'karyotype')
  
  if(all(unique(design_mat$label) %in% colnames(data_filt)) == TRUE & length(unique(design_mat$PDO)) > 2){
    
    # create the dep object and normalise intensities
    dep = dep_preprocesing(data_filt, design_mat)
    
    if(length(unique(dep$condition)) > 1 & any(unique(dep$condition) == 'X1.1')) {    
      # run dep 
      dep_res = differential_expression(dep, ctrl = 'X1.1')
      
      list('gene' = g, 'dep' = dep_res, 'design' = design_mat)
    }
  }
}, mc.cores = 4)

# removing null elements
dep_by_karyotype = Filter(Negate(is.null), dep_by_karyotype)
names(dep_by_karyotype) = lapply(dep_by_karyotype, function(x) x$gene) %>% unlist

saveRDS(dep_by_karyotype, 'data/dep_by_karyotype_results_v2.rds')


# smad2= assay(dep_by_karyotype[[1]]$dep)
# kras= assay(dep_by_karyotype[[2]]$dep)
# 
# 
# x1 = tibble(
#   v1 = smad2['SMAD2', ], 
#   sample = smad2 %>% colnames)
# 
# x2 = tibble(v2 = kras['SMAD2', ], 
#             sample = colnames(kras))
# 
# full_join(x1, x2) %>% 
#   ggplot(aes(v1, v2)) + 
#   geom_point()

# dep_by_ploidy = lapply(genes, function(g) {
#   
#   print(g) 
#   # create the design matrix
#   design_mat = create_design(ann, genes_cna_status, gene = g, cond = 'ploidy')
#   
#   if(all(unique(design_mat$label) %in% colnames(data_filt)) == TRUE){
#     
#     # create the dep object and normalise intensities
#     dep = dep_preprocesing(data_filt, design_mat)
#     
#     if(length(unique(dep$condition > 1))) {
#       # run dep
#       dep_res = differential_expression(dep, ctrl = 'X2')
#       
#       list('gene' = g, 'dep' = dep_res, 'design' = design_mat)
#     }
#   }
# })
# 
# saveRDS(dep_by_karyotype, 'data/dep_by_ploidy_proteomics.rds')
# 
# plot_single(dep_by_karyotype[[2]]$dep, proteins = dep_by_karyotype[[2]]$gene, plot = T)
# plot_single(dep_by_ploidy[[2]]$dep, proteins = dep_by_ploidy[[2]]$gene, plot = T)
# 
# plot_single(dep_by_karyotype[[1]]$dep, proteins = dep_by_karyotype[[1]]$gene, plot = T)
# plot_single(dep_by_ploidy[[1]]$dep, proteins = dep_by_ploidy[[1]]$gene, plot = T)
# 
# SMAD2_karyo = get_df_long(dep_by_karyotype[[1]]$dep) %>% 
#   filter(PG.Genes %in% dep_by_karyotype[[1]]$gene)  
# 
# KRAS_karyo = get_df_long(dep_by_karyotype[[2]]$dep) %>% 
#   filter(PG.Genes %in% dep_by_karyotype[[2]]$gene)  
#   
# SMAD2_karyo %>% 
#   select(ends_with('CI.L'), ends_with('CI.R'), ends_with('diff'), ends_with('p.adj'), ends_with('p.val')) %>% 
#   distinct() %>% 
#   




