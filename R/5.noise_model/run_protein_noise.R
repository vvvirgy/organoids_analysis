rm(list=ls())

library(tidyverse)
library(DEP)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')

source("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/constants.R")
source("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/make_groups.R")
source('/orfeo/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/proteomics_de_utils.R')

IMG_PATH = file.path(res_path, 'noise_model')
RES_PATH = file.path(data_path, 'noise_model')

dir.create(IMG_PATH, recursive = T)
dir.create(RES_PATH, recursive = T)

# using metadata with removed genes that are highly fragmented
# genes_cna_status = readRDS(file.path(data_path, 'processed_data/dna/genes_filtered_karyo_mut_status_filt_ccf_08.rds'))

# get correct sample names
dict = readRDS(file.path(data_path, 'full_dict_dna_rna_prot.rds'))

design_matrix = readRDS(file.path(RES_PATH, 'filtered_design_matrix.rds'))

# get the pair proteomics-genomics correct sample name
samples_check = dict %>% 
  dplyr::select(proteomics_code, DNA, fixed_name) %>%
  dplyr::distinct()

data = readxl::read_excel(file.path(data_path, 'proteomics/Results_Organoids_NoIsoforms.xlsx'), sheet = 1) %>% 
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


# selecting only some genes to run the dge --> removing genes without or with low number of diploids
tetraploid_genes = names(design_matrix)


dep_by_karyotype = parallel::mclapply(tetraploid_genes, function(g) {

    message(g)
    des = design_matrix[[g]] %>% 
      dplyr::rename(hgnc_symbol = gene)
    
    # create the design matrix
    iterations = des$iteration %>% unique
    
    # run for n iterations 
    res = lapply(iterations, function(i) {
      
      message(i)
      
      des = des %>% 
        filter(iteration == i)
      
      des_matrix = full_join(des, ann, 
                             by = join_by('sample' == 'fixed_name')) %>% 
        filter(!is.na(group)) %>% 
        filter(!is.na(PDO)) 
      
      design_mat = des_matrix %>% 
        group_by(group) %>%
        rename(replicate_name = replicate) %>% 
        filter(!is.na(replicate_name)) %>%
        mutate(replicate = seq_len(n())) %>% 
        dplyr::rename(label = replicate_name) %>%
        rename(condition = group) 
      
      data_gene = data %>% 
        filter(PG.Genes == g) %>% 
        dplyr::select(design_mat$label)
      
      if(nrow(data_gene) > 1 | sum(is.na(data_gene)) < (ncol(data_gene)/2)) {
        
        data_filt = data %>% 
          dplyr::select(PG.ProteinGroups, PG.Genes, design_mat$label)
        
        # create the dep object and normalise intensities
        dep = dep_preprocesing(data_filt, design_mat)
        
        if(length(unique(dep$condition)) > 1) {    
          # run dep 
          dep_res = differential_expression(dep, ctrl = 'A')
          
          list('dep' = dep_res, 'design' = design_mat, 'iteration' = i)
        }
      } 
      
    })
    names(res) = paste('iteration', lapply(res, function(r) {r$iteration}) %>% unlist, sep = '_')
    res = c('iterations' = list(res), 'gene' = g)
    return(res)
    
  }, mc.cores = 4)

# removing null elements
dep_by_karyotype = Filter(Negate(is.null), dep_by_karyotype)
names(dep_by_karyotype) = lapply(dep_by_karyotype, function(x) x$gene) %>% unlist

# dir.create(file.path(data_path, 'noise_model/protein'), recursive = T)

saveRDS(dep_by_karyotype, file.path(data_path, 'noise_model/protein/dep_tetraploid_genes_results.rds'))




