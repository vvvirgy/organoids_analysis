.libPaths()
rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
library(ggpmisc)
library("factoextra")
library(ggExtra)
library(patchwork)
library(jtools)

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')
# load results of the fits
rna_fit = readRDS('data/glm_fit_v2_rna.rds')
prot_fit = readRDS('data/glm_fit_v2_prot.rds')



get_clean_fit_res = function(fit_list) {
  lapply(names(fit_list), function(x) {
    fit_list[[x]] %>% 
      broom::tidy() %>% 
      dplyr::mutate(gene = x) %>% 
      dplyr::mutate(R2 = summary(fit_list[[x]])$r.squared)
      # dplyr::mutate(R2 = get_r2(fit_list[[x]]))
  }) %>% 
    dplyr::bind_rows()
}

rna_fit_all = get_clean_fit_res(rna_fit)
prot_fit_all = get_clean_fit_res(x)

