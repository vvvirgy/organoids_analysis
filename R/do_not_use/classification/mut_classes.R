rm(list=ls())
.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4/')
library(tidyverse)
library(ggExtra)
library(ComplexHeatmap)
library(ggalluvial)
library(patchwork)
library(ggh4x)
library(pals)
source('organoids_analysis/R/classification/utils.R')

test_expr_rna = readRDS('data/test_expr_rna_bootstrap_v2.rds') %>% 
  mutate(Assay = 'RNA')
test_expr_prot = readRDS('data/test_expr_prot_bootstrap_v2.rds') %>% 
  mutate(Assay = 'Protein')

# checking the effect of mutation type on the class 

x = test_expr_rna
x$mut_consequence %>% unique 
cols = setNames(nm = unique(x$mut_consequence) , 
                pals::alphabet(n = length(unique(x$mut_consequence))))
x %>% 
  # filter(is_mutated == TRUE) %>% 
  filter(!is.na(cls_dosage)) %>% 
  filter(cls_dosage == 'Mutation_sensitive') %>% 
  filter(mut_consequence != 'intron_variant') %>% 
  ggplot(aes(cls_dosage, fill = mut_consequence)) + 
  geom_bar(stat = 'count', position = 'dodge')+ 
  theme_bw()+
  # light_theme() + 
  scale_fill_manual(values = cols)

dosage_colors
  

