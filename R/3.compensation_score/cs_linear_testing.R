library(tidyverse)
library(ggpmisc)

source('organoids_analysis/R/functions_utils/constants.R')

multi_omics = "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/CS_scores_prot_and_rna.rds"
multi_omics = readRDS(multi_omics)

multi_omics_red = multi_omics %>% 
  group_by(name, omic) %>% 
  filter(n() > 2)

multi_omics_red %>% 
  # filter(name == 'SERBP1') %>% 
  mutate(karyotype = factor(karyotype, levels = c('1:0', '2:0', '2:1', '2:2')),
         CS = as.numeric(CS)) %>% 
  mutate(karyo_num = as.numeric(karyotype)) %>% 
  ggplot(aes(x = karyo_num, 
             y = CS, 
             color = omic)) + 
  geom_point() + 
  geom_line() +
  facet_wrap(~name) + 
  theme_light() + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey70') + 
  scale_x_continuous(labels = setNames(object = c('1:0', '2:0', '2:1', '2:2'), seq(1:4))) + 
  ggpmisc::stat_poly_eq(use_label('eq'))
  
