library(tidyverse)
library(CNAqc)
.libPaths()

rm(list=ls())

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj')
cnas = readRDS('data/cnaqc_v2/cnas_list_v2.rds')

sequenced = 3.1e9

FGA = lapply(cnas, function(x) {
  fga = x$basepairs_by_karyotype %>% 
    filter(karyotype != "1:1") %>%
    pull(n) %>% 
    sum()/sequenced
  tb = tibble(sample = x$sample, FGA = fga)
  return(tb)
}) %>% bind_rows()

FGA = FGA %>% 
  filter(sample != '11')

summary(FGA$FGA)

karyo_prop = lapply(cnas, function(x) {
  tb = x$basepairs_by_karyotype %>% 
    filter(karyotype != "1:1") %>%
    group_by(karyotype) %>% 
    summarise(tot = sum(n)) %>% 
    mutate(prop = tot/sequenced) %>% 
    filter(karyotype %in% c('1:0', '2:0', '2:1', '2:2')) %>% 
    mutate(sample = x$sample)
  return(tb)
}) %>% bind_rows()

karyo_prop = karyo_prop %>% 
  filter(sample != '11')

karyo_prop %>% 
  group_by(karyotype) %>% 
  mutate(prop = prop * 100) %>% 
  summarise(mean = mean(prop), min = min(prop), max = max(prop), sd = sd(prop)) 



