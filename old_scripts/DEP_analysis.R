library(DEP)
library(tidyverse)

data = readxl::read_excel('~/Desktop/cdslab_scratch/organoids_prj/data/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
  as.data.frame()

data[which(is.na(data[,'PG.Genes'])), 'PG.Genes'] = 'Unknown'

colnames(data)

exp_design = tibble(
  label = colnames(data)[-c(1,2)]
) %>% 
  mutate(condition = gsub('_a|_b', '', label)) %>% 
  mutate(replicate = str_extract(label, 'a$|b$')) 

data = data %>% 
  dplyr::mutate(across(-c(PG.Genes, PG.ProteinGroups), as.numeric)) %>% 
  mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

data = make_unique(data, names = 'PG.Genes', ids = 'PG.ProteinGroups', delim = '_')
samples_index = which(!colnames(data) %in% c('PG.ProteinGroups', 'PG.Genes', 'name', 'ID'))

dep = DEP::make_se_parse(proteins_unique = data, columns = samples_index)

data_filt <- filter_missval(dep, thr = 0)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

plot_normalization(data_filt, data_norm)

data_imp <- impute(data_norm, fun = "MLE")
plot_imputation(data_norm, data_imp)

imputed_data = DEP::get_df_long(data_imp)
imputed_data = imputed_data %>% 
  dplyr::mutate(label = gsub('X', '', label)) 

DEP::get_df_long(data_imp) %>% 
  filter(name == 'KRAS') %>% 
  ggplot(aes(intensity))+
  geom_histogram(binwidth = 0.1)

# remap sample names to the original ones!

