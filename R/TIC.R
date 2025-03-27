# new normalization method
library(limma)

data = readxl::read_excel('~/Desktop/cdslab_scratch/organoids_prj/data/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
  as.data.frame()

data[which(is.na(data[,'PG.Genes'])), 'PG.Genes'] = 'Unknown'
raw_data = data %>% 
  dplyr::select(everything(), -PG.ProteinGroups) %>% 
  tibble::column_to_rownames('PG.Genes') %>% 
  # dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric())) %>% 
  dplyr::mutate(across(where(is.character), as.numeric)) %>% 
  dplyr::mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

log_transformed = log(raw_data, base = 2)

# try to compute tic
# Compute TIC per sample
tic_factors <- colSums(log_transformed, na.rm = TRUE)

# Normalize using TIC
normalized_data <- sweep(log_transformed, 2, tic_factors, FUN = "/")

# checking and imputing NAs -- using DEP 

# creating experimental table design
exp_design = tibble(
  label = colnames(data)[-c(1,2)]
) %>% 
  mutate(condition = gsub('_a|_b', '', label)) %>% 
  mutate(replicate = str_extract(label, 'a$|b$')) 

