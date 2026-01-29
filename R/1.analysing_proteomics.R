rm(list=ls())

library(tidyverse)
library(prolfqua)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
source('organoids_analysis/R/functions_utils/proteomics_utils.R')
# source('organoids_analysis/R/functions_utils/fit_plots.R')
# source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

# get correct sample names
dict = readRDS('data/full_dict_dna_rna_prot.rds')

# get the pair proteomics-genomics correct sample name
samples_check = dict %>% 
  dplyr::select(proteomics_code, DNA) %>%
  dplyr::distinct()

data = readxl::read_excel('data/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
  as.data.frame()

data[which(is.na(data[,'PG.Genes'])), 'PG.Genes'] = 'Unknown'

# processing data
data = data %>% 
  dplyr::select(-PG.ProteinGroups) %>% 
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

# add annotations to the data
data = data %>% 
  dplyr::full_join(., ann, by = join_by('sample' == 'replicate')) %>%
  filter(!is.na(sample))


# check on raw intensities 
data %>% 
  ggplot(aes(Intensity, fill = DNA)) + 
  geom_histogram() +
  # geom_density() + 
  theme_bw()

# creating AnalysisTableAnnotation obj
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "sample"
atable$workIntensity = "Intensity"
atable$hierarchy[['PG.Genes']] = 'PG.Genes'
atable$factors[['PDO']] = 'PDO'

config <- prolfqua::AnalysisConfiguration$new(atable)
analysis_data <- prolfqua::setup_analysis(data, config)
lfqdata <- prolfqua::LFQData$new(analysis_data, config)

# remove small intensities values 
lfqdata$remove_small_intensities()

# visualize the intensities
smrz <- lfqdata$get_Summariser()
smrz$plot_hierarchy_counts_sample()

lfqplotter <- lfqdata$get_Plotter()
density_nn <- lfqplotter$intensity_distribution_density()
density_nn + 
  theme(legend.position = 'bottom')

lfqplotter$pca()

# visualize the NAs
lfqplotter$NA_heatmap()
lfqdata$get_Summariser()$plot_missingness_per_group()

# visualize the coefficient of variation
stats <- lfqdata$get_Stats()
stats$violin()
stats$density_median()

stdm_raw <- stats$stdv_vs_mean(size = 10000) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10()
stdm_raw

# perform normalization 
lt <- lfqdata$get_Transformer()
# transformed <- lt$log2()$robscale()$lfq
transformed <- lt$robscale()$lfq
transformed$config$table$is_response_transformed

pl <- transformed$get_Plotter()
pl$intensity_distribution_density() + 
  xlim(c(10,23))

normalized_proteomics = transformed$data
saveRDS(normalized_proteomics, 'data/proteomics_normalized.rds')

prot_data = lfqdata$data %>% 
  mutate(log_intensity = log(Intensity+.1))


x = transformed$data %>% filter(PG.Genes == 'KRAS') 


hist((x$Intensity - mean(x$Intensity))/sd(x$Intensity))

genes_cna_status = readRDS('data/karyotypes_mutations_all_genes_qc_ccf_v4.rds')

# get correct sample names
dict = readRDS('data/full_dict_dna_rna_prot.rds')

# get the pair proteomics-genomics correct sample name
samples_check = dict %>% 
  dplyr::select(proteomics_code, DNA, fixed_name) %>%
  dplyr::distinct()

genes_to_check = genes_cna_status$hgnc_symbol %>% unique



genes_cna_status = genes_cna_status %>% 
  dplyr::mutate(mut_consequence = ifelse(is.na(mut_consequence), 'wild-type', mut_consequence))

kras_genes_cna_status = genes_cna_status %>% 
  filter(hgnc_symbol == 'KRAS')

kras_genes_cna_status = full_join(dict, kras_genes_cna_status, by = join_by('fixed_name' == 'sample')) %>% 
  mutate(proteomics_code = gsub('_', '', proteomics_code))

x = x %>% 
  mutate(scaled_int = (Intensity - mean(Intensity))/sd(Intensity))

x = x %>% 
  full_join(., kras_genes_cna_status, by = join_by('PDO' == 'proteomics_code'))

x %>% 
  select(PG.Genes, karyotype, PDO, scaled_int) %>% 
  distinct() %>% 
  # separate(karyotype, into = c('Major', 'minor'), sep = ':', convert = T) %>% 
  # mutate(tot_cna = Major+minor) %>% 
  ggplot(aes(x = karyotype, y = scaled_int, colour = PDO)) + 
  # geom_point() +
  geom_boxplot()
  
