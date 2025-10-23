library(tidyverse)
library(prolfqua)

source('organoids_analysis/R/functions_utils/proteomics_utils.R')
# source('organoids_analysis/R/functions_utils/fit_plots.R')
# source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

# get correct sample names
new_mapping_experiment = readRDS("data/mapping_samples.rds")

# get the pair proteomics-genomics correct sample name
samples_check = new_mapping_experiment %>% 
  dplyr::select(proteomics_code, fixed_name) %>%
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
  dplyr::full_join(., ann, by = join_by('sample' == 'replicate'))

# creating AnalysisTableAnnotation obj
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "sample"
atable$workIntensity = "Intensity"
atable$hierarchy[['PG.Genes']] = 'PG.Genes'
atable$factors[['PDO']] = 'PDO'

config <- prolfqua::AnalysisConfiguration$new(atable)
analysis_data <- prolfqua::setup_analysis(data, config)
lfqdata <- prolfqua::LFQData$new(analysis_data, config)

# visualize the intensities
smrz <- lfqdata$get_Summariser()
smrz$plot_hierarchy_counts_sample()

lfqplotter <- lfqdata$get_Plotter()
density_nn <- lfqplotter$intensity_distribution_density()
density_nn + 
  theme(legend.position = 'bottom')

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
transformed <- lt$log2()$robscale()$lfq
transformed$config$table$is_response_transformed

pl <- transformed$get_Plotter()
pl$intensity_distribution_density() + 
  xlim(c(10,23))

normalized_proteomics = transformed$data
saveRDS(normalized_proteomics, 'data/proteomics_normalized.rds')


