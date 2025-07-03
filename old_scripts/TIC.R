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

# try to compute tic
# Compute TIC per sample
tic_factors <- colSums(raw_data, na.rm = TRUE)

apply(raw_data, 2, function(x) median(x, na.rm = TRUE)) %>% median

# Normalize using TIC
normalized_data <- sweep(raw_data, 2, tic_factors, FUN = "/") * mean(tic_factors)

# checking and imputing NAs -- using DEP 

# creating experimental table design
exp_design = tibble(
  label = colnames(data)[-c(1,2)]
) %>% 
  mutate(condition = gsub('_a|_b', '', label)) %>% 
  mutate(replicate = str_extract(label, 'a$|b$')) 

normalized_data = normalized_data %>% 
  tibble::rownames_to_column('proteins') %>% 
  full_join(., (data %>% dplyr::select(PG.ProteinGroups, PG.Genes)), by = join_by('proteins' == 'PG.Genes'))

unique_data = make_unique(normalized_data, names = 'proteins', ids = 'PG.ProteinGroups', delim = '_')
samples_index = which(!colnames(unique_data) %in% c('PG.ProteinGroups', 'proteins', 'name', 'ID'))

dep = DEP::make_se(proteins_unique = unique_data, columns = samples_index, expdesign = exp_design)

# rownames(unique_data) = unique_data$proteins
# data = unique_data[,samples_index]
# colnames(data) = colnames(assay(dep))
# rownames(data) = rownames(dep)
# assay(dep) = data

data_imp <- impute(dep, fun = "knn")

data_reconstruct = DEP::get_df_long(data_imp)

drivers = readRDS("drivers.rds")

# map sample names among the experiments (hopefully)
mapping = readxl::read_excel("mapping_ICR_names.xlsx")

mapping = mapping %>% 
  select(`Patient code (HSR)`, `Patient code (ICR)`, fixed_name) %>% 
  dplyr::rename(proteomics_code = `Patient code (HSR)`) %>% 
  dplyr::rename(genomics_code = `Patient code (ICR)`) %>% 
  dplyr::mutate(proteomics_code = gsub("#", "", proteomics_code))

mapping_bis = data_reconstruct %>% 
  dplyr::select(label) %>% 
  dplyr::mutate(genomics_code = gsub('_a$|_b$','', label)) %>% 
  dplyr::rename(sample_proteomics_code = label) %>%
  distinct()

tt = full_join(mapping, mapping_bis, by = join_by("proteomics_code" == "genomics_code"))

new_mapping = tt %>% 
  mutate(fixed_name = ifelse(is.na(fixed_name), gsub("HSR", "", proteomics_code), fixed_name)) %>% 
  mutate(fixed_name = gsub("-", "_", fixed_name)) 

new_mapping_experiment = new_mapping %>% 
  dplyr::filter(!is.na(sample_proteomics_code)) %>% 
  dplyr::mutate(fixed_name = gsub("_seg", "", fixed_name)) %>% 
  dplyr::mutate(fixed_name = case_when(fixed_name == "65_VI" ~ "65_IV", 
                                       .default = fixed_name))

# select the top genes to check first
genes_to_check = c("APC", "KRAS", "TP53", "PIK3CA")

drivers_to_check = drivers %>% 
  dplyr::filter(gene %in% genes_to_check)

drivers_to_check_correct_samples = drivers_to_check %>% 
  filter(PDO %in% new_mapping_experiment$fixed_name %>% unique)

drivers_not_mutated = read.table("cnas_important_genes.csv", sep = ",", header = T)
drivers_not_mutated = drivers_not_mutated %>% 
  dplyr::mutate(gene = case_when(chr == "chr17" ~ "TP53", 
                                 chr == "chr5" ~ "APC",
                                 chr == "chr3" ~ "PIK3CA", 
                                 chr == "chr12" ~ "KRAS")) %>% 
  dplyr::mutate(from = NULL, X = NULL, to = NULL)

drivers_to_check_correct_samples = full_join(drivers_to_check_correct_samples, drivers_not_mutated, by = join_by("chr" == "chr", 
                                                                                                                 "gene" == "gene", 
                                                                                                                 "PDO" == "sample", 
                                                                                                                 "Major" == "Major", 
                                                                                                                 "minor" == "minor"))

drivers_to_check_correct_samples = drivers_to_check_correct_samples %>% 
  dplyr::mutate(karyotype = ifelse(is.na(karyotype), paste(Major, minor, sep = ":"), karyotype)) %>% 
  dplyr::mutate(multiplicity = ifelse(is.na(multiplicity), 0, multiplicity)) %>% 
  dplyr::mutate(Consequence = ifelse(multiplicity == 0, "not mutated", Consequence))

proteomics = data_reconstruct %>% 
  filter(name %in% genes_to_check) %>% 
  select(label, name, intensity) %>% 
  mutate(PDO = gsub("_a|_b", "", label)) %>% 
  dplyr::group_by(PDO, name)

# dplyr::summarise(mean_expr = mean(value), sd_expr = sd(value)) 

samples_check = new_mapping_experiment %>% 
  select(proteomics_code, fixed_name) %>%
  distinct

x = drivers_to_check_correct_samples %>%
  select(gene, karyotype, PDO, multiplicity, Consequence, Major, minor) %>% 
  # mutate(state = paste0(karyotype, "-", multiplicity)) %>% 
  distinct(.keep_all = F) %>% 
  left_join(., samples_check, by = join_by("PDO" == "fixed_name")) %>% 
  # filter(!is.na(gene)) %>% 
  full_join(., proteomics, by = join_by("proteomics_code" == "PDO", "gene" == "name")) %>% 
  # mutate(state = ifelse(is.na(state), "no alteration known", state)) %>% 
  # mutate(multiplicity = ifelse(is.na(multiplicity), 0, multiplicity)) %>% 
  # mutate(karyotype = ifelse(is.na(karyotype), "no alteration", karyotype))
  dplyr::filter(!is.na(proteomics_code))


# first visualization
x = x %>% 
  mutate(karyotype = factor(karyotype,
                            levels = c("1:0", "1:1", "2:0", "2:1", "3:0", "2:2", "3:1", "4:0","3:2", "5:0", "3:3"))) %>% 
  mutate(multiplicity = factor(multiplicity, levels = 0:5))

multihit_identifiers = x %>% 
  group_by(gene, karyotype, PDO, multiplicity) %>% 
  dplyr::count()

x = full_join(x, multihit_identifiers, by = join_by("gene" == "gene", "karyotype" == 'karyotype', "PDO" == 'PDO', 'multiplicity' == 'multiplicity'))
x = x %>% 
  dplyr::mutate(Consequence = ifelse(n > 1, "multihit", Consequence)) %>% 
  arrange(desc(n)) %>% 
  distinct() %>% 
  separate(Consequence, into = "Consequence", sep = "&") %>% 
  mutate(Consequence = case_when(Consequence %in% c("frameshift_variant","inframe_deletion") ~ "frame_alteration", 
                                 .default = Consequence))


x = x %>% 
  dplyr::mutate(tot_cna = Major+minor)

colors = setNames(ggthemes::calc_pal()(11), nm = x$karyotype %>% unique)

x %>% 
  # filter(gene=='KRAS') %>%
  # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
  mutate(expr_by_cna = intensity/tot_cna) %>% 
  ggplot(aes(x = expr_by_cna, y = tot_cna, colour = karyotype)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~gene, scales = 'free') +
  # scale_color_brewer(palette = 'Paired')
  # ggsci::scale_color_futurama()
  scale_color_manual(values = colors)

