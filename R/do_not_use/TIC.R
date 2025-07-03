# new normalization method
library(limma)
library(tidyverse)
library(DEP)

data = readxl::read_excel('data/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
  as.data.frame()

data[which(is.na(data[,'PG.Genes'])), 'PG.Genes'] = 'Unknown'
raw_data = data %>% 
  dplyr::select(everything(), -PG.ProteinGroups) %>% 
  tibble::column_to_rownames('PG.Genes') %>% 
  # dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric())) %>% 
  dplyr::mutate(across(where(is.character), as.numeric)) %>% 
  dplyr::mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

# log_transformed = log(raw_data, base = 2)

# try to compute tic
# Compute TIC per sample
tic_factors <- colSums(raw_data, na.rm = TRUE)

# Normalize using TIC
normalized_data <- sweep(raw_data, 2, tic_factors, FUN = "/") * median(tic_factors)

proteomics_data_norm = normalized_data %>% 
  tibble::rownames_to_column('Genes') %>% 
  reshape2::melt() %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(value = log(value, base = 2)) 
# checking and imputing NAs -- using DEP 

# # creating experimental table design
# exp_design = tibble(
#   label = colnames(data)[-c(1,2)]
# ) %>% 
#   mutate(condition = gsub('_a|_b', '', label)) %>% 
#   mutate(replicate = str_extract(label, 'a$|b$')) 




# normalized_data = normalized_data %>% 
#   tibble::rownames_to_column('genes') %>% 
#   full_join(., data[,c(1, 2)], by = join_by('genes' == 'PG.Genes'))

# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
# data_unique <- make_unique(normalized_data, "Genes", "PG.ProteinGroups", delim = ";")
# 
# # columns
# samples_columns = which(!colnames(normalized_data) %in% c('genes', 'PG.ProteinGroups'))
# 
# dep = make_se(data_unique, samples_columns, exp_design)
# 
# proteomics_data_norm = dep %>% DEP::get_df_long()


# map the sample names
cancer_genes  = read.table('~/Google Drive/Il mio Drive/PhD/organoids/Cancer_gene_census_tier1.csv', sep = ",", header = T) %>% 
  dplyr::as_tibble() %>% 
  # dplyr::filter(Germline != 'yes') %>% 
  dplyr::filter(grepl('*colorectal*', Tumour.Types.Somatic.)) %>% 
  dplyr::select('Gene.Symbol', 'Role.in.Cancer')

new_mapping_experiment = readRDS("~/Google Drive/Il mio Drive/PhD/organoids/mapping_samples.rds")
# drivers = readRDS('drivers_info/drivers.rds')
# genes_to_check = c("APC", "KRAS", "TP53", "PIK3CA")
genes_cna_status = readRDS('data/karyotypes_full_cohort_cgs.rds')
genes_to_check = cancer_genes$Gene.Symbol %>% unique




drivers_to_check_correct_samples = genes_cna_status %>% 
  # dplyr::mutate(karyotype = paste(Major, minor, sep = ":")) %>% 
  dplyr::filter(sample %in% (new_mapping_experiment$fixed_name %>% unique)) %>%
  dplyr::filter(region %in% genes_to_check) #%>% 
  # dplyr::full_join(., genes_phasing, by = join_by('sample' == 'sample', 'driver' == 'VEP.SYMBOL')) %>% # aggiungi info su mutazione e phasing
  # dplyr::mutate(multiplicity = ifelse(is.na(multiplicity), 0, multiplicity)) %>% 
  # dplyr::rename(Consequence = VEP.Consequence) %>% 
  # dplyr::mutate(Consequence = ifelse(multiplicity == 0, "not mutated", Consequence))

# drivers_not_mutated = readRDS('~/Google Drive/Il mio Drive/PhD/organoids/genes_phasing.rds')
drivers_not_mutated = readRDS('~/Google Drive/Il mio Drive/PhD/organoids/cgs_pos_cna_state.rds')
# drivers_not_mutated = lapply(names(drivers_not_mutated), function(x) {
#   drivers_not_mutated[[x]] %>%
#     dplyr::mutate(sample = x) %>%
#     dplyr::select(c(chr, from, to, ref, alt, DP, NV, VAF, VEP.SYMBOL, driver_label, sample))
# }) %>%
#   dplyr::bind_rows()

# drivers_not_mutated = read.table("~/Google Drive/Il mio Drive/PhD/organoids/cnas_important_genes.csv", sep = ",", header = T)
# # 
# drivers_not_mutated = drivers_not_mutated %>%
#   dplyr::mutate(gene = case_when(chr == "chr17" ~ "TP53",
#                                  chr == "chr5" ~ "APC",
#                                  chr == "chr3" ~ "PIK3CA",
#                                  chr == "chr12" ~ "KRAS")) %>%
#   dplyr::mutate(from = NULL, X = NULL, to = NULL)
# 
drivers_to_check_correct_samples = full_join(drivers_to_check_correct_samples, drivers_not_mutated, by = join_by("chr" == "chr",
                                                                                                                 "gene" == "gene",
                                                                                                                 "PDO" == "sample",
                                                                                                                 "Major" == "Major",
                                                                                                                 "minor" == "minor"))

drivers_to_check_correct_samples = drivers_to_check_correct_samples %>%
  dplyr::mutate(karyotype = ifelse(is.na(karyotype), paste(Major, minor, sep = ":"), karyotype)) %>%
  dplyr::mutate(multiplicity = ifelse(is.na(multiplicity), 0, multiplicity)) %>%
  dplyr::mutate(Consequence = ifelse(multiplicity == 0, "not mutated", Consequence))


# Filter the data to keep only selected genes
proteomic_to_check = proteomics_data_norm %>%
  # tidyr::separate(gene, into = c("uniprot", "gene_name"), sep = "_") %>%
  # dplyr::filter(gene_name %in% (cancer_genes_somatic_colon$driver %>% unique))
  dplyr::filter(Genes %in% genes_to_check) 

# proteomic_genes = proteomic_to_check %>% 
#   dplyr::select(everything(), -uniprot) %>%
#   reshape2::melt()

proteomic_genes = proteomic_to_check %>% 
  # dplyr::rename(sample_id = variable) %>% 
  dplyr::mutate(PDO = gsub("_a|_b", "", variable)) %>% 
  dplyr::group_by(PDO, Genes) %>% 
  dplyr::mutate(PDO = gsub('_HSR', 'HSR', PDO))
# dplyr::summarise(mean_expr = mean(value), sd_expr = sd(value)) 

# proteomic_genes %>% 
#   filter(gene_name == "GAPDH") %>%
#   ggplot(aes(value)) + 
#   geom_histogram() + 
#   facet_wrap(vars(gene_name))

samples_check = new_mapping_experiment %>% 
  dplyr::select(proteomics_code, fixed_name) %>%
  dplyr::distinct()

x = drivers_to_check_correct_samples %>%
  dplyr::select(gene, karyotype, PDO, multiplicity, Consequence, Major, minor) %>% 
  # mutate(state = paste0(karyotype, "-", multiplicity)) %>% 
  distinct(.keep_all = F) %>% 
  left_join(., samples_check, by = join_by("PDO" == "fixed_name")) %>% 
  # filter(!is.na(gene)) %>% 
  full_join(., proteomic_genes, by = join_by("proteomics_code" == "PDO", "gene" == "Genes")) %>% 
  # mutate(state = ifelse(is.na(state), "no alteration known", state)) %>% 
  # mutate(multiplicity = ifelse(is.na(multiplicity), 0, multiplicity)) %>% 
  # mutate(karyotype = ifelse(is.na(karyotype), "no alteration", karyotype))
  dplyr::filter(!is.na(proteomics_code)) 

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
                                 .default = Consequence)) %>% 
  dplyr::mutate(tot_cna = Major+minor)

# colors = setNames(Polychrome::kelly.colors(length(x$karyotype %>% unique)), nm = x$karyotype %>% unique)
colors = setNames(ggthemes::calc_pal()(length(x$karyotype %>% unique)), nm = x$karyotype %>% unique)

x = x %>% 
  dplyr::mutate(karyotype = factor(karyotype, 
                                   levels = c("1:0", "1:1", "2:0", "2:1", "3:0", "2:2", "3:1", "4:0", "3:2", "5:0", "3:3")))

x %>% 
  filter(gene %in% genes_to_check) %>% 
  # full_join(., cancer_genes, by = join_by('gene' == 'Gene.Symbol')) %>% 
  # dplyr::rename(Role = Role.in.Cancer) %>% 
  # filter(driver %in% c("APC", "KRAS", "TP53", "PIK3CA")) %>% 
  # filter(gene=='KRAS') %>%
  # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
  # mutate(expr_by_cna = value/tot_cna) %>% 
  ggplot(aes(x = tot_cna, y = value, colour = karyotype)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~gene, scales = 'free') +
  xlab('total cnas') + 
  scale_color_manual(values = colors) +
  ggtitle('Normalized data, NAs not imputed')
  # geom_smooth(method=lm , se=TRUE)

x %>% 
  filter(gene %in% genes_to_check) %>% 
  # full_join(., cancer_genes, by = join_by('gene' == 'Gene.Symbol')) %>% 
  # dplyr::rename(Role = Role.in.Cancer) %>% 
  # filter(driver %in% c("APC", "KRAS", "TP53", "PIK3CA")) %>% 
  # filter(gene=='KRAS') %>%
  # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
  mutate(expr_by_cna = value/tot_cna) %>% 
  ggplot(aes(x = tot_cna, y = expr_by_cna, colour = karyotype)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~gene, scales = 'free') +
  xlab('protein_expr/tot_cna') + 
  scale_color_manual(values = colors) + 
  ggtitle('Normalized data, NAs not imputed')

x %>% 
  filter(gene %in% genes_to_check) %>% 
  # full_join(., cancer_genes, by = join_by('gene' == 'Gene.Symbol')) %>% 
  # dplyr::rename(Role = Role.in.Cancer) %>% 
  # filter(driver %in% c("APC", "KRAS", "TP53", "PIK3CA")) %>% 
  # filter(gene=='KRAS') %>%
  # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
  mutate(expr_by_cna = value/tot_cna) %>%
  ggplot(aes(x = expr_by_cna, y = tot_cna, colour = karyotype)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~gene, scales = 'free') +
  xlab('protein_expr/tot_cna') + 
  scale_color_manual(values = colors) + 
  ggtitle('Normalized data, NAs not imputed')


x %>% 
  # filter(gene=='KRAS') %>%
  filter(!is.na(karyotype)) %>% 
  ggplot(aes(intensity, fill = karyotype)) +
  geom_histogram(binwidth = 0.1) +
  # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
  # mutate(expr_by_cna = value/tot_cna) %>% 
  # ggplot(aes(x = expr_by_cna, y = tot_cna, colour = karyotype)) +
  # geom_point() +
  theme_bw() +
  facet_wrap(karyotype~gene, scales = 'free') +
  # scale_color_brewer(palette = 'Paired')
  # ggsci::scale_color_futurama()
  scale_fill_manual(values = colors) + 
  ggtitle('Normalized data, NAs imputed')

