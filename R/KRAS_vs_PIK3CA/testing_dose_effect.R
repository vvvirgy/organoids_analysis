.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4')
rm(list=ls())
library(CNAqc)
library(tidyverse)
library(gprofiler2)
library(enrichR)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
cnas = readRDS('data/cnaqc/cnas_list_rough_ccf.rds')

# using only KRAS mutated samples with associated CCF
samples_to_use = lapply(cnas, function(x) {
  CCF(x) %>% 
    dplyr::filter(VEP.SYMBOL == 'KRAS') %>% 
    dplyr::filter(is_driver)
})
samples_to_use = samples_to_use[which(sapply(samples_to_use, nrow) > 0)] %>% names

genes_list = lapply(cnas[samples_to_use], function(x) {
  CCF(x) %>% 
    dplyr::filter(VAF > 0) %>% 
    filter(VEP.SYMBOL != '.') %>% 
    dplyr::select(VEP.SYMBOL, sample) %>% 
    distinct()
}) %>% 
  bind_rows() %>% 
  # pull(VEP.SYMBOL) %>% 
  # grep('^PIK',., value = T) %>% 
  # unique
  group_by(VEP.SYMBOL) %>%
  count %>%
  dplyr::filter(n >= 5) %>%
  dplyr::filter(VEP.SYMBOL != 'KRAS')
genes_list = genes_list$VEP.SYMBOL

get_gene_multiplicity = function(x, genes) {
  CNAqc::CCF(x) %>% 
    dplyr::filter(VEP.SYMBOL %in% genes) %>% 
    # dplyr::filter(is_driver) %>% 
    dplyr::select(chr, from, to, ref, alt, VEP.SYMBOL, VAF, mutation_multiplicity, CCF, karyotype, sample, is_driver)
}

# re-do everything using all genes!!!
muts_of_interest = lapply(cnas[samples_to_use], function(x) {
  get_gene_multiplicity(x, genes = c('KRAS', genes_list))
}) %>% 
  bind_rows() 
muts_of_interest = muts_of_interest %>% 
  mutate(is_driver = ifelse(VEP.SYMBOL != 'KRAS', TRUE, is_driver))

drivers = muts_of_interest %>% 
  filter(is_driver)

# run for each copy
genes_list = gsub('-', '_', genes_list)
drivers = drivers %>% 
  dplyr::mutate(VEP.SYMBOL = gsub('-', '_', VEP.SYMBOL))

# prepare data to run Rediscover
# for KRAS --> m = 1 is 0, m > 1 is 1
# for the other genes --> wild-type is 0, mutated is 1

tests_all_genes = lapply(genes_list, function(x) {
  df = drivers %>% 
    dplyr::filter(VEP.SYMBOL %in% c(x, 'KRAS')) %>% 
    dplyr::select(VEP.SYMBOL, mutation_multiplicity, sample) %>% 
    group_by(VEP.SYMBOL, sample) %>% 
    mutate(mutation_multiplicity = case_when(
      (VEP.SYMBOL == 'KRAS' & mutation_multiplicity > 1) ~ 'm > 1', 
      (VEP.SYMBOL == 'KRAS' & mutation_multiplicity == 1) ~ 'm = 1', 
      (VEP.SYMBOL != 'KRAS' & !is.na(mutation_multiplicity)) ~ 'mutated', 
    )) %>% 
    distinct() %>% 
    tidyr::pivot_wider(names_from = VEP.SYMBOL, values_from = mutation_multiplicity, values_fill = 'wild-type') %>% 
    mutate(!!x := factor(.data[[x]], levels = c('wild-type', 'mutated'))) %>% 
    mutate(KRAS = factor(KRAS, levels = c('m = 1', 'm > 1'))) %>% 
    mutate(!!x := ifelse(.data[[x]] == 'wild-type', 0, 1)) %>% 
    mutate(KRAS = ifelse(KRAS == 'm = 1', 0, 1)) 
  return(df)
})
names(tests_all_genes) = genes_list
saveRDS(tests_all_genes, 'data/pik3_complexes_mut_status_kras_multiplicity.rds')

tt = readRDS('data/pik3_complexes_mut_status_kras_multiplicity.rds')
# fisher test
tests_all_genes = lapply(genes_list, function(x) {
  form = as.formula(paste0('n ~ KRAS + ', x))
  print(form)
  df = drivers %>% 
    dplyr::filter(VEP.SYMBOL %in% c(x, 'KRAS')) %>% 
    dplyr::select(VEP.SYMBOL, mutation_multiplicity, sample) %>% 
    group_by(VEP.SYMBOL, sample) %>% 
    mutate(mutation_multiplicity = case_when(
      (VEP.SYMBOL == 'KRAS' & mutation_multiplicity > 1) ~ 'm > 1', 
      (VEP.SYMBOL == 'KRAS' & mutation_multiplicity == 1) ~ 'm = 1', 
      (VEP.SYMBOL != 'KRAS' & !is.na(mutation_multiplicity)) ~ 'mutated', 
    )) %>% 
    distinct() %>% 
    tidyr::pivot_wider(names_from = VEP.SYMBOL, values_from = mutation_multiplicity, values_fill = 'wild-type') %>% 
    group_by(pick(c(x, 'KRAS'))) %>% 
    count() %>% 
    # dplyr:: mutate(across(all_of(x), ~ factor(.x, levels = c('mutated', 'wild-type'))))
    mutate(!!x := factor(.data[[x]], levels = c('mutated', 'wild-type'))) %>% 
    mutate(KRAS = factor(KRAS, levels = c('m = 1', 'm > 1'))) %>% 
    xtabs(form, data = ., drop.unused.levels = F) %>% 
    fisher.test()
})
  
names(tests_all_genes) = genes_list
saveRDS(tests_all_genes, 'data/tests_all_genes_with_kras_multiplicity.rds')

sign_associations = lapply(tests_all_genes %>% names, function(x){
  tests_all_genes[[x]] %>% 
    broom::tidy() %>% 
    mutate(gene = x)
}) %>% 
  bind_rows() %>% 
  filter(p.value <= 0.05)

sign_associations %>% view()
write.table(sign_associations, 'res/tests_all_genes_with_kras_multiplicity.csv', sep = ',', quote = F, col.names = T, row.names = F)

genes_sign_associations = sign_associations$gene %>% unique

pathways = gprofiler2::gost(genes_sign_associations, 
                            organism = 'hsapiens', 
                            ordered_query = F,
                            multi_query = F,
                            significant = F,
                            exclude_iea = T, 
                            measure_underrepresentation = T, 
                            evcodes = TRUE, 
                            correction_method = 'fdr', 
                            domain_scope = 'annotated', 
                            sources = c('GO', 'KEGG', 'REAC', 'WP', 'CORUM', 'HP', 'HPA'),
                            highlight = FALSE
                            ) 
pathways$result %>% filter(significant) %>% dim


# trying enrichr since gprofiler2 did not give any significant result

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

if (websiteLive) {
  dbs <- listEnrichrDbs()
  head(dbs)
}

dbs_small = c('GO_Biological_Process_2025', 'GO_Cellular_Component_2025', 'GO_Molecular_Function_2025')
enriched_res <- enrichr(genes_sign_associations, dbs_small)

enriched_res_all = enriched_res %>% 
  bind_rows()

enriched_res_all %>% filter(P.value <= 0.05) %>% view

plotEnrich(enriched_res_all, showTerms = 40, numChar = 50, 
           y = "Count", orderBy = "P.value") 

enriched_res_all %>% 
  filter(P.value <= 0.05) %>% 
  separate(Overlap, into = c('n_genes', 'set_dim'), sep = '/', convert = T) %>% 
  mutate(ratio = n_genes/set_dim) %>% 
  arrange(desc(ratio)) %>% head(30)

# tomorrow --> run an enrichement analysis to see which are the functions that are altered
###############################################################################################################
# 
# drivers_muts = drivers %>%
#   dplyr::select(VEP.SYMBOL, mutation_multiplicity, sample) %>%
#   tidyr::pivot_wider(names_from = VEP.SYMBOL, values_from = mutation_multiplicity) %>%
#   filter(!is.na(KRAS)) %>%
#   mutate(KRAS = ifelse(KRAS > 1, 'm > 1', 'm = 1')) %>%
#   mutate(PIK3CA = ifelse(is.na(PIK3CA), 'wild-type', 'mutated'))
# 
# test = drivers_muts %>%
#   group_by(KRAS, PIK3CA) %>%
#   count() %>%
#   xtabs(n ~ KRAS + PIK3CA, data = .) %>%
#   # as.data.frame() %>%
#   # mutate(Freq = Freq + 0.5) %>%
#   # xtabs(Freq ~ KRAS + PIK3CA, data = .) %>%
#   fisher.test()
# broom::tidy(test)
# 
# ggstatsplot::ggbarstats(data = drivers_muts, PIK3CA, KRAS)
# ggsave('res/kras_pik3ca_association.pdf', width = 6, height = 6)
# ,
#                         results.subtitle = FALSE,
#                         title = paste(x$info$region %>% unique, x$info$CNA %>% unique),
#                         subtitle = paste0(
#                           "Fisher's exact test", ", p-value = ",
#                           ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
#                         ),
#                         package = "RColorBrewer",
#                         palette = 'Dark2')


# kras = lapply(cnas, function(x) {
#   Mutations(x) %>% 
#     filter(VEP.SYMBOL == 'KRAS') %>% 
#     select(chr, from, to, ref, alt, VEP.SYMBOL,VAF, karyotype, sample, is_driver) %>% 
#     filter(is_driver == TRUE)
# }) 
# 
# setdiff( ( muts_of_interest %>% filter(VEP.SYMBOL == 'KRAS') %>% filter(is_driver) %>% pull(sample)),  (which(sapply(kras, nrow) > 0) %>% names))
