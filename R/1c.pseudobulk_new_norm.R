library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Matrix.utils)
library(DESeq2)
# library(harmony)
# library(RUVSeq)
setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
source('organoids_analysis/R/functions_utils/scRNA_utils.R')
# source('../TLS/CLL_Cu/copper/CuRes/R/cu_vs_atn/remove_batch_effect.R')

data_path = 'data/scRNA'
data = lapply(list.files(data_path, full.names = T), function(x) {
  readRDS(x)
})

samples = list.files(data_path) %>% 
  gsub('_filtered.rds', '', .) %>% 
  gsub('Sample_', '', .)

names(data) = samples

data_v2 = lapply(data, function(x) {
  x$batch = x@meta.data %>% 
    dplyr::mutate(batch = case_when(
      str_detect(sample, 'PDO') ~ 'HSR', 
      .default = 'ICR')
    ) %>% 
    dplyr::pull(batch)  
  
  return(x)
})

# generate the pseudobulk
sces = lapply(data_v2, create_sc_expr)

aggregated_data = lapply(sces, aggregate_cells)

pseudobulk = generate_pseudobulk_count_matrix(aggregated_data) 

pseudobulk = pseudobulk %>% 
  dplyr::filter(if_all(everything(), ~ !is.na(.))) %>% 
  filter(!grepl("^MT-", genes)) %>% 
  tibble::column_to_rownames('genes')

# now treat it as common deseq analysis with batch effect correction --- test using RUVseq
coldata = tibble(
  organoid = colnames(pseudobulk)
) %>%
  as.data.frame()
rownames(coldata) = coldata$organoid

coldata = lapply(data_v2, function(x) {
  x@meta.data$batch %>% unique
}) %>% 
  as_tibble %>% 
  t %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('organoid') %>% 
  full_join(coldata, .) %>% 
  dplyr::rename(batch = 'V1')
rownames(coldata) = coldata$organoid  

filtered = apply(pseudobulk, 1, function(x) length(x[x > 5]) >= 4)
filtered_counts = pseudobulk[filtered,]

coldata = coldata %>% 
  mutate(batch = factor(batch, levels = c('ICR', 'HSR')))

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData   = coldata,
  design    = ~ batch
)

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
norm_rna <- assay(vsd)

dict = readRDS('data/full_dict_dna_rna_prot.rds')
genes_karyo_muts = readRDS('data/karyotypes_mutations_all_genes_qc_ccf_v4.rds')

samples <- colnames(norm_rna)

norm_rna_long = norm_rna %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('gene') %>% 
  pivot_longer(cols = all_of(samples), names_to = 'sample', values_to = 'expression')

norm_rna_long = norm_rna_long %>% 
  full_join(., dict, by = join_by('sample' == 'RNA')) %>% 
  filter(!is.na(expression)) 

# not accounting yet for multihits and multi-segments
norm_rna_long_v2 = norm_rna_long %>% 
  # filter(!is.na(expression)) %>% 
  left_join(., genes_karyo_muts, by = join_by('fixed_name' == 'sample', 
                                     'gene' == 'hgnc_symbol')) %>% 
  filter(!grepl('^LINC', gene))

saveRDS(norm_rna_long_v2, 'data/rna_new_norm_merged.rds')


# 
# p1 = norm_rna_long_v2 %>% 
#   filter(gene == 'KRAS') %>% 
#   group_by(sample) %>% 
#   filter(n() > 1) %>% 
#   tidyr::separate(karyotype, sep = ':', into = c('Major', 'minor'), convert = T) %>% 
#   mutate(tot_cna = Major+minor) %>% 
#   ggplot(aes(x = tot_cna, y = expression)) + 
#   geom_point() + 
#   theme_bw()
# 
# 
# p2 = x %>% 
#   filter(hgnc_symbol == 'KRAS') %>% 
#   # tidyr::separate(karyotype, sep = ':', into = c('Major', 'minor'), convert = T) %>% 
#   # mutate(tot_cna = Major+minor) %>% 
#   ggplot(aes(x = tot_cna, y = value)) + 
#   geom_point() + 
#   theme_bw()
# 
# p1/p2

# expr_set = newSeqExpressionSet(as.matrix(filtered_counts), 
#                                phenoData = coldata)
# 
# expr_set <- betweenLaneNormalization(expr_set, which = "upper")
# 
# design <- model.matrix(~batch, data=pData(expr_set))
# y <- DGEList(counts=counts(expr_set))
# y <- calcNormFactors(y, method="upperquartile")
# y <- estimateGLMCommonDisp(y, design)
# y <- estimateGLMTagwiseDisp(y, design)
# 
# fit <- glmFit(y, design)
# res <- residuals(fit, type="deviance")
# 
# ruvRes = RUVr(expr_set, cIdx = rownames(expr_set), k = 1, res)
# 
# dds <- DESeqDataSetFromMatrix(countData = counts(expr_set),
#                               colData = pData(ruvRes),
#                               design = ~ batch + W_1)
# dds <- DESeq(dds)
# 
# # get normalized vst data
# # normalized_res = assay(vst(dds, blind = FALSE))
# normalized_res = counts(dds,normalized=TRUE)
# saveRDS(normalized_res, 'data/normalized_res_pseudobulk_v2.rds')
saveRDS(normalized_res, 'data/normalized_res_pseudobulk_counts.rds')

new_mapping_experiment = readRDS("data/mapping_samples.rds")
all_genes = rownames(normalized_res) %>% unique

saveRDS(all_genes, 'data/genes_to_check.rds')




