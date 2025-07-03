library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Matrix.utils)
library(DESeq2)
# library(harmony)
library(RUVSeq)
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

# generate the pseudobulk
sces = lapply(data, create_sc_expr)

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

filtered = apply(pseudobulk, 1, function(x) length(x[x > 5]) >= 2)
filtered_counts = pseudobulk[filtered,]
expr_set = newSeqExpressionSet(as.matrix(filtered_counts), 
                               phenoData = coldata)

expr_set <- betweenLaneNormalization(expr_set, which = "upper")

design <- model.matrix(~1, data=pData(expr_set))
y <- DGEList(counts=counts(expr_set))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

ruvRes = RUVr(expr_set, cIdx = rownames(expr_set), k = 1, res)

dds <- DESeqDataSetFromMatrix(countData = counts(expr_set),
                              colData = pData(ruvRes),
                              design = ~ W_1)
dds <- DESeq(dds)

# get normalized vst data
normalized_res = assay(vst(dds, blind = FALSE))
saveRDS(normalized_res, 'data/normalized_res_pseudobulk.rds')

new_mapping_experiment = readRDS("data/mapping_samples.rds")
all_genes = rownames(normalized_res) %>% unique

saveRDS(all_genes, 'data/genes_to_check.rds')

# # testing normalization and using harmony to remove batch effect
# total_reads<- colSums(pseudobulk)
# final_mat<- t(t(pseudobulk)/total_reads)
# final_mat<- log2(final_mat + 1)
# 
# library(genefilter)
# 
# # choose the top 1000 most variabel genes 
# top_genes<- genefilter::rowVars(final_mat) %>% 
#   sort(decreasing = TRUE) %>%
#   names() %>%
#   head(1000)
# 
# # subset only the top 1000 genes
# expression_mat_sub<- final_mat[top_genes, ]
# 
# # calculate the PCA
# pca<- prcomp(t(expression_mat_sub),center = TRUE, scale. = TRUE) 
# 
# PC1_and_PC2<- data.frame(PC1=pca$x[,1], PC2= pca$x[,2])
# PC1_and_PC2 = PC1_and_PC2 %>% 
#   tibble::rownames_to_column('PDO')
# 
# 
# p1<- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
#   geom_point(aes(color = PDO)) +
#   theme_bw(base_size = 14) 
# 
# set.seed(123)
# 
# harmony_embeddings <- harmony::HarmonyMatrix(
#   expression_mat_sub, 
#   meta_data = coldata, 
#   vars_use = 'organoid'
# )
# 
# rownames(harmony_embeddings)<- rownames(coldata)
# harmony_pc<- data.frame(harmony1=harmony_embeddings[,1], harmony2= harmony_embeddings[,2])
# 
# harmony_pc<- cbind(harmony_pc, final_meta)



##########################################

# normalized_res = counts(dds, normalized = TRUE)

# saveRDS(normalized_res, 'data/normalized_res_pseudobulk.rds')




