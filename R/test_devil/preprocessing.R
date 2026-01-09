rm(list=ls())

library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
# library(DESeq2)

data_path = 'data/scRNA'
# data = lapply(list.files(data_path, full.names = T), function(x) {
#   readRDS(x)
# })

matching_samples = readxl::read_excel('data/scRNA_samples.xlsx', sheet = 1) %>% 
  as.data.frame() %>% 
  dplyr::mutate(genomics_code = ifelse(genomics_code == 'NA', NA, genomics_code)) %>% 
  dplyr::mutate(fixed_name = ifelse(fixed_name == 'NA', NA, fixed_name))

samples_check = matching_samples %>% 
  dplyr::select(fixed_name, scRNA_sample) %>%
  dplyr::distinct()

new_mapping_experiment = readRDS("data/mapping_samples.rds")

all_samples_dict = full_join(matching_samples, new_mapping_experiment, 
                             by = join_by('genomics_code' == 'genomics_code', 
                                          'fixed_name' == 'fixed_name'), na_matches = 'na')

genes_karyo = readRDS('/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_all_genes_qc_v3.rds')
genes_karyo = genes_karyo %>% 
  mutate(karyotype = gsub(' ', '', karyotype))

# data =  lapply(list.files(data_path, full.names = T), function(p) {
#   x = readRDS(p)
# })
# data = Reduce('merge', data)

data = lapply(list.files(data_path, full.names = T), function(p) {
  
  x = readRDS(p)
  sample = x$sample %>% unique
  
  # x = x@assays$RNA$counts
  # x = as.matrix(x)
  
  # smad2 = x[which(rownames(x) == 'SMAD2'), ]
  
  genomic_sample = all_samples_dict %>% 
    filter(scRNA_sample == sample)
  
  if(nrow(genomic_sample > 0)) {
    genomic_sample = genomic_sample$fixed_name %>% unique
    
    karyo = genes_karyo %>% 
      filter(sample == genomic_sample) %>% 
      filter(hgnc_symbol == 'SMAD2') %>% 
      dplyr::select(sample, karyotype)
    
    meta = do.call("rbind", replicate(
      length(colnames(x)), karyo, simplify = FALSE))
    meta = meta %>% 
      mutate(cell_id = colnames(x@assays$RNA$counts)) %>% 
      mutate(new_id = paste(sample, cell_id, sep = '_'))
    
    # # colnames(x) = paste(genomic_sample, colnames(x), sep = '_')
    # colnames(x@assays$RNA$counts) = paste(genomic_sample, colnames(x@assays$RNA$counts), sep = '_')
    
    print(sample)
    renamed = RenameCells(x[['RNA']], new.names = paste(genomic_sample, colnames(x), sep = '_'))
    # 
    # x = x %>% 
    #   as.data.frame %>% 
    #   tibble::rownames_to_column('gene')
    
    return(list(meta = meta, data = renamed))
  }
})

data = Filter(Negate(is.null), data)
counts_all = lapply(data, function(x) {x$data})
merge(counts_all)


counts = counts %>% 
  tibble::column_to_rownames('gene')
counts = as.matrix(counts)
counts = Matrix(counts, sparse = TRUE)

metadata = lapply(data, function(x) {x$meta}) %>% 
  Reduce('full_join', .)

scexp = SingleCellExperiment(list(counts = counts))

# counts = data %>% 
#   do.call('full_join', .)


expr = setNames(nm = data$cell_id, object = data$smad2_expr)
metadata = data %>% 
  select(-smad2_expr)

saveRDS(expr, 'data/smad2_scRNA_expression.rds')
saveRDS(metadata, 'data/smad2_metadata.rds')

# single_geme

gene = "SMAD2"

dg = data %>% dplyr::filter(hgnc_symbol == gene)
meta = data %>% 
  filter(hgnc_symbol == gene) %>% 
  dplyr::select(sample, karyotype) %>%  
  distinct() %>% 
  mutate(karyotype = factor(karyotype, levels =  c("1:1", "1:0", "2:0", "2:1")))


Y = data %>% 
  filter(sample %in% meta$sample) %>% 
  dplyr::select(expr, hgnc_symbol, sample) %>% 
  distinct() %>% 
  tidyr::pivot_wider(values_from = expr, names_from = hgnc_symbol, values_fill = 0) %>% 
  tibble::column_to_rownames('sample') %>% 
  as.matrix()

dds = DESeqDataSetFromMatrix(countData = t(Y), colData = meta, design = ~karyotype)
dds = DESeq(dds)

res = results(dds)
hist(res$pvalue)





