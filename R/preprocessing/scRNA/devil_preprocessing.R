rm(list=ls())
.libPaths()
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(readxl)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')

data_path = 'data/scRNA'
list.files(data_path)


dict = readRDS('data/full_dict_dna_rna_prot.rds')

# genes_karyo = readRDS('/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_all_genes_qc_v3.rds')
# genes_karyo_muts = readRDS('data/karyotypes_mutations_all_genes_qc_ccf_v4.rds')
# genes_karyo_muts = genes_karyo_muts %>% 
#   mutate(karyotype = gsub(' ', '', karyotype))

# loading table filtered by ccf 
genes_karyo_muts = readRDS('data/processed_data/genes_filtered_karyo_mut_status_filt_ccf_08.rds')
genes_karyo_muts = genes_karyo_muts %>% 
  filter(sample != '11')

# genes_karyo_muts = readRDS('data/genes_cna_mut_status_filtered.rds')

# excluding pdo11 from the analysis
sc_samples = grep('Sample_PDO11_filtered', 
                  list.files(data_path, full.names = T), 
                  value =T, 
                  invert = T)

data = lapply(sc_samples, function(p) {
  
  x = readRDS(p)
  sample = x$sample %>% unique
  
  # x = x@assays$RNA$counts
  # x = as.matrix(x)
  
  # smad2 = x[which(rownames(x) == 'SMAD2'), ]
  
  genomic_sample = dict %>% 
    filter(RNA == sample)
  
  if(nrow(genomic_sample > 0)) {
    genomic_sample = genomic_sample$fixed_name %>% unique
    
    karyo = genes_karyo_muts %>% 
      filter(sample == genomic_sample) %>% 
      # filter(hgnc_symbol == 'SMAD2') %>% 
      dplyr::select(sample, hgnc_symbol, karyotype) %>% 
      distinct() 
    
    if(nrow(karyo) > 0){
      
      # meta = do.call("rbind", replicate(
      #   length(colnames(x)), karyo, simplify = FALSE))
      # meta = meta %>% 
      #   mutate(cell_id = colnames(x@assays$RNA$counts)) %>% 
      #   mutate(new_id = paste(sample, cell_id, sep = '_'))
      
      # # colnames(x) = paste(genomic_sample, colnames(x), sep = '_')
      # colnames(x@assays$RNA$counts) = paste(genomic_sample, colnames(x@assays$RNA$counts), sep = '_')
      
      # change the cell ids 
      print(sample)
      renamed = RenameCells(x[['RNA']], new.names = paste(genomic_sample, colnames(x), sep = '_'))
      # 
      # x = x %>% 
      #   as.data.frame %>% 
      #   tibble::rownames_to_column('gene')
      
      return(list(meta = karyo, data = renamed))}
  }
})

data = Filter(Negate(is.null), data)

common_genes = Reduce('intersect', lapply(data, function(x) {x$data@counts %>% rownames}))
counts_all_common_genes = lapply(data, function(x) {subset(x$data, features = common_genes)})
counts = merge(counts_all_common_genes[[1]], counts_all_common_genes[2:length(counts_all_common_genes)])
# 
# saveRDS(counts, 'data/scexp_karyo_all_pdo_counts_filtered.rds')
# counts = readRDS('data/scexp_karyo_all_pdo_counts_filtered.rds')

barcodes = colnames(counts)

samples_list = lapply(data, function(x) {x$meta %>% pull(sample) %>% unique}) %>% unlist
pattern_samples = paste(samples_list, collapse = "|")

metadata = tibble(cell_id = barcodes)
# pattern <- paste(unique(dict$fixed_name), collapse = "|")

# metadata_v1 <- metadata %>%
#   mutate(sample_id = str_extract(cell_id, pattern))

metadata <- metadata %>%
  mutate(sample_id = str_extract(cell_id, pattern_samples))

metadata = metadata %>%
  tibble::column_to_rownames('cell_id')

mapping_ICR_names <- read_excel("data/mapping_ICR_names.xlsx")

metadata = metadata %>% 
  dplyr::mutate(batch = case_when(
    sample_id %in% mapping_ICR_names$fixed_name ~ 'ICR', 
    .default = 'HSR')
  )
saveRDS(metadata, 'data/new_metadata_sc.rds')

counts_sc = counts$counts
counts_sc = counts_sc[,rownames(metadata)]
dim(counts_sc) 

scexp = SingleCellExperiment(list(counts = counts_sc), colData = metadata)
saveRDS(scexp, 'data/scexp_karyo_all_organoids_filt.rds')

samples = metadata$sample_id %>% unique
karyo = genes_karyo_muts %>% 
  filter(sample %in% samples) %>% 
  # filter(hgnc_symbol == 'SMAD2') %>% 
  dplyr::select(sample, hgnc_symbol, karyotype) %>% 
  distinct() 
saveRDS(karyo, 'data/karyotypes_genes_filtered_scrna.rds')

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





