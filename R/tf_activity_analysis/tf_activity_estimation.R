rm(list=ls())
.libPaths()
library(SingleCellExperiment)
library(Seurat)
# BiocManager::install('saezlab/decoupleR')
library(decoupleR)
library(OmnipathR)
library(clusterProfiler)
library(org.Hs.eg.db)

data_path = 'data/scRNA'
list.files(data_path)

dict = readRDS('data/full_dict_dna_rna_prot.rds')

genes_karyo_muts = readRDS('data/processed_data/genes_filtered_karyo_mut_status_filt_ccf_08.rds')
genes_karyo_muts = genes_karyo_muts %>% 
  filter(sample != '11')

sc_samples = grep('Sample_PDO11_filtered', 
                  list.files(data_path, full.names = T), 
                  value =T, 
                  invert = T)

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)

# aggregate all the counts
aggregated_counts_cohort = lapply(sc_samples, function(p) {
  
  data = readRDS(p)
  sample = data$sample %>% unique
  
  ggenes = grep('ENS', rownames(data), value = T, invert = T)
  
  data = data[ggenes, ]
  # data = filter_sce(data)
  
  data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  
  mat = AggregateExpression(data, assays = 'RNA', slot = 'data')
  mat = mat[[1]]
  colnames(mat) = sample
  
  mat = mat %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('gene')
  
  return(mat)
})

# now merge all the aggregated counts by common genes
common_genes = Reduce(f = 'intersect', lapply(aggregated_counts_cohort, function(x){x$gene}))

counts = lapply(aggregated_counts_cohort, function(x) {
  x %>% 
    filter(gene %in% common_genes)
}) 
counts = Reduce(f = 'full_join', counts) 
counts = counts %>% 
  tibble::column_to_rownames('gene') %>% 
  as.matrix()

# Run ulm
acts <- decoupleR::run_ulm(mat = counts, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor='mor', 
                           minsize = 5)

tf_acts <- acts %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') 

saveRDS(tf_acts, 'data/tf_activity/TF_activity_by_sample')

pdos = dict %>% 
  filter(RNA %in% unique(tf_acts$condition)) %>% 
  dplyr::select(RNA, fixed_name) %>% 
  distinct()
  
tf_activity_cohort = tf_acts %>% 
  left_join(., pdos, by = join_by('condition' == 'RNA')) %>% 
  dplyr::select(-condition) %>% 
  relocate(fixed_name, .before = ABL1) %>% 
  rename(sample = fixed_name) %>% 
  pivot_longer(cols = -sample, values_to = 'TF_activity', names_to = 'TF')

tf_activity_karyo = genes_karyo_muts %>% 
  filter(sample %in% unique(tf_activity_cohort$sample)) %>%
  right_join(., tf_activity_cohort, 
             by = join_by(
               'hgnc_symbol' == 'TF', 
               'sample' == 'sample'
               )) %>% 
  filter(!is.na(sample))  %>% 
  filter(!is.na(karyotype)) %>% 
  rename(TF = hgnc_symbol)

saveRDS(tf_activity_karyo, 'data/tf_activity/tf_activity_with_karyotypes.rds')
tf_activity_karyo = readRDS('data/tf_activity/tf_activity_with_karyotypes.rds')

tf_activity_karyo %>% 
  filter(TF == 'RUNX3') %>% 
  ggplot(aes(
    y = TF_activity, 
    x = karyotype
  )) + 
  geom_boxplot()

FOXM1 = tf_acts %>% 
  dplyr::select(condition, FOXM1)
  
FOXM1_targets = net %>% 
  filter(source == 'FOXM1') %>% 
  pull(target)

FOXM1_targets_karyo = genes_karyo_muts %>% 
  filter(hgnc_symbol %in% FOXM1_targets)

tf_activity_karyo %>% 
  filter(TF == 'AHRR')



# tf_activity_karyo %>% 
#   ggplot(aes(
#     x = sample, 
#     y = reorder(hgnc_symbol, +mean), 
#     fill = mean
#   )) + 
#   geom_tile() + 
#   # facet_wrap(~karyotype, scales = 'free', ncol = 1)
#   facet_grid(vars(karyotype), scales = 'free_y', 
#              space = "free_y") + 
#   scale_fill_gradient2(
#     low = 'dodgerblue4', 
#     mid = 'white', 
#     high = 'firebrick', 
#     midpoint = 0,
#     guide = 'colourbar', 
#     name = 'TF activity'
#     ) + 
#   theme_bw() + 
#   labs(
#     y = 'TF', 
#     x = 'Sample'
#   ) 



tf_activity_karyo %>% 
  ggplot(aes(
    x = sample, 
    y = reorder(TF, +TF_activity), 
    fill = TF_activity
  )) + 
  geom_tile() + 
  # facet_wrap(~karyotype, scales = 'free', ncol = 1)
  facet_grid(vars(karyotype), scales = 'free_y', 
             space = "free_y") + 
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 10, name = 'RdBu')),
    # colours = c('dodgerblue4', 'white', 'firebrick'), 
    rescaler = ~ scales::rescale_mid(.x, mid = 0), 
    name = 'TF activity'
  ) + 
  theme_bw() + 
  labs(
    y = 'TF', 
    x = 'Sample'
  ) 


# # checking the expression of JUN 
# jun_up_targets = net %>% 
#   filter(source == 'JUN') %>% 
#   filter(mor == 1) %>% 
#   pull(target) %>%
#   unique
# 
# jun_targets_activity = mat[rownames(mat) %in% jun_up_targets, ] %>% 
#   as.data.frame()
# 
# jun_targets_activity %>% 
#   tibble::rownames_to_column('gene') %>% 
#   reshape2::melt() %>% 
#   group_by(gene) %>% 
#   summarise(expr = sum(value)) %>% 
#   filter(expr != 0) %>% 
#   slice_max(expr, n = 20) %>% 
#   pull(gene) -> ggenes
# 
#   
# genes_karyo_muts %>% 
#   filter(sample == pdo) %>% 
#   filter(hgnc_symbol %in% ggenes)

