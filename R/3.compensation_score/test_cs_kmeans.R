library(tidyverse)

df_groups = readRDS('data/compensation/cs_classification_gene_labels.rds')

df = df_groups %>% 
  dplyr::select(name, RNA, Protein) %>% 
  distinct() %>% 
  tibble::column_to_rownames('name')

k_clustering = kmeans(df, centers = 5)

cluster_res = tibble(
  gene = k_clustering$cluster %>% names, 
  cluster = k_clustering$cluster %>% unname
)

df_clustered = df_groups %>% 
  dplyr::select(-c(class, label)) %>% 
  distinct() %>% 
  full_join(., cluster_res, by = join_by('name' == 'gene'))

df_clustered %>% 
  mutate(cluster = as.character(cluster)) %>% 
  ggplot(aes(
    RNA, 
    Protein, 
    color = cluster
  )) + 
  geom_point() + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  facet_wrap(~reg_group)

