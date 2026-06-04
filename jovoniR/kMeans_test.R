
library(VIM)

df = readRDS("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/compensation/cs_karyotype.rds")

df_wide = df %>% 
  dplyr::select(name, karyotype, omic, CS) %>% 
  tidyr::pivot_wider(names_from = c(omic, karyotype), values_from = CS)

df_imputed <- kNN(df_wide, k = 5)
df_imputed <- df_imputed[, colnames(df_wide)]

library(factoextra)

# 1. Keep only numeric columns (drop 'name' identifier)
df_numeric <- df_imputed %>% 
  tibble::column_to_rownames("name") %>%  # use name as row label
  dplyr::select(where(is.numeric))

# 2. Scale
df_scaled <- scale(df_numeric)

# 3. Choose k
fviz_nbclust(df_scaled, kmeans, method = "silhouette")

# 4. Fit k-means
set.seed(42)
km <- kmeans(df_scaled, centers = 2, nstart = 25)

# 5. Attach clusters back
df_imputed$cluster <- km$cluster

# 6. Visualize
fviz_cluster(km, data = df_scaled)

# 7. Inspect cluster means
aggregate(df_numeric, by = list(cluster = km$cluster), FUN = mean)

df_clusters = dplyr::tibble(name = names(km$cluster), cluster = km$cluster)
df_classes = readRDS("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/compensation/cs_classification.rds")

df_classes %>% 
  dplyr::left_join(df_clusters) %>% 
  ggplot(mapping = aes(x = RNA, y = Protein, col = as.factor(cluster))) +
  geom_point(alpha = 0.2) +
  geom_density2d()
