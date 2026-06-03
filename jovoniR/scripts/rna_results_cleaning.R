
rm(list = ls())
library(tidyverse)

df = readRDS("results/RNA/lfc_res_psinorm_stable_FALSE.rds")
df$min_sample_mean <- sapply(df$sample_means, min)

df %>%
  ggplot(mapping = aes(x = karyotype, y = lfc)) +
  geom_boxplot()

df %>%
  ggplot(mapping = aes(x = mean_expr)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap(~karyotype)

df %>%
  ggplot(mapping = aes(x = lfc, y = mean_expr)) +
  geom_point() +
  scale_y_log10() +
  facet_grid(n_samples~karyotype)

# df_sub = df %>%
#   dplyr::filter(n_samples >= 5) %>%
#   dplyr::filter(mean_expr >= .5) %>%
#   dplyr::filter(min_sample_mean >= .05) %>%
#   dplyr::filter(non_zero_percent > 1)

df_sub = df %>%
  dplyr::filter(n_samples >= 4) %>%
  dplyr::filter(mean_expr >= .1) %>%
  dplyr::filter(min_sample_mean >= .05) %>%
  dplyr::filter(non_zero_percent > 1)

df_sub %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(m = mean(lfc))

df_sub %>%
  ggplot(mapping = aes(x = karyotype, y = lfc)) +
  geom_boxplot()

df_sub %>%
  ggplot(mapping = aes(x = lfc, y = mean_expr)) +
  geom_point() +
  scale_y_log10() +
  facet_grid(n_samples~karyotype) +
  geom_vline(xintercept = 0, colour = "red")

unique(df_sub$name) %>% length()

df_sub %>% saveRDS("results/RNA/lfc_res_clean.rds")
