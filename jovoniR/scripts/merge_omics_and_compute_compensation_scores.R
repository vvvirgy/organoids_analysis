
rm(list = ls())
library(tidyverse)
library(dplyr)
library(boot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
source("scripts/constants.R")
source("utils.R")

# Default values if not provided
sf_method  <- SF_METHOD
use_stable <- USE_STABLE

IMG_PATH = paste0("img/sf_", sf_method, "_stable_", use_stable)
RES_PATH = paste0("results/multiOmic")
#RNA_PATH = paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/lfc_res_",sf_method,"_stable_",use_stable,".rds")
RNA_PATH = "results/RNA/lfc_res_clean.rds"

dir.create(IMG_PATH, recursive = T)
dir.create(RES_PATH, recursive = T)

# Parameters
# Get gene/karyotypes with at least N samples
karyotypes_df_all = readRDS(META_PATH)
karyotypes_df_good = karyotypes_df_all %>%
  dplyr::group_by(hgnc_symbol, karyotype) %>%
  dplyr::distinct() %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n >= MIN_SAMPLES) %>%
  dplyr::rename(name = hgnc_symbol)

res_rna = readRDS(RNA_PATH) %>% dplyr::mutate(omic = "RNA") %>%
  dplyr::mutate(lfc = ifelse(abs(lfc) > 10, sign(lfc) * 10, lfc))

res_rna = res_rna %>% dplyr::select(!c(mean_expr, non_zero_percent, min_sample_mean, sample_means, n_samples))

res_prot = readRDS("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/processed_data/protein/fc_tb_clean_v2.rds") %>%
  dplyr::rename(coef = condition, name = PG.Genes, lfc = diff, pval = p.val, adj_pval = p.adj) %>%
  dplyr::mutate(omic = "Protein") %>%
  tidyr::separate(coef, sep = "X", into = c(".", "karyotype")) %>%
  dplyr::select(lfc, name, karyotype, omic, pval, adj_pval) %>%
  mutate(karyotype = str_replace(karyotype, "\\.", ":")) %>%
  dplyr::group_by(karyotype)

df_dna = readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/DNA_lfc.rds") %>% dplyr::rename(DNA_lfc = lfc)
# df_dna = df_dna %>%
#   dplyr::mutate(DNA_lfc = ifelse(karyotype == "2:2", 1, ifelse(karyotype == "1:0", -1, ifelse(karyotype == "2:0", 0, 0.5))))

df = dplyr::bind_rows(res_rna, res_prot) %>%
  dplyr::filter(!is.na(lfc))
df$omic = factor(df$omic, levels = c("RNA", "Protein"))

if (USE_ONLY_COMMON_GENES) {
  df = df %>%
    dplyr::group_by(karyotype, name) %>%
    dplyr::filter(n() == 2)
}

df %>% saveRDS(file.path(RES_PATH, "lfc_prot_and_rna.rds"))

df = df %>%
  dplyr::left_join(karyotypes_df_good)

df_dna = df_dna %>%
  dplyr::left_join(karyotypes_df_good)

df_with_compensation_scores = df %>%
  dplyr::filter(omic != "DNA") %>%
  ungroup() %>%
  dplyr::left_join(df_dna %>% dplyr::select(DNA_lfc, name, karyotype)) %>%
  dplyr::mutate(CS = ifelse(DNA_lfc > 0, DNA_lfc - lfc, lfc - DNA_lfc)) %>%
  dplyr::filter(!is.na(CS))

df_with_compensation_scores %>% saveRDS(file.path(RES_PATH, "CS_scores_prot_and_rna.rds"))

df_with_compensation_scores %>%
  ggplot(mapping = aes(x = karyotype, y = lfc, col = omic)) +
  geom_boxplot()


df_with_compensation_scores %>%
  ggplot(mapping = aes(x = DNA_lfc, y = lfc)) +
  geom_point() +
  facet_wrap(~omic)
