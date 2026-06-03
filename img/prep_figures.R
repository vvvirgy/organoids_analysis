
library(tidyverse)
library(patchwork)

# MY_THEME = theme(text = element_text(size = 14))

# Main 1 ####
pA = readRDS("../jovoniR/img/sf_psinorm_stable_FALSE/DNA_ratio.rds")
pB = readRDS("../jovoniR/img/sf_psinorm_stable_FALSE/compensated_frac.rds")

des = "AABBB"

main_1 = pA + pB + plot_layout(design = des)
ggsave("figs/main1_bottom.pdf", plot = main_1, width = 8, height = 4, units = "in")


# Main 2 ####

pA = readRDS("../jovoniR/img/sf_psinorm_stable_FALSE/omics_trends.rds")
pB = readRDS("../jovoniR/img/sf_psinorm_stable_FALSE/tsg_og_reg_group_dist.rds")
pC = readRDS("../jovoniR/img/sf_psinorm_stable_FALSE/reg_groups_plot_annotated.rds")
pD = readRDS("../jovoniR/img/sf_psinorm_stable_FALSE/semantic_enrichment.rds")

des = "
AAAA
AAAA
BBDD
CCDD
CCDD
"

main_2 = pA + pB + pC + pD + plot_layout(design = des) + plot_annotation(tag_levels = c("A"))
ggsave("figs/main2.pdf", plot = main_2, width = 12, height = 10, units = "in")
ggsave("figs/main2.png", plot = main_2, width = 12, height = 10, units = "in", dpi = 450)
