
rm(list = ls())
library(patchwork)
library(tidyverse)

# Main lower part ####
pB = readRDS("jovoniR/img/sf_psinorm_stable_FALSE/DNA_ratio.rds") + labs(x = "Karyotype")
pC = readRDS("jovoniR/img/sf_psinorm_stable_FALSE/compensated_frac.rds")

pLower = pB + pC + 
  plot_layout(des = "AABBB") +
  plot_annotation(tag_levels = list(c("B", "C"))) &
  theme(plot.tag = element_text(face = "bold"))

ggsave("img/main1_lower.pdf", plot = pLower, width = 8, height = 3.8, units = "in")
