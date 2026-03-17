
rm(list = ls())
library(patchwork)
library(tidyverse)

# Main lower part ####
p1 = readRDS("jovoniR/img/sf_psinorm_stable_FALSE/DNA_ratio.rds") + labs(x = "Karyotype") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p1@layers$geom_boxplot$aes_params$fill = alpha("#66c2a5", alpha = .5)
p2 = readRDS("jovoniR/img/sf_psinorm_stable_FALSE/volcano_sub.rds")
p3 = readRDS("jovoniR/img/sf_psinorm_stable_FALSE/dge_barplot_sub.rds")
p4 = readRDS("jovoniR/img/sf_psinorm_stable_FALSE/compensated_frac.rds")

n_annotation = p3@data %>% dplyr::filter(FC_class != "Not differential") %>% 
  dplyr::select(FC_class, karyotype, omic, n) %>% 
  dplyr::mutate(label = paste0("n=", n)) %>% 
  dplyr::mutate(x = 7.5, y = ifelse(omic == "RNA", 280, 12)) %>% 
  dplyr::mutate(y = ifelse(FC_class == "Up", y * 0.9, y))

p2 = p2 +
  geom_text(data = n_annotation, mapping = aes(x = x, y = y, col = FC_class, label = label), inherit.aes = F, show.legend = F)

des = "
AAAA
BBBC
BBBC
DDDD
"

des = "
ABBB
ABBB
CDDD
CDDD
"

des = "ABBCC"


p2
p3
p4

pLower = p1 + free(p2) + free(p4) +
  plot_layout(des = des) +
  plot_annotation(tag_levels = list(c("B", "C", "D"))) &
  theme(plot.tag = element_text(face = "bold"))
pLower
ggsave("img/omic_abundances_and_compensation.pdf", plot = pLower, width = 12, height = 4.7, units = "in")
ggsave("img/omic_abundances_and_compensation.png", plot = pLower, width = 12, height = 4.7, units = "in", dpi = 450)




p_volcan_and_bar = p_volcano + p_bar + 
  plot_layout(ncol = 1, nrow = 2, heights = c(2,1.1)) +
  plot_annotation(tag_levels = c("A")) & theme(plot.tag = element_text(face = "bold"))

ggsave("../img/volcano_and_bar.pdf", plot = p_volcan_and_bar, width = 10, height = 7, units = "in")
ggsave("../img/volcano_and_bar.png", plot = p_volcan_and_bar, width = 10, height = 7, units = "in", dpi = 450)