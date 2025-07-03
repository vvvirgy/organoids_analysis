library(tidyverse)
library(msqrob2)

raw_data = readxl::read_excel('~/Downloads/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
  as.matrix

raw_data[which(is.na(raw_data[,'PG.Genes'])), 'PG.Genes'] = 'Unknown'

rownames(raw_data) = raw_data[,'PG.Genes']

median()

ecols = colnames(raw_data)[-c(1,2)]

raw_data_v2 = raw_data[,-c(1,2)]

raw_data_v2 = apply(raw_data_v2, 1, function(s) {
  as.numeric(s)
}) %>% t

colnames(raw_data_v2) = colnames(raw_data)[-c(1,2)]
raw_data_v2 = as.data.frame(raw_data_v2)

# x2 = x2 %>% 
#   mutate(Gene = rownames(x2)) %>%
#   dplyr::relocate(Gene, .before = '2L_a')
  
# pe = readQFeatures(assayData = x2, 
#               fnames = 1, 
#               quantCols = ecols)
# 
# rowData(pe[["quants"]])$nNonZero <- rowSums(assay(pe[["quants"]]) > 0)
# pe <- zeroIsNA(pe, "quants") 
# 
# # MSnbase::plotNA(assay(pe[["quants"]])) +
# #   xlab("Peptide index (ordered by data completeness)")
# 
# pe <- logTransform(pe, base = 2, i = "quants", name = "logData")
# pe2 <- filterFeatures(pe, ~ nNonZero >= 9628)
# 
# nrow(pe2[['logData']])

# pe2 <- normalize(pe2,
#                 i = "logData",
#                 name = "normData",
#                 method = "center.median"
# )
# 
# assay(pe2[['normData']])

all_genes = apply(raw_data_v2, 1, function(s) {
  # any(
    any(is.nan(as.numeric(s))) # == FALSE)
})

ok_genes = all_genes[which(all_genes == FALSE)] %>% names

reduced_genes = raw_data_v2[which(rownames(raw_data_v2) %in% ok_genes), ]

log2transformed = apply(reduced_genes, 1, function(x) {log(as.numeric(x), base = 2)})
rownames(log2transformed) = colnames(reduced_genes)

median_samples = apply(log2transformed, 1, median)
hist(median_samples - median(median_samples))

tt = log2transformed %>% 
  ggplot(aes(KRAS), alpha = 0.5) +
  geom_histogram(binwidth = 0.1) + 
  geom_vline(xintercept = median(log2transformed[,'KRAS']), colour = 'red') + 
  theme_bw() +
  xlim(c(14, 18))

tt/pp

hist(reduced_genes['KRAS', ] %>% as.numeric())

hist((log2transformed['57-VII_b',] - median(log2transformed['57-VII_b',]))/mad(log2transformed['57-VII_b',]))


# ppe <- PhosphoExperiment(assays = list(Quantification = t(log2transformed)),
#                          # Site = as.numeric(gsub("[A-Z]","", sapply(strsplit(rownames(ppe), ";"), "[[", 3))),
#                          GeneSymbol = sapply(strsplit(rownames(ppe), ";"), "[[", 1),
#                          # Residue = gsub("[0-9]","", sapply(strsplit(rownames(ppe), ";"), "[[", 3)),
#                          # Sequence = sapply(strsplit(rownames(ppe), ";"), "[[", 4)
#                          )

