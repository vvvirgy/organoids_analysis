
rm(list = ls())
library(tidyverse)

df = read.delim("results/sf_psinorm_stable_FALSE/gene_CS_groups.csv", sep = ",")

groups = unique(df$reg_group)
g = groups[3]

gene_list = df %>% dplyr::filter(reg_group == g) %>% dplyr::pull(name)

library(clusterProfiler)
library(org.Hs.eg.db)  # for human

# Gene Ontology enrichment
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,
                ont = "ALL", keyType = "SYMBOL",
                pAdjustMethod = "BH",
                readable = TRUE)

# KEGG pathway enrichment
kegg <- enrichKEGG(gene = gene_list, keyType = "SYMBOL",
                   organism = 'hsa',
                   pvalueCutoff = 0.05)

# Reactome pathways
library(ReactomePA)
reactome <- enrichPathway(gene = gene_list)


library(biomaRt)
library(GenomicRanges)

# Get gene coordinates
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_coords <- getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                                    "start_position", "end_position"),
                     filters = "hgnc_symbol",
                     values = gene_list,
                     mart = mart)

# Test for chromosomal clustering
library(regioneR)
gr <- makeGRangesFromDataFrame(gene_coords, 
                               seqnames.field = "chromosome_name",
                               start.field = "start_position",
                               end.field = "end_position")

# Test against random expectation
pt <- permTest(A = gr, ntimes = 1000, 
               evaluate.function = numOverlaps,
               randomize.function = randomizeRegions)



library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get sequences and analyze
sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)

# GC content
gc_content <- letterFrequency(sequences, "GC", as.prob = TRUE)

# Secondary structure propensity (requires RNAfold)
# Look for highly structured regions that might escape degradation

# miRNA binding sites
library(multiMiR)
mirna_targets <- get_multimir(org = "hsa", 
                              target = gene_list,
                              table = "all")



library(STRINGdb)

string_db <- STRINGdb$new(version="11", species=9606, 
                          score_threshold=400, input_directory="")

# Map and get network
string_mapped <- string_db$map(data.frame(gene = gene_list), 
                               "gene", removeUnmappedRows = TRUE)
string_interactions <- string_db$get_interactions(string_mapped$STRING_id)

# Network clustering
library(igraph)
g <- graph_from_data_frame(string_interactions)
clusters <- cluster_louvain(g)