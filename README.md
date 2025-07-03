# Organoids analysis

Analysis of scRNA, WGS and bulk proteomics data from organoids

The main idea is to be able to define a relationship between CNA alterations and differences in the expression

Analysis overview (divided by thematic parts)

-   Part 1 - **preparation of expression data**
    -   mapping of sample names between proteomics/transcriptomics samples and genomic samples
    -   TIC normalization of proteomics data
    -   Generation of pseudobulk and vst normalization (via DESEq2) from scRNA data
-   Part 2 - **preparation of genomic data**
    -   Extraction of genomic positions for each gene present in both datasets
    -   Extraction of ploidy for each gene (from CNAqc objects)
    -   Extraction of mutational status (both mutation type and if mutated/wild type) from mutational data
-   Part 3 - **First integration**
    -   Integration of expression and mutational data (CNA + mut)
    -   linear model fitting
    -   generalised linear model fitting (expression \~ cna + mut)
    -   visualization
-   Part 4 - **second integration**
    -   Integration of expression data
    -   model comparison
