# comparing the fits between RNA and proteomics - good candidates
library(tidyverse)
library(UpSetR)

proteomics_genes_fit = readRDS('res/protein_vs_CNA/interesting_genes_fit.rds')
proteomics_genes_fit = proteomics_genes_fit[-good_best_up_prot]

transcriptomics_genes_fit = readRDS('res/RNA_vs_CNA/interesting_genes_fit.rds')

names(proteomics_genes_fit) = paste(names(proteomics_genes_fit), 'prot', sep = '_')
names(transcriptomics_genes_fit) = paste(names(transcriptomics_genes_fit), 'rna', sep = '_')

good_candidates = c(proteomics_genes_fit, transcriptomics_genes_fit)

good_candidates = lapply(good_candidates, function(x) {
  x$gene %>% unique
})

good_candidates = good_candidates[-c(3,4,10,11)]

 "cgc_colon_prot"    "cgc_somatic_prot"  "drivers_coad_prot" "good_down_rna"    
[7] "good_up_rna"       "cgc_colon_rna"     "cgc_somatic_rna"   "drivers_coad_rna" 


names(good_candidates) = c(
  'Proteomics_vs_CNA_positive', 
  'Proteomics_vs_CNA_negative',
  'Proteomics_vs_CNA_CGC_COAD', 
  'Proteomics_vs_CNA_CGC_somatic',
  'Proteomics_vs_CNA_COAD_Intogen',
  
  'Transcriptomics_vs_CNA_negative',
  'Transcriptomics_vs_CNA_positive', 
  'Transcriptomics_vs_CNA_CGC_COAD', 
  'Transcriptomics_vs_CNA_CGC_somatic',
  'Transcriptomics_vs_CNA_COAD_Intogen'
)

metadata = tibble(
  set = names(good_candidates)
) %>% 
  mutate(Assay = case_when(str_detect(set, 'Transcriptomics') ~ 'RNAseq', .default = 'Proteomics'))

pdf('res/intersection_fit_res_clean.pdf', width = 12, height = 8, onefile = F)
upset(data = UpSetR::fromList(good_candidates), 
      nsets = 10, 
      order.by = 'freq', 
      text.scale = c(1.3, 1.3, 1.3, 1.3, 1.5, 1.3),
      queries = list(
        list(
          query = intersects,
          params = list("Proteomics_vs_CNA_CGC_somatic", "Transcriptomics_vs_CNA_CGC_somatic", "Proteomics_vs_CNA_positive", 'Transcriptomics_vs_CNA_positive'), 
          color = "orange", active = T),
        list(query = intersects,
             params = list("Proteomics_vs_CNA_positive", 'Transcriptomics_vs_CNA_positive'), color = "forestgreen", active = T),
        list(query = intersects, 
             params = list("Transcriptomics_vs_CNA_negative", 'Proteomics_vs_CNA_negative'), color = '#A62C2C', active = T)#, 
        # list(query = intersects,
        #      params = list("good_down_prot", "good_down_rna"), color = "forestgreen", active = T)
      ),
      set.metadata = list(
        data = metadata, 
        plots = list(list(type = 'heat', 
                          column = 'Assay', 
                          assign = 10, 
                          colors = c('RNAseq' = '#129990', 
                                        'Proteomics' = '#27548A')))))
dev.off()
