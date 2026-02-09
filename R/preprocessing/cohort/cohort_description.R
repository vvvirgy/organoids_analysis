library(tidyverse)
library(ggupset)

dict = readRDS('data/full_dict_dna_rna_prot.rds')

# dna = dict %>% 
#   filter(!is.na(genomics_code)) %>% 
#   pull(genomics_code) %>% 
#   unique
# 
# prot = dict %>% 
#   filter(!is.na(sample_proteomics_code)) %>% 
#   pull(genomics_code) %>% 
#   unique
# 
# rna = dict %>% 
#   filter(!is.na(RNA)) %>% 
#   pull(genomics_code) %>% 
#   unique

# samples = list(DNA = dna, Protein = prot, RNA = rna)

presence = dict %>% 
  mutate(seq_dna = 
           ifelse(!is.na(DNA), 1, 0)) %>% 
  mutate(seq_rna = 
           ifelse(!is.na(RNA), 1, 0)) %>% 
  mutate(proteomics = 
           ifelse(!is.na(sample_proteomics_code), 1, 0)) %>% 
  select(fixed_name, seq_dna, seq_rna, proteomics) %>% 
  distinct() %>% 
  rename(Proteomics = proteomics) %>% 
  rename(scRNA = seq_rna) %>% 
  rename(WGS = seq_dna)


pdf('res/cohort_upset_different_omics.pdf', width = 8, height = 5)
upset(presence %>% as.data.frame(), 
      order.by = 'freq', 
      matrix.color = '#088395',
      main.bar.color = '#088395', 
      sets.bar.color = '#088395', 
      mainbar.y.label = 'Number of samples'
      ) 
dev.off()

png('res/cohort_upset_different_omics.png', width = 13, height = 10, units = 'cm', res = 300)
upset(presence %>% as.data.frame(),
      order.by = 'freq',
      matrix.color = '#088395',
      main.bar.color = '#088395',
      sets.bar.color = '#088395',
      mainbar.y.label = 'Number of samples'
)
dev.off()



presence = dict %>% 
  mutate(seq_dna = 
           ifelse(!is.na(DNA), 'DNA', NA)) %>% 
  mutate(seq_rna = 
           ifelse(!is.na(RNA), 'RNA', NA)) %>% 
  mutate(proteomics = 
           ifelse(!is.na(sample_proteomics_code), 'Proteomics', NA)) %>% 
  select(fixed_name, seq_dna, seq_rna, proteomics) %>% 
  distinct() %>% 
  rename(Proteomics = proteomics) %>% 
  rename(scRNA = seq_rna) %>% 
  rename(WGS = seq_dna)

presence %>% 
  pivot_longer(values_to = 'omic', cols = c(WGS, scRNA, Proteomics)) %>% 
  filter(!is.na(omic)) %>% 
  select(fixed_name, name) %>% 
  group_by(fixed_name) %>%
  summarise(name = list(name), .groups = "drop") %>% 
  rename(Omic = name) %>% 
  ggplot(aes(Omic)) + 
  geom_bar(fill = '#088395', width = .8)+ 
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset() + 
  theme_light() +
  theme_combmatrix(combmatrix.panel.point.color.fill = "#088395",
                   combmatrix.label.extra_spacing = 12,
                   # combmatrix.panel.line.size = 0,
                   combmatrix.panel.line.color = '#088395',
                   combmatrix.label.make_space = TRUE) + 
  labs(y = 'Number of samples')
  


