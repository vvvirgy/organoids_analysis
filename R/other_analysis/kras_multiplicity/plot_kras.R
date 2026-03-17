library(tidyverse)

# prot = readRDS('data/kras_multiplicity/fc_tb_clean_prot.rds')

prot = fc_tb_clean
rna = readRDS('data/kras_multiplicity/KRAS_RES/lfc.rds')

prot = prot %>% 
  filter(condition != 'X0') %>%
  # mutate(m_status = ifelse(condition == 'single', 'm = 1', 'm > 1')) %>% 
  mutate(mut_ratio = gsub('X', '', condition),
         mut_ratio = gsub('\\.', '/', mut_ratio)) %>% 
  dplyr::select(-c(condition)) %>% 
  mutate(p.adj = p.adjust(p.val, method = 'BH')) %>% 
  mutate(omic = 'Proteomics') %>% 
  rename(name = 'PG.Genes') %>% 
  rename(ci_low = "CI.L") %>% 
  rename(ci_high = "CI.R") %>% 
  rename(lfc = diff) %>% 
  rename(adj_pval = p.adj) %>% 
  rename(pval = p.val) %>% 
  mutate(significant = ifelse(adj_pval <= 0.05, 'significant (adj_pval <= 0.05)', 'ns'))

rna = rna %>% 
  mutate(omic = 'Transcriptomics') %>% 
  select(-se_lfc) %>% 
  rename(mut_ratio = coef) %>% 
  mutate(significant = ifelse(adj_pval <= 0.05, 'significant (adj_pval <= 0.05)', 'ns'))
  
setdiff(colnames(rna), colnames(prot))
rna = rna[,colnames(prot)]
lv = c('0/1', '0/2', '0/3', '0/4', '0/5', '1/2', "2/2", '1/3', '2/3', '3/3', '1/4', '2/4', '4/4', '5/5', '3/6', '6/6')

df = bind_rows(rna, prot) %>% 
  mutate(mut_ratio = factor(mut_ratio, levels = lv)) %>% 
  mutate(omic = factor(omic, levels = c('Transcriptomics', 'Proteomics')))

cols = setNames(
  nm = df$significant %>% unique, 
  object = c('#C4DAD2', '#018790')
)

df %>%
  ggplot(aes(
    y = mut_ratio, 
    x = lfc, 
    fill = significant
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  labs(
    x = 'log2FC', 
    y = 'm/k KRAS'
  ) + 
  facet_wrap(~omic, scales = 'free_x') + 
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#740A03') + 
  scale_fill_manual(values = cols) + 
  theme(legend.position = 'bottom') 
ggsave('res/kras_mk/img.png', width = 8, height = 7, dpi = 600)

saveRDS(df, 'data/kras_multiplicity/multi_omics_lfc.rds')

prot %>% 
  mutate(mut_ratio = factor(mut_ratio, levels = lv)) %>% 
  ggplot(aes(
    y = mut_ratio, 
    x = lfc, 
    fill = significant
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  labs(
    x = 'log2FC', 
    y = 'm/k KRAS'
  ) + 
  # facet_wrap(~omic, scales = 'free_x') + 
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#740A03') + 
  scale_fill_manual(values = cols) + 
  theme(legend.position = 'bottom') 
