library(tidyverse)
library(ggforestplot)

rna_fit = readRDS('data/glm_fit_v2_rna.rds')
rna_fit_coeffs = lapply(rna_fit %>% names, function(s) {
  broom::tidy(rna_fit[[s]]) %>% 
    dplyr::mutate(gene = s)
}) %>% bind_rows()

coeffs = rna_fit_coeffs %>% 
  dplyr::filter(term != '(Intercept)')

rr = ggforestplot::forestplot(
  df = coeffs,
  name = term,
  estimate = estimate,
  se = std.error, 
  # pvalue = p.value,
  # psignif = 0.002
  xlab = 'Estimated effect of mutation multiplicity 
  on the RNA expression',
  title = 'Genes multiplicity - RNA'
)
ggsave('res/rna_forest_plot.png', rr, width = 6, height = 5)

rna_msh6 = rna_fit_coeffs %>% 
  filter(gene == 'MSH6')


prot_fit = readRDS('data/glm_fit_v2_prot.rds')
prot_fit = lapply(prot_fit %>% names, function(s) {
  broom::tidy(prot_fit[[s]]) %>% 
    dplyr::mutate(gene = s)
}) %>% bind_rows()

coeffs_prot = prot_fit %>% 
  dplyr::filter(term != '(Intercept)')

pp = ggforestplot::forestplot(
  df = coeffs_prot,
  name = term,
  estimate = estimate,
  se = std.error, 
  # pvalue = p.value,
  # psignif = 0.002
  xlab = 'Estimated effect of mutation multiplicity 
  on the protein expression',
  title = 'Genes multiplicity - protein'
)
ggsave('res/prot_forest_plot.png', pp, width = 6, height = 5)

rr + pp

smad2_rna = coeffs %>% 
  # tidyr::pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::filter(gene == 'SMAD2') %>% 
  dplyr::mutate(assay = 'RNA')
  #%>% 
  # ggplot(aes(x = term, y = estimate, fill = term)) + 
  # geom_bar(stat = 'identity') + 
  # theme_light() + 
  # scale_fill_brewer(palette = 'GnBu') + 
  # ggtitle('SMAD2 betas - RNA')

smad2_prot = coeffs_prot %>% 
  # tidyr::pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::filter(gene == 'SMAD2') %>% 
  mutate(assay = 'Protein')
  # ggplot(aes(x = term, y = estimate, fill = term)) + 
  # geom_bar(stat = 'identity') + 
  # theme_light() + 
  # scale_fill_brewer(palette = 'GnBu') + 
  # ggtitle('SMAD2 betas - protein')

smad2 = bind_rows(smad2_rna, smad2_prot)

smad2 %>% 
  mutate(assay = factor(assay, levels = c('RNA','Protein'))) %>% 
  ggplot(aes(x = term, y = estimate, fill = assay)) +
  geom_bar(stat = 'identity') +
  # theme_light() +
  theme_bw()+
  # scale_fill_brewer(palette = '') +
  # ggsci::scale_fill_frontiers() +
  scale_fill_manual(values = c('RNA' = 'steelblue', 'Protein' = 'goldenrod')) +
  ggtitle('SMAD2 betas') + 
  facet_wrap(vars(assay), scales = 'free') + 
  xlab('predictor') + 
  ylab('beta')


smad2_rna %>% 
  # mutate(assay = factor(assay, levels = c('RNA','Protein'))) %>% 
  ggplot(aes(x = term, y = estimate, fill = assay)) +
  geom_bar(stat = 'identity') +
  # theme_light() +
  theme_bw()+
  # scale_fill_brewer(palette = '') +
  # ggsci::scale_fill_frontiers() +
  scale_fill_manual(values = c('RNA' = 'steelblue')) + #, 'Protein' = 'goldenrod')) +
  # ggtitle('RNA') + 
  facet_wrap(vars(assay), scales = 'free') + 
  xlab('predictor') + 
  ylab('beta')
ggsave('res/smad2_beta_rna.png', width = 5, height = 3, bg = 'white')

smad2_prot %>% 
  # mutate(assay = factor(assay, levels = c('RNA','Protein'))) %>% 
  ggplot(aes(x = term, y = estimate, fill = assay)) +
  geom_bar(stat = 'identity') +
  # theme_light() +
  theme_bw()+
  # scale_fill_brewer(palette = '') +
  # ggsci::scale_fill_frontiers() +
  scale_fill_manual(values = c('Protein' = 'goldenrod')) +
  # ggtitle('Protein') + 
  facet_wrap(vars(assay), scales = 'free') + 
  xlab('predictor') + 
  ylab('beta')
ggsave('res/smad2_beta_prot.png', width = 5, height = 3, bg = 'white')
# (smad2_rna + smad2_prot) + 
#   plot_layout(guides = 'collect') & 
#   theme(legend.position = 'bottom')
