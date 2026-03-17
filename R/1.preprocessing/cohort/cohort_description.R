library(tidyverse)
library(ggupset)
# install.packages('ggupset')
library(UpSetR)
library(gt)
dict = readRDS('data/full_dict_dna_rna_prot.rds')

presence = dict %>% 
  mutate(seq_dna = 
           ifelse(!is.na(DNA), 1, 0)) %>% 
  mutate(seq_rna = 
           ifelse(!is.na(RNA), 1, 0)) %>% 
  mutate(proteomics = 
           ifelse(!is.na(sample_proteomics_code), 1, 0)) %>% 
  dplyr::select(fixed_name, seq_dna, seq_rna, proteomics) %>% 
  distinct() %>% 
  rename(Proteomics = proteomics) %>% 
  rename(scRNA = seq_rna) %>% 
  rename(WGS = seq_dna)


nsamples = tibble(
  omic = c('WGS', 'scRNA-seq', 'Proteomics'), 
  n_samples = c(
    dict$DNA %>% unique %>% length, 
    dict %>% 
      filter(!is.na(RNA)) %>% 
      pull(RNA) %>% 
      unique %>% 
      length(), 
    dict %>% 
      filter(!is.na(sample_proteomics_code)) %>% 
      pull(proteomics_code) %>% 
      unique %>% 
      length()
  )
)

tb = nsamples %>% 
  gt(., rowname_col = "omic") %>% 
  tab_header(
    title = "ICR-HSR cohort") %>% 
  tab_stubhead(label = "Omic") %>% 
  cols_label(
    n_samples = 'Number of samples'
  ) %>% 
  cols_width(everything() ~ px(150)) %>% 
  tab_options(heading.align = 'center', 
              table.width = 8, page.height = '8') %>% 
  cols_align(align = 'center', columns = 'n_samples') %>% 
  gt::as_gtable(. , text_grob = gridtext::richtext_grob, plot = T) %>% 
  ggplotify::as.ggplot()

ggsave('res/tb_cohort_description.pdf', width = 8)

tb_plot = function(data) {
  data %>% 
    gt(., rowname_col = "omic") %>% 
    tab_header(
      title = "ICR-HSR cohort") %>% 
    tab_stubhead(label = "Omic") %>% 
    cols_label(
      n_samples = 'Number of samples'
    ) %>% 
    cols_width(everything() ~ px(150)) %>% 
    tab_options(heading.align = 'center', 
                table.width = 8, page.height = '8') %>% 
    cols_align(align = 'center', columns = 'n_samples') %>% 
    gt::as_gtable(. , text_grob = gridtext::richtext_grob, plot = T) %>% 
    ggplotify::as.ggplot()
}

pdf('res/cohort_upset_different_omics.pdf', width = 8, height = 6)
upset(presence %>% as.data.frame(), 
      order.by = 'freq', 
      matrix.color = '#088395',
      main.bar.color = '#088395', 
      sets.bar.color = '#088395', 
      mainbar.y.label = 'Number of samples', 
      point.size = 4, 
      line.size =1, 
      sets.x.label = 'Samples per omic', 
      text.scale = c(2, 2, 1.5, 1.5, 2, 2)#, 
      # attribute.plots = list(gridrows = 50, plots = list(list(plot = tb_plot(nsamples))))
#        c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
      
      ) 

dev.off()

png('res/cohort_upset_different_omics.png', width = 15, height = 12, units = 'cm', res = 300)
upset(presence %>% as.data.frame(), 
      order.by = 'freq', 
      matrix.color = '#088395',
      main.bar.color = '#088395', 
      sets.bar.color = '#088395', 
      mainbar.y.label = 'Number of samples', 
      point.size = 4, 
      line.size =1, 
      sets.x.label = 'Samples per omic', 
      text.scale = c(2, 2, 1.5, 1.5, 2, 2), 
      #        c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
      
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
  # theme_light() +
  theme_bw() + 
  theme_combmatrix(combmatrix.panel.point.color.fill = "#088395",
                   combmatrix.label.extra_spacing = 12,
                   # combmatrix.panel.line.size = 0,
                   combmatrix.panel.line.color = '#088395',
                   combmatrix.label.make_space = TRUE) + 
  labs(y = 'Number of samples')
  


