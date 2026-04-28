rm(list = ls())
gc()
.libPaths()
library(tidyverse)
library(ggsurvfit)
library(survival)
library(wesanderson)
library(survminer)
library(patchwork)
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/extended_analysis/utils.R')

data_path = '/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/TCGA'
cna = readRDS(file.path(data_path, 'cna_data.rds'))
# protein = readRDS('data/TCGA/protein_data.rds')
# rna = readRDS('data/TCGA/rna_data.rds')
muts = readRDS(file.path(data_path, 'muts_data.rds'))
clinical = readRDS(file.path(data_path, 'clinical_data.rds'))

# select specific samples for the survival 
# kras mutants

kras_mutant_samples = muts %>% 
  filter(Hugo_Symbol == 'KRAS') %>% 
  dplyr::select(Tumor_Sample_Barcode) %>% 
  tidyr::separate(Tumor_Sample_Barcode, into = c('Project', 'TSS', 'Participant', 'Sample'), sep = '-') %>% 
  distinct() %>% 
  mutate(id = paste(Project, TSS, Participant, Sample, sep = '-')) %>% 
  pull(id) %>% 
  unique

kras_pik3ca_mutant = muts %>% 
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>% 
  filter(Hugo_Symbol %in% c('PIK3CA', 'KRAS')) %>% 
  distinct() %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  mutate(n = n()) %>% 
  filter(n == 2) %>% 
  dplyr::select(Tumor_Sample_Barcode) %>% 
  tidyr::separate(Tumor_Sample_Barcode, into = c('Project', 'TSS', 'Participant', 'Sample'), sep = '-') %>% 
  distinct() %>% 
  mutate(id = paste(Project, TSS, Participant, Sample, sep = '-')) %>% 
  pull(id) %>% 
  unique

# clinical$clinical_follow_up_v1.0_coad

# treat cnas
sample_cols = colnames(cna)[6:ncol(cna)] 

kras_cna = cna %>% 
  filter(gene_name == 'KRAS') %>% 
  pivot_longer(names_to = 'sample',cols = sample_cols, values_to = 'CNA_state') %>%
  tidyr::separate(sample, into = 'id', sep = ';') %>% 
  mutate(id = gsub('_copy_number|_min_copy_number|_max_copy_number', '', id)) %>% 
  distinct() %>% 
  mutate(state = case_when(
    CNA_state > 2 ~ 'amp', 
    CNA_state < 2 ~ 'del', 
    .default = 'none'))

kras_amp = kras_cna %>% 
  dplyr::filter(state == 'amp') %>% 
  pull(id) %>% 
  unique

kras_mutant_samples = setdiff(kras_mutant_samples, kras_pik3ca_mutant)

kras_mut_amp = intersect(kras_mutant_samples, kras_amp)

# test
muts %>% 
  tidyr::separate(Tumor_Sample_Barcode, into = c('Project', 'TSS', 'Participant', 'Sample'), sep = '-') %>% 
  distinct() %>% 
  mutate(id = paste(Project, TSS, Participant, Sample, sep = '-')) %>% 
  filter(id %in% kras_mut_amp) %>% 
  filter(Hugo_Symbol %in% c('KRAS', 'PIK3CA')) %>% 
  dplyr::select(Hugo_Symbol, id) %>% 
  view

kras_2n = kras_cna %>% 
  dplyr::filter(state == 'none') %>% 
  pull(id) %>% 
  unique

kras_mut_2n_pik3ca_mut = intersect(kras_pik3ca_mutant, kras_2n)
kras_mut_2n_pik3ca_mut = tibble(
  id = kras_mut_2n_pik3ca_mut
) %>%
  tidyr::separate(id, into = c('Project', 'TSS', 'Participant', 'Sample'), sep = '-') %>% 
  distinct() %>% 
  mutate(patient_id = paste(Project, TSS, Participant, sep = '-')) %>% 
  mutate(group = 'KRAS mut 2n + PIK3CA mut')

kras_mut_amp = tibble(
  id = kras_mut_amp
) %>%
  tidyr::separate(id, into = c('Project', 'TSS', 'Participant', 'Sample'), sep = '-') %>% 
  distinct() %>% 
  mutate(patient_id = paste(Project, TSS, Participant, sep = '-')) %>% 
  mutate(group = 'KRAS mut, amp')

sample_metadata = bind_rows(kras_mut_amp, kras_mut_2n_pik3ca_mut)

# get clinical data and merge with the mutation class
clinical_data = clinical$clinical_patient_coad

clinical_data = clinical_data[-c(1,2), ]
clinical_data = clinical_data %>% 
  right_join(., sample_metadata, by = join_by('bcr_patient_barcode' == 'patient_id'))

# run survival
clinical_data = clinical_data %>% 
  dplyr::select(bcr_patient_barcode, group, death_days_to, vital_status, lost_follow_up)


clinical_followup = clinical$clinical_follow_up_v1.0_coad
clinical_followup = clinical_followup[-c(1,2), ]
clinical_followup = clinical_followup %>% 
  right_join(., sample_metadata, by = join_by('bcr_patient_barcode' == 'patient_id')) %>% 
  group_by(bcr_patient_barcode) %>% 
  slice_max(form_completion_date) %>% 
  dplyr::select(bcr_patient_barcode, group, death_days_to, vital_status, last_contact_days_to) %>% 
  distinct() %>% 
  filter(vital_status != '[Not Available]') %>% 
  filter(death_days_to != '[Discrepancy]') %>% 
  mutate(OS_days = case_when(
    vital_status == 'Alive' ~ as.numeric(last_contact_days_to), 
    vital_status == 'Dead' ~ as.numeric(death_days_to)
  )) %>% 
  filter(OS_days > 0) %>% 
  mutate(OS_years = OS_days/365) %>% 
  dplyr::mutate(status = case_when(vital_status == 'Dead' ~ 1, 
                                   vital_status == 'Alive' ~ 0, 
                                   .default = NA))

clinical_followup %>% 
  filter(group == 'KRAS mut 2n + PIK3CA mut') %>% 
  filter(OS_years == 0)

test_survival = survival_analysis(clinical_followup, 
                  time = 'OS_years', 
                  event = 'status', 
                  what = 'group')
dev.off()

dir.create('res/survival_analysis', recursive = T)
ggsave(plot = plot(test_survival$km_plot), filename = 'res/survival_analysis/tcga_km_kras_amp_vs_kras_pik3ca_comut.png', width = 7, height = 6, dpi = 300, units = 'in')
ggsave(plot = plot(test_survival$km_plot), filename = 'res/survival_analysis/tcga_km_kras_amp_vs_kras_pik3ca_comut.pdf', width = 7, height = 6)





