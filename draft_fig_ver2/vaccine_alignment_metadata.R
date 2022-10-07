library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(writexl)
library(patchwork)
source('/Users/prashant/R_themes/theme_publication.R')
source('/Users/prashant/R_themes/theme_publication_rnaseq.R')
source('epitope_alignment_p5_funs.R')

#######
#fig1a
all_data_dv = read_xlsx('india_dv_strains_ncbi_combined_108_strains_9feb.xlsx', 
                        sheet = 1)

#remove QVV22270.1 as it is showing no alignment with dengue3 epitopes
all_data_dv = all_data_dv[!all_data_dv$ncbiId %in% 'QVV22270.1',]

iedb_cd4 = read.csv('cd4_epitopes_pos_annotated_9th_feb.csv')
iedb_cd4$X = NULL
colnames(iedb_cd4)[colnames(iedb_cd4) == "Antigen.Name...11"] = "Antigen.Name"
iedb_cd4$Protein_name = gsub('Envelope', 'Env', iedb_cd4$Protein_name)
iedb_cd8 = read.csv('cd8_epitopes_pos_annotated.csv')
iedb_cd8$Protein_name = gsub('Envelope', 'Env', iedb_cd8$Protein_name)

iedb_cd4$group = 'CD4'
iedb_cd8$group = 'CD8'

all_data_dv_list = split(all_data_dv, f = all_data_dv$virusType)
iedb_cd4_list = split(iedb_cd4, f = iedb_cd4$DV_status)
iedb_cd8_list = split(iedb_cd8, f = iedb_cd8$DV_status)

align_info_iedb_cd4 = readRDS('rds_data/align_info_iedb_cd4_10feb.rds')
align_info_iedb_cd8 = readRDS('rds_data/align_info_iedb_cd8_10feb.rds')
#########################################
#save alignment metadata
cd4_metadata = lapply(align_info_iedb_cd4, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd4)
})
cd8_metadata = lapply(align_info_iedb_cd8, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd8)
})

experimental_data_cd4 = read.csv('tcell_table_export_1644390265_9_feb.csv')
colnames(experimental_data_cd4) = experimental_data_cd4[1,]
experimental_data_cd4 = experimental_data_cd4[-1,]
experimental_data_cd4 = experimental_data_cd4[c("Epitope ID", "Method/Technique", 
                                                "Assay Group", "Qualitative Measure")]
colnames(experimental_data_cd4)[1] = 'Epitope_ID'
experimental_data_cd4 = experimental_data_cd4[!duplicated(
  experimental_data_cd4$Epitope_ID),]

##
#add vaccines info
vac_files = list.files('draft_fig_ver2/supply/', pattern = 'metadata')
vac_files = c("denvax_metadata.xlsx", "dpiv_metadata.xlsx", 
              "nih_metadata.xlsx", "sen_metadata.xlsx", 
              "tden_metadata.xlsx", "tvdv_metadata.xlsx")

setwd('draft_fig_ver2/supply/')
vac_dat_cd4 = lapply(vac_files, function(x){
  dat = read_xlsx(x, sheet = 'CD4')
  dat = dat[c('Epitope_ID', 'value')]
  dat$value = as.numeric(dat$value)
  colnames(dat)[2] = gsub('_.*', '', x)
  return(dat)
})
setwd('..')
setwd('..')

temp = iedb_cd4[c('Epitope.ID', 'DV_status')]
colnames(temp)[1] = 'Epitope_ID'
vac_dat_cd4 = c(list(temp), vac_dat_cd4)
vac_dat_cd4 = Reduce(function(x, y)merge(x, y, all = T), vac_dat_cd4)
vac_dat_cd4[vac_dat_cd4 == 0] = 1
vac_dat_cd4[is.na(vac_dat_cd4)] = 0
colnames(vac_dat_cd4) = c('Epitope_ID', 'DV_status','DENVax', 'DPIV', 'TV003',
                          'Dengvaxia', 'TDEN', 'TVDV')
vac_dat_cd4 = vac_dat_cd4[c('Epitope_ID', 'DV_status','TV003', 'TDEN', 'DPIV', 
                            'Dengvaxia', 'DENVax', 'TVDV')]
vac_dat_cd4 = split(vac_dat_cd4, f = vac_dat_cd4$DV_status)

cd4_align_metadata = lapply(seq_along(align_info_iedb_cd4), function(i){
  align_dat = align_info_iedb_cd4[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd4_metadata[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd4[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd4, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd4_align_metadata) = c('CD4_DENV1', 'CD4_DENV2', 'CD4_DENV3', 'CD4_DENV4')

experimental_data_cd8 = read.csv('tcell_table_export_1655118742_cd8_9feb.csv')
colnames(experimental_data_cd8) = experimental_data_cd8[1,]
experimental_data_cd8 = experimental_data_cd8[-1,]
experimental_data_cd8 = experimental_data_cd8[c("Epitope ID", "Method/Technique", 
                                                "Assay Group", "Qualitative Measure")]
colnames(experimental_data_cd8)[1] = 'Epitope_ID'
experimental_data_cd8 = experimental_data_cd8[!duplicated(
  experimental_data_cd8$Epitope_ID),]

#add vaccines info
setwd('draft_fig_ver2/supply/')
vac_dat_cd8 = lapply(vac_files, function(x){
  dat = read_xlsx(x, sheet = 'CD8')
  dat = dat[c('Epitope_ID', 'value')]
  dat$value = as.numeric(dat$value)
  colnames(dat)[2] = gsub('_.*', '', x)
  return(dat)
})
setwd('..')
setwd('..')

temp = iedb_cd8[c('Epitope.ID', 'DV_status')]
colnames(temp)[1] = 'Epitope_ID'
vac_dat_cd8 = c(list(temp), vac_dat_cd8)
vac_dat_cd8 = Reduce(function(x, y)merge(x, y, all = T), vac_dat_cd8)
vac_dat_cd8[vac_dat_cd8 == 0] = 1
vac_dat_cd8[is.na(vac_dat_cd8)] = 0
colnames(vac_dat_cd8) = c('Epitope_ID', 'DV_status','DENVax', 'DPIV', 'TV003',
                          'Dengvaxia', 'TDEN', 'TVDV')
vac_dat_cd8 = vac_dat_cd8[c('Epitope_ID', 'DV_status','TV003', 'TDEN', 'DPIV', 
                            'Dengvaxia', 'DENVax', 'TVDV')]
vac_dat_cd8 = split(vac_dat_cd8, f = vac_dat_cd8$DV_status)

cd8_align_metadata = lapply(seq_along(align_info_iedb_cd8), function(i){
  align_dat = align_info_iedb_cd8[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd8_metadata[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd8[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd8, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd8_align_metadata) = c('CD8_DENV1', 'CD8_DENV2', 'CD8_DENV3', 'CD8_DENV4')

###
#strain info
epi_path = list.files('country_wise_data/', pattern = 'json')
countries = str_to_title(gsub('_.*', '', epi_path))

setwd('country_wise_data/')
all_seqs = lapply(epi_path, read.and.filter.strain)
all_seqs = lapply(all_seqs, function(x){
  x['Dengue virus'] = NULL
  out = lapply(x, function(df){
    df = df[!duplicated(df$aaSeq),]
  })
  return(out)
})
setwd('..')
names(all_seqs) = countries
all_seqs = all_seqs[c('Thailand', 'Brazil', 'Mexico')]
all_seqs = lapply(all_seqs, function(x)bind_rows(x))
all_seqs = bind_rows(all_seqs)
all_seqs = all_seqs[intersect(colnames(all_seqs), colnames(all_data_dv))]
india_seq = all_data_dv[intersect(colnames(all_seqs), colnames(all_data_dv))]
all_seqs = bind_rows(india_seq, all_seqs)
all_seqs = all_seqs[with(all_seqs, order(country, virusType)),]
all_seqs$id_year = NULL
all_seqs$date = as.character(all_seqs$date)
all_seqs_india = all_seqs[all_seqs$country == 'India',]

out = c(list('sequence_data' = all_seqs_india),cd4_align_metadata, cd8_align_metadata)
write_xlsx(out, 'draft_fig_ver2/supply/aligmnent_metadata_india.xlsx')

######################################################
#Thailand
align_info_all = readRDS('rds_data/country_wise_alignment_21mar.rds')

setwd('country_wise_data/')
all_seqs = lapply(epi_path, read.and.filter.strain)
all_seqs = lapply(all_seqs, function(x){
  x['Dengue virus'] = NULL
  out = lapply(x, function(df){
    df = df[!duplicated(df$aaSeq),]
  })
  return(out)
})
setwd('..')
names(all_seqs) = countries
all_seqs = all_seqs[c('Thailand', 'Brazil', 'Mexico')]
all_seqs = lapply(all_seqs, function(x)bind_rows(x))

align_info_thailand = align_info_all[['Thailand']]
align_info_thailand_cd4 = align_info_thailand[['CD4_align']]
align_info_thailand_cd8 = align_info_thailand[['cd8_alignment']]

all_seqs_thailand = all_seqs[['Thailand']]
all_seqs_thailand = all_seqs_thailand[c('ncbiId', 'strainname',	'host',
                                        'date',	'continent', 'country',
                                        'protName',	'virusType', 'aaSeq')]
#save alignment metadata
cd4_metadata_thailand = lapply(align_info_thailand_cd4, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd4)
})
cd8_metadata_thailand = lapply(align_info_thailand_cd8, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd8)
})

cd4_align_meta_thailand = lapply(seq_along(align_info_thailand_cd4), function(i){
  align_dat = align_info_thailand_cd4[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd4_metadata_thailand[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd4[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd4, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd4_align_meta_thailand) = c('CD4_DENV1', 'CD4_DENV2', 'CD4_DENV3', 'CD4_DENV4')

cd8_align_meta_thailand = lapply(seq_along(align_info_thailand_cd8), function(i){
  align_dat = align_info_thailand_cd8[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd8_metadata_thailand[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd8[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd8, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd8_align_meta_thailand) = c('CD8_DENV1', 'CD8_DENV2', 'CD8_DENV3', 'CD8_DENV4')

out = c(list('sequence_data' = all_seqs_thailand), cd4_align_meta_thailand, 
        cd8_align_meta_thailand)
write_xlsx(out, 'draft_fig_ver2/supply/aligmnent_metadata_thailand.xlsx')

##############################
###Brazil
align_info_all = readRDS('rds_data/country_wise_alignment_21mar.rds')

align_info_brazil = align_info_all[['Brazil']]
align_info_brazil_cd4 = align_info_brazil[['CD4_align']]
align_info_brazil_cd8 = align_info_brazil[['cd8_alignment']]

all_seqs_brazil = all_seqs[['Brazil']]
all_seqs_brazil = all_seqs_brazil[c('ncbiId', 'strainname', 'host',
                                    'date', 'continent', 'country',
                                    'protName', 'virusType', 'aaSeq')]
#save alignment metadata
cd4_metadata_brazil = lapply(align_info_brazil_cd4, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd4)
})
cd8_metadata_brazil = lapply(align_info_brazil_cd8, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd8)
})

cd4_align_meta_brazil = lapply(seq_along(align_info_brazil_cd4), function(i){
  align_dat = align_info_brazil_cd4[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd4_metadata_brazil[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd4[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd4, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd4_align_meta_brazil) = c('CD4_DENV1', 'CD4_DENV2', 'CD4_DENV3', 'CD4_DENV4')

cd8_align_meta_brazil = lapply(seq_along(align_info_brazil_cd8), function(i){
  align_dat = align_info_brazil_cd8[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd8_metadata_brazil[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd8[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd8, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd8_align_meta_brazil) = c('CD8_DENV1', 'CD8_DENV2', 'CD8_DENV3', 'CD8_DENV4')

out = c(list('sequence_data' = all_seqs_brazil), cd4_align_meta_brazil, 
        cd8_align_meta_brazil)
write_xlsx(out, 'draft_fig_ver2/supply/aligmnent_metadata_brazil.xlsx')

####################
#Mexico
align_info_all = readRDS('rds_data/country_wise_alignment_21mar.rds')

align_info_mexico = align_info_all[['Mexico']]
align_info_mexico_cd4 = align_info_mexico[['CD4_align']]
align_info_mexico_cd8 = align_info_mexico[['cd8_alignment']]

all_seqs_mexico = all_seqs[['Mexico']]
all_seqs_mexico = all_seqs_mexico[c('ncbiId', 'strainname', 'host',
                                    'date', 'continent', 'country',
                                    'protName', 'virusType', 'aaSeq')]
#save alignment metadata
cd4_metadata_mexico = lapply(align_info_mexico_cd4, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd4)
})
cd8_metadata_mexico = lapply(align_info_mexico_cd8, function(x){
  x = epi.metadata(x, epitope_data = iedb_cd8)
})

cd4_align_meta_mexico = lapply(seq_along(align_info_mexico_cd4), function(i){
  align_dat = align_info_mexico_cd4[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd4_metadata_mexico[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd4[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd4, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd4_align_meta_mexico) = c('CD4_DENV1', 'CD4_DENV2', 'CD4_DENV3')

cd8_align_meta_mexico = lapply(seq_along(align_info_mexico_cd8), function(i){
  align_dat = align_info_mexico_cd8[[i]]
  align_dat[align_dat == 0] = 1
  align_dat[is.na(align_dat)] = 0
  align_dat$Epitope_ID = rownames(align_dat)
  meta_dat = cd8_metadata_mexico[[i]]
  meta_dat = meta_dat[c('Epitope_ID', 'aa_seq', 'protein_name', 'start', 'end',
                        'Allele', 'step', 'Serotype')]
  vaccine_dat = vac_dat_cd8[[i]]
  vaccine_dat$DV_status = NULL
  out = Reduce(function(x, y)merge(x, y, by = 'Epitope_ID'), 
               list(meta_dat, experimental_data_cd8, align_dat,
                    vaccine_dat))
  proteins = c('C', 'PreM', 'E', 'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5')
  out = out %>% arrange(factor(out$protein_name, levels = proteins))
  print(paste(nrow(align_dat), nrow(out), nrow(vaccine_dat), sep = ', '))
  return(out)
})
names(cd8_align_meta_mexico) = c('CD8_DENV1', 'CD8_DENV2', 'CD8_DENV3')

out = c(list('sequence_data' = all_seqs_mexico), cd4_align_meta_mexico, 
        cd8_align_meta_mexico)
write_xlsx(out, 'draft_fig_ver2/supply/aligmnent_metadata_mexico.xlsx')
