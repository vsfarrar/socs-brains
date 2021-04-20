### qPCR DATA PROCESSING - FUNCTIONIZED SCRIPT
#created 12/5/2020

#imported data should be in long format with the following lowercase columns: 
  #sample,gene,tissue,cq

###################### EDIT HERE ############################
#import data
#replace abc in whole script with your tissue / dataframe name 

setwd("~/Documents/projects/socs-brains/")
hyp <- read.csv("data/2020-12-12_PRLvsHPG_qPCR_all_PVN_POA_raw_data.csv", stringsAsFactors = F)
key <-read.csv("~/Documents/projects/PRL_vs_HPG/data/PRLvsHPG_sample_key.csv")

#source packages, functions
source("~/Documents/GitHub/PRL-vs-HPG/qpcr/prlvshpg_qpcr_packages.R")
source("~/Documents/GitHub/PRL-vs-HPG/qpcr/prlvshpg_qpcr_functions.R")

#if you want to normalize to specific reference genes, list them here
custom_refgenes <- c("hprt1", "gapdh")
current_date <-Sys.Date()
######################################################################


#clean and join data 
hyp <- 
  hyp %>%
  mutate(gene = tolower(gene), #all lowercase genes
         sample = as.numeric(sample) #sample to numeric
  ) %>%
  replace_with_na(replace = list(gene = "")) %>% #blank to NA
  drop_na(gene) 

#join with project key 
#adds treatment group, sex, year etc. 

hyp_join <- hyp %>%
  inner_join(., key, by = c("sample"="band_no")) %>%
  mutate(treatment = factor(treatment, levels = c("vehicle", "PRL")))

#clean triplicates
hyp_clean <- clean_triplicates(hyp_join)

#average cq across triplicates (returns 1 avg cq value for each sample:gene)
hyp_avg <- avg_triplicates(hyp_clean)


#get average reference gene and dct
hyp_dct <- 
  hyp_avg %>%
  group_by(tissue,sample) %>%
  mutate(is_ref = ifelse(gene %in% reference_genes, 1, 0),
         ref_gene = mean(mean_cq[gene %in% custom_refgenes], na.rm = T)) %>% #get average reference gene
  mutate(dct = mean_cq - ref_gene) #calculate dct

#get dataframe of controls 
hyp_controls <- 
  hyp_dct %>%
  filter(is_ref == 0) %>% #do not calculate control dct for reference genes
  group_by(tissue,gene) %>%
  summarise(control_dct = mean(dct[treatment == "vehicle"], na.rm = T)) # control = average dct for vehicle group

#join with control data by genes
hyp_dct <- left_join(hyp_dct, hyp_controls) #joins by TISSUE and gene

# FINAL PRODUCT: calculate ddct, fold change, and log fold change
hyp_qpcr <- 
  hyp_dct %>%
  group_by(tissue,sample, gene) %>%
  mutate(ddct = dct-control_dct, 
         fold_change = 2^-(ddct),
         log_fold = log2(fold_change)) %>%
  ungroup()

#####################################################

#filter for SOCS-CISH project 
socs_qpcr <-
hyp_qpcr %>%
  filter(gene %in% c("cish","socs1","socs2","socs3", "prlr"))

#optional: save processed data file 
setwd("~/Documents/projects/socs-brains/data/")
write.csv(socs_qpcr, paste0(current_date,"_socs_cish_qpcr_ddct.csv"), row.names = F)
