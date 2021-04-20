#import RNAseq full dataset 
  #large, will take a moment to load
rnaseq <- read.csv("~/Documents/projects/parental_care_RNAseq/full_RNAseq_dataset/parentalcare_rnaseq_allgenes_vsc_RMH.csv")

#select CISH SOCS genes
cish_socs <-
rnaseq %>%
  filter(gene %in% c("CISH", "PRLR", "SOCS1", "SOCS2", "SOCS3"))

#hormones data for characterization only 
hormones  <- read.csv("~/Documents/projects/parental_care_RNAseq/archived_raw_data/RMH_hormonesandPRL.csv")

#convert prolactin plasma data from Rayna to joinable format 
prolactin_plasma_rnaseq_format <- 
  hormones %>%
  select(bird_id,sex,tissue,stage = treatment, hormone, plasma_conc) %>%
  filter(hormone == "prolactin") %>% #select only prolactin
  select(-hormone, - tissue) #no longer need the "hormone" column 


#join all tissue data with left_join 

all_tissues <- 
  prolactin_plasma_rnaseq_format %>%
  select(-stage) %>%  #remove stage from this original file so there are not duplicated stages
  left_join(cish_socs, ., by = c("bird_id", "sex")) %>% 
  rename( prolactin_conc = plasma_conc)  %>%   #rename plasma conc to prolactin conc
  distinct_all() #avoids duplication issue, only keep unique entire rows. 


#trim and clean up factors and stages

all_tissues$stage <- str_remove(all_tissues$stage, ".NYNO") #remove extra .NYNO from some stages
all_tissues$stage <- as.factor(all_tissues$stage) #convert to factor
all_tissues$stage <- recode(all_tissues$stage, "manip.d8" = "m.inc.d8")  #recode manip.d8 to m.inc.d8

levels(all_tissues$stage) 

#select characterizations only 
#select hypothalamus only
characterization <- c("bldg", "lay", "inc.d3", "inc.d9","inc.d17","hatch","n5", "n9")

cish_char <- 
  all_tissues %>%
  filter(stage %in% characterization) %>%
  mutate(stage = factor(stage, levels = characterization)) %>%
  filter(tissue == "hypothalamus")

levels(cish_char$stage)


#EXPORT and save as .csv
write.csv(cish_char, "~/Documents/projects/socs-brains/data/socs_cish_RNAseq_hyp_hormone_data.csv")
