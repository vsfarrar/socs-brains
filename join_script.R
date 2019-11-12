###This is the script used for the original join of prolactin data and CISH/SOCS RNAseq data
#import data
#"hyp" for hypothalamus gene counts. Data from Rayna from transcriptome parental care study.
hyp<- read.csv(file="~/Desktop/data/10_hypcandidatecounts.csv")

#"prl" = prolactin RIA data from Fred Angelier for parental care study birds
prl<-read.csv(file="~/Desktop/data/prolactin_ria_data.csv")

#code from PRLR_MS_full_analyses.Rmd  (HPC crops project) for joins. Used to create a join that I then filled in manually. 

#clean data using stringr

#create a clean_id variable for band combinations in  RNAseq dataset 
hyp$clean_id<-str_replace_all(hyp$bird, "x", "") #remove meaningless "x" from band combos
hyp$clean_id<-str_replace_all(hyp$clean_id, "[[:punct:]]", "") #remove punctuation

#create clean_id variable in prolactin RIA dataset
prl$clean_id<-str_replace_all(prl$ColorBands, "x", "")
prl$clean_id<-str_replace_all(prl$clean_id,"[[:punct:]]", "" )


#join prolactin plasma data with prl.long RNAseq data
prl.join<-left_join(hyp, prl, by ="clean_id")

#filter the join by those that actually joined 

socs_join<- prl.join %>%
  filter(!is.na(Prolactin.ng.mL))

#sample size: 
socs_join %>%
  group_by(treatment,sex) %>%
  summarize(n=n())