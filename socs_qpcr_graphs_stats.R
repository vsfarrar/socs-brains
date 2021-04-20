library(EnvStats)
library(tidyverse)

setwd("~/Documents/projects/socs-brains/")
socs <- read.csv("data/2021-04-06_socs_cish_qpcr_ddct.csv")

#set vehicle to the reference level 
socs$treatment <- relevel(as.factor(socs$treatment), ref = "vehicle")

#create dummy variable for sex*treatment 

socs$treatsex <- paste(socs$treatment, socs$sex, sep = "_")
socs$treatsex <- factor(socs$treatsex, levels = c("vehicle_m", "vehicle_f", "PRL_m", "PRL_f"))


#by treatment only 
socs %>%
  filter(gene == "socs2") %>%
  ggplot(aes(x = treatment, y = log_fold, fill = treatment, label = sample)) + 
  geom_point(aes(shape = sex), position = position_jitter(0.1), size = 2, alpha =0.7,  fill = "white", color = "black") + 
  stat_summary(fun.y = mean, geom = "col", color = "black", alpha = 0.6) + 
  stat_summary(fun.data = "mean_se", color = "black",geom = "errorbar", width = 0.5) + 
  scale_shape_manual(values= c(21,24)) + 
  scale_fill_manual(values = c("white", "black")) + 
  geom_hline(yintercept = 0) + 
  stat_n_text() + 
  labs(title = "SOCS2", x = NULL, y = "log(fold_change)") + 
  theme_classic() + 
  facet_wrap(~tissue) + 
  theme(legend.position = "none")


##STATISTICS

socs %>% filter(gene == "cish" ) %>%
  t.test(.$log_fold ~ .$treatment, data = .)


socs %>% filter(gene == "socs3") %>%
  aov(log_fold ~ treatment + tissue + sex, data = .) %>%
  car::Anova(., type = "III")

#pivot wider for correlations. 

socs %>%
  select(treatment:gene,log_fold) %>%
  pivot_wider(names_from = gene,
              values_from = log_fold)

#ALTERNATIVE GRAPH

#by sex*treatment
socs %>%
  filter(gene == "cish") %>%
  ggplot(aes(x = treatsex, y = log_fold, fill = treatment)) + 
  geom_point(aes(shape = sex), position = position_jitter(0.05), size = 2,  color = "gray") + 
  stat_summary(fun.y = mean, geom = "col", color = "black", alpha = 0.5) + 
  stat_summary(fun.data = "mean_se", color = "black",geom = "errorbar", width = 0.4) + 
  scale_fill_manual(values = c("white", "black")) + 
  geom_hline(yintercept = 0) + 
  stat_n_text() + 
  theme_classic()