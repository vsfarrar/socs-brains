#RNAseq dataset
socs_rna <- read.csv("~/Documents/projects/socs-brains/data/socs_cish_RNAseq_hyp_hormone_data.csv")

#select characterization and hypothalamus only 
characterization <- c("bldg", "lay", "inc.d3", "inc.d9","inc.d17","hatch","n5", "n9")

socs_rna <- 
  socs_rna %>%
  filter(stage %in% characterization) %>%
  mutate(stage = factor(stage, levels = characterization)) %>%
  filter(tissue == "hypothalamus")

#correlation plots 

socs_rna %>%
  ggplot(aes(x = log(prolactin_conc), y = vsd_counts)) + 
  geom_point(aes(color = stage, shape = sex)) + 
  geom_smooth(method = "lm", se = TRUE) + 
  #scale_color_manual(values = stage_colors) + 
  facet_wrap(~gene)

#convert to wider for correlation matrix

socs_wide <- 
socs_rna %>%
  select(-X, -X.1) %>% #remove extraneous columns
  pivot_wider(names_from = gene, values_from = vsd_counts) 

#only numerical variables for correlation functions 
socs_corr <- 
  socs_wide %>%
  select(-bird_id, -sex, -tissue, -stage) %>%
  drop_na(prolactin_conc) %>%
  mutate(log_prolactin = log(prolactin_conc))

### CORRELATIONS ### 

cor(socs_corr)

library("Hmisc")
cor_2 <- rcorr(as.matrix(socs_corr))

### MODELS ###

cish <- socs_rna %>% filter(gene == "CISH")

m_cish <- lm(vsd_counts ~ stage + sex, data = cish)

cish_cld <- as.data.frame(cld(emmeans(m_cish, ~ stage),
                              Letters = letters ) ) 


socs2 <- socs_rna %>% filter(gene == "SOCS2")

m_socs2 <- lm(vsd_counts ~ stage + sex, data = socs2)

socs2_cld <- as.data.frame(cld(emmeans(m_socs2, ~ stage),
                              Letters = letters ) ) 


socs_rna %>%
  ggplot(aes(x = stage, y = vsd_counts)) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar" ) + 
  #scale_color_manual(values = stage_colors) + 
  facet_wrap(~gene)



stage_colors <- c("bldg" = "#C0C0C0",
                  "hatch" = "	#800000",
                  "inc.d17" = "#FF0000",
                  "inc.d3" = "#008080",
                  "inc.d9" = "	#FFFF00",
                  "lay"  = "#1E90FF", 
                  "n5" = "#EE82EE",
                  "n9" = "#FF00FF")
