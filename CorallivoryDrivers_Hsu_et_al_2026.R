#' ---
#' title: "Hsu et al. (2026): Corallivore assemblage, feeding behaviour, and cover of Acropora corals drive fish corallivory on tropical reefs"
#' author:
#'   - name: Tsai-Hsuan Tony Hsu
#'     orcid: 0000-0002-9362-9769
#' date: "`r format(Sys.Date(), '%d %B %Y')`"
#' format: 
#'   html:
#'     toc: true
#'     toc-depth: 4
#'     number-sections: true
#'     toc-location: left
#'     embed-resources: true
#' editor: visual
#' ---
#' 
#' ## **Description**
#' 
#' Quarto document to replicate analyses in Hsu et al. (2026). If used in full or in part, please cite the data and the original publication:
#' 
#' **Publication**: Hsu T-H.T, Hoey A.S, Ferrari R, Toor M, Gordon S, Remmers T, Smallhorn-West J, Piccaluga A, Figueira W.F (in press) Corallivore assemblage, feeding behaviour, and *Acropora* corals drive fish corallivory on tropical reefs. *Ecosphere*. [doi](doi)
#' 
#' **Data**: Hsu T-H.T, Hoey A.S, Ferrari R, Toor M, Gordon S, Remmers T, Smallhorn-West J, Piccaluga A, Figueira W.F (2026) Data from "Corallivore assemblage, feeding behaviour, and *Acropora* corals drive fish corallivory on tropical reefs.". \[Data set\]. Sydney eScholarship, [doi](doi)
#' 
#' Data and **R script** also available on [GitHub](https://github.com/THTonyHsu/Drivers-of-corallivory.git)
#' 
#' ## **Required packages**
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# Data wrangling and plotting
library(tidyverse)
library(cowplot)
library(ggrepel)
library(qgraph) # network plot
# Data analyses
library(lmerTest) # LMM
library(emmeans) # Tukey's HSD test
library(vegan) # nMDS and PERMANOVA
library(pairwiseAdonis) # pairwise PERMANOVA
library(car) # vif
library(glmulti) # glm with AICc
library(ggeffects) # model prediction
library(piecewiseSEM) # SEM

#' 
#' ## **Dataset**
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false

# Corallivory data
df_cor = read.csv("Data/Hsuetal_dataset_Corallivory.csv")
df_cor$Reef = factor(df_cor$Reef, levels = c("CBHE", "OCDA", "ONLI", "TSMA"))

# List of obligate corallivores (as defined by Cole et al., 2008)
ObCor = read.csv("Data/Hsuetal_dataset_Obligate_corallivore.csv")

# Benthic communities: 72 taxa
Ben_cov = read.csv("Data/Hsuetal_dataset_benthic_taxa.csv")

# Taxa to morpho-functional groups (MFGs)
Taxa_MFG = read.csv("Data/Hsuetal_dataset_benthic_MFG.csv")

# Reference list for classifying prey items into MFG
Ben_food = read.csv("Data/Hsuetal_dataset_prey_MFG.csv")
colnames(Ben_food)[1] = "Comment"

# Complexity metrics
ComplexP = read.csv("Data/Hsuetal_dataset_Complexity.csv")

# Fish assemblages from UVC
Fishcomm = read.csv("Data/Hsuetal_dataset_FishUVC.csv")

#' 
#' ## **Data wrangling**
#' 
#' ### Create a summary table
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 150)
CF_smry = df_cor[, 1:5] %>% unique()
row.names(CF_smry) = NULL
CF_smry$Habitat = as.factor(CF_smry$Habitat)
CF_smry$Reef = factor(CF_smry$Reef, levels = levels(df_cor$Reef))
CF_smry$Site = factor(CF_smry$Site)
CF_smry$Region = ifelse(CF_smry$Reef == "TSMA", "TS", "GBR")
CF_smry = CF_smry[,c(6,1:5)]

# Check the no. of plots at each site
table(CF_smry$Site) # No corallivory activities at CBHE_FR1S_P1 & OCDA_FR1S_P2

# Add empty rows for the plots with zero feeding events
CF_smry = rbind(CF_smry[1:8,], CF_smry[9,], CF_smry[9:24,], CF_smry[24:62,])
row.names(CF_smry) = NULL
CF_smry[9,6] = 1
CF_smry[26,6] = 2

# Add transect info
CF_smry$TR = ifelse(CF_smry$Plot<3, "T1", "T2")

#' 
#' ### Corallivory
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false

# Calculate the event duration per event
df_event = df_cor %>% filter(Activity == "Passing")  %>% 
  pivot_wider(id_cols = c(Reef, Site, Plot, SN, TL, ID), 
              names_from = Comment, values_from = Time) %>% 
  mutate(Event_duration = `time out`-`time in`)

# Calculate bites and bite rate per event
df_event = df_cor %>% filter(Activity == "Feeding") %>% 
  group_by(Reef, Site, Plot, SN, TL, ID) %>%
  summarise(Bite = n()) %>% left_join(df_event, .) %>% na.omit() %>% 
  mutate(Bite_rate = Bite/Event_duration)

# Assemblage-level corallivory (Corallivorous bites by species across plots)
CF_sp = df_cor %>% filter(Activity == "Feeding") %>% left_join(., Ben_food) %>% na.omit() %>% 
  group_by(Reef, Site, Habitat, Position, Plot, SN) %>% 
  summarise(Bite = n()) %>% 
  pivot_wider(id_cols = c(Reef, Site, Habitat, Position, Plot), 
              names_from = SN, values_from = Bite, values_fill = 0)
CF_sp = cbind(data.frame(Region = ifelse(CF_sp$Reef == "TSMA", "TS", "GBR")), CF_sp)

#' 
#' ### Fish assemblages
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
Fishcomm = Fishcomm[is.na(Fishcomm$Biomass_g_per100m2)*-1,]
## Fish trophic groups
FishTG = Fishcomm %>% group_by(Reef, Habitat, Position, Site, TR, Trophic_group) %>% 
  summarise(Biomass = sum(Biomass_g_per100m2), Density = sum(Density_per100m2))

#' 
#' ### Benthic communities
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
Ben_MFG = pivot_longer(Ben_cov, cols = Acropora.Arborescent.Staghorn:Zoanthid, names_to = "taxa") %>% 
  left_join(., Taxa_MFG) %>% group_by(Region, Reef, Site, Habitat, Position, Plot, MFG) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(id_cols = c(Region, Reef, Site, Habitat, Position, Plot), 
              names_from = MFG, values_from = value, values_fill = 0) %>% 
  select(-c(oth.mob, oth.ben, Abiotic, UnID)) # remove other taxa, unidentified, and abiotoc

# Only keep abundant MFG with a max >5%
MFG15 = colnames(Ben_MFG[,-(1:6)])[apply(Ben_MFG[,-(1:6)], 2, max)>5] # 15 morpho-functional groups
Ben_com = cbind(Ben_MFG[,1:6], Ben_MFG[,MFG15])

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
# Find the most common 6 coral preys which were later used as predictors of corallivory
df_cor %>% filter(Activity == "Feeding") %>% pull(Comment) %>% 
  table() %>% sort(., decreasing = T) %>% .[1:6]

df_cor %>% filter(Activity == "Feeding") %>% pull(Comment) %>% 
  table() %>% sort(., decreasing = T) %>% .[1:6] %>% sum()/4329*100 # ~87%

#' 
#' ## **Metrics for analyses**
#' 
#' ### Total corallivory (no. of bites per 30 mins)
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry$Total_corallivory = df_event %>% group_by(Reef, Site, Plot) %>% summarise(Total_corallivory = sum(Bite)) %>% 
  full_join(CF_smry, .) %>% pull(Total_corallivory) %>% replace_na(., 0)

#' 
#' ### Total events
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry$Total_events = df_event %>% group_by(Reef, Site, Plot) %>% summarise(Total_events = n()) %>% 
  full_join(CF_smry, .) %>% pull(Total_events) %>% replace_na(., 0)

#' 
#' ### Feeding intensity (bite rate per event)
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry$Feeding_intensity = df_event %>% group_by(Reef, Site, Plot) %>% summarise(Feeding_intensity = mean(Bite_rate)) %>% 
  full_join(CF_smry, .) %>% pull(Feeding_intensity) %>% replace_na(., 0)

#' 
#' ### Event duration
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry$Event_duration = df_event %>% group_by(Reef, Site, Plot) %>% summarise(Event_duration = mean(Event_duration)) %>% 
  full_join(CF_smry, .) %>% pull(Event_duration) %>% replace_na(., 0)

#' 
#' ### Richness of corallivores on the videos
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry$Richness_Corallivores = df_cor[,c("Reef", "Site", "Plot", "SN")] %>% unique() %>% 
  group_by(Reef, Site, Plot) %>% summarise(Richness_Corallivores = n()) %>% 
  left_join(CF_smry, .) %>% pull(Richness_Corallivores) %>% replace_na(., 0)


#' 
#' ### Biomass of corallivores (100 g per 100 m-2)
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false

CF_smry$Biomass_Corallivores = Fishcomm %>% filter(Valid_as %in% ObCor$Valid_as) %>% 
  group_by(Reef, Site, TR, Trophic_group) %>% 
  summarise(Biomass = sum(Biomass_g_per100m2)) %>% 
  pivot_wider(id_cols = c(Reef, Site, TR), 
              names_from = Trophic_group, values_from = Biomass, values_fill = 0) %>% 
  mutate(Biomass_Corallivores = Corallivores) %>% left_join(CF_smry[,1:7], .) %>% 
  pull(Biomass_Corallivores) %>% replace_na(., 0)/100

#' 
#' ### Biomass of piscivores (kg per 100 m-2)
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry$Biomass_Piscivores = FishTG %>% filter(Trophic_group == "Piscivores") %>% 
  pivot_wider(id_cols = c(Reef, Site, TR), 
              names_from = Trophic_group, values_from = Biomass, values_fill = 0) %>% 
  left_join(CF_smry[,1:7], .) %>% pull(Piscivores) %>% replace_na(., 0)/1000


#' 
#' ### Complexity: plane rugosity (plot level)
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry$Rugosity = left_join(CF_smry, ComplexP) %>% arrange(Reef, Position) %>% pull(Rugosity_plane)

#' 
#' ### Top 6 abundant coral prey
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry = CF_smry %>% mutate(Acropora.arb = Ben_com$Acro.arb,
                             Acropora.cory = Ben_com$Acro.cory,
                             Acropora.tab = Ben_com$Acro.tab,
                             Porites.mas = Ben_com$Porites.mas,
                             Montipora = Ben_com$Montipora,
                             Pocillopora = Ben_com$Pocillopora)


#' 
#' ### Scale summary ts
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# All variables scaled
CF_smry_s = CF_smry
CF_smry_s[,-c(1:7)] = scale(CF_smry[,-c(1:7)])

# Separate GBR and TS data
CF_smry_GBR = CF_smry[1:48,]
CF_smry_TS = CF_smry[49:64,]

# Normalize variables
CF_smry_GBR_s = cbind(CF_smry_GBR[,1:7], scale(CF_smry_GBR[,-c(1:7)]))
CF_smry_TS_s = cbind(CF_smry_TS[,1:7], scale(CF_smry_TS[,-c(1:7)]))

#' 
#' ## **Data analysis**
#' 
#' ### Summary of corallivory
#' 
#' #### Total corallivory
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
sum(CF_smry$Total_corallivory) # 4329 bites

#' 
#' #### Total events
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
sum(CF_smry$Total_events) # 842 events

#' 
#' #### Feeding intensity
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
mean(CF_smry$Feeding_intensity) # Mean Feeding intensity: 12.31
sd(CF_smry$Feeding_intensity)/sqrt(nrow(CF_smry)) # SE of Feeding intensity: 0.55

#' 
#' #### Event duration
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
mean(CF_smry$Event_duration) # Mean event duration: 0.40 (min)
sd(CF_smry$Event_duration)/sqrt(nrow(CF_smry)) # SE of event duration: 0.03

#' 
#' #### Species richness
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
unique(df_cor$SN) # 12 species

#' 
#' ### Spatial variations in total corallivory, corallivore behaviour, and corallivore assemblage
#' 
#' #### Total corallivory
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
mod = lmer(sqrt(Total_corallivory) ~ Reef * Habitat + (1 | Reef:Site), data = CF_smry)
summary(mod)
anova(mod) # Habitat**
shapiro.test(resid(mod))

#' 
#' #### Feeding intensity
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
mod = lmer(Feeding_intensity ~ Reef * Habitat + (1 | Reef:Site), data = CF_smry)
summary(mod)
anova(mod) # N.S.
shapiro.test(resid(mod))

#' 
#' #### Event duration
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
mod = lmer(Event_duration ~ Reef * Habitat + (1 | Reef:Site), data = CF_smry)
summary(mod)
anova(mod) # N.S.
shapiro.test(resid(mod))

#' 
#' #### Richness of corallivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
mod = lmer(Richness_Corallivores ~ Reef * Habitat + (1 | Reef:Site), data = CF_smry)
summary(mod)
anova(mod) # N.S.
shapiro.test(resid(mod))

#' 
#' #### Biomass of corallivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
mod = lm(Biomass_Corallivores ~ Reef*Habitat, CF_smry[CF_smry$Plot %in% c(1,3),])
summary(mod)
anova(mod) # Reef*, Habiat*
shapiro.test(resid(mod))
emmeans(mod, list(pairwise ~ Reef), adjust = "tukey")

#' 
#' ### Spatial variation in assemblage-level corallivory and benthic communities
#' 
#' #### Assemblage-level corallivory
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# nMDS of Assemblage-level corallivory (corallivorous bites by species)]
nmds_bite_sp = metaMDS(as.data.frame(CF_sp[,-c(1:6)]), distance='bray',trymax=999)

# extract data from nMDS for plotting
df_nmds_bite_sp = as.data.frame(scores(nmds_bite_sp)$sites) # extract nMDS value of the 60 plots (another 4 plots were absent due to no feeding activity)
df_nmds_bite_sp$Reef = CF_sp$Reef
df_nmds_bite_sp$Habitat = CF_sp$Habitat
df_nmds_bite_sp$Region = CF_sp$Region
sps_bite_sp = as.data.frame(scores(nmds_bite_sp,"species")) # extract nMDS value of all species
sps_bite_sp$SN = rownames(sps_bite_sp)
sps_bite_sp = left_join(sps_bite_sp,ObCor[,4:5])

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
# PERMANOVA on assemblage-level corallivory
PERM_C = adonis2(CF_sp[,-c(1:6)]~ Reef*Habitat, data = CF_sp, 
                 permutations = 999, method="bray", by = "terms")
PERM_C

# pairwise PERMANOVA on assemblage-level corallivory
CF_sp_RxH = CF_sp
CF_sp_RxH$Reef_Habitat = interaction(CF_sp_RxH$Reef, CF_sp_RxH$Habitat)
CF_sp_RxH = CF_sp_RxH[, c(18, 1:17)]
pairwise.adonis2(CF_sp_RxH[,-c(1:7)] ~ Reef_Habitat, 
                 method="bray", p.adjust.method = "bonferroni", nperm = 999, data = CF_sp_RxH)
# Higher bites at back than front sites only at Heron

# Similarity Percentages analysis on assemblage-level corallivory
simper_c = simper(comm = CF_sp_RxH[,-(1:7)], group = CF_sp_RxH$Reef_Habitat) %>% 
  summary()
# Heron: Back vs Fore
simper_c$CBHE.BA_CBHE.FR # C. rainfordi* and C. trifascialis**

#' 
#' #### Benthic communities
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# nMDS of benthic community
nmds_ben = metaMDS(Ben_com[,-1:-6], distance='bray',trymax=999)

# extract data from nMDS for plotting
df_nmds_ben = as.data.frame(scores(nmds_ben)$sites) # extract nMDS value of all 64 plots
df_nmds_ben$Reef = Ben_com$Reef
df_nmds_ben$Reef = factor(df_nmds_ben$Reef, levels = levels(df_cor$Reef))
df_nmds_ben$Habitat = Ben_com$Habitat
df_nmds_ben$Region = Ben_com$Region
sps_ben = as.data.frame(scores(nmds_ben,"species")) # extract nMDS value of all benthos
sps_ben$Benthos = rownames(sps_ben)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
# PERMANOVA on benthic communities
PERM_B = adonis2(Ben_com[,-(1:6)]~ Reef*Habitat, data = Ben_com, 
                 permutations = 999, method="bray", by = "terms")
PERM_B

## pairwise PERMANOVA on benthic communities
Ben_com_RxH = Ben_com
Ben_com_RxH$Reef_Habitat = interaction(Ben_com_RxH$Reef, Ben_com_RxH$Habitat)
Ben_com_RxH = Ben_com_RxH[, c(22, 1:21)]
pairwise.adonis2(Ben_com_RxH[,-c(1:7)] ~ Reef_Habitat, 
                 method="bray", p.adjust.method = "bonferroni", nperm = 999, data = Ben_com_RxH)
# All significantly differ at back than front sites

# Similarity Percentages analysis on benthic communities
simper_b = simper(comm = Ben_com_RxH[,-(1:7)], group = Ben_com_RxH$Reef_Habitat) %>% 
  summary()
# Heron: Back vs Fore
simper_b$CBHE.BA_CBHE.FR
# Davies: Back vs Fore
simper_b$OCDA.BA_OCDA.FR
# Lizard: Back vs Fore
simper_b$ONLI.BA_ONLI.FR
# Masig: Back vs Fore
simper_b$TSMA.BA_TSMA.FR

#' 
#' #### Table S3: Benthic composition between back- and fore-reef habitats among reefs
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# Function to compute mean and SE together
calc_mean_se = function(x) {
  mean_val = round(mean(x), 2)
  se_val = round(sd(x) / sqrt(length(x)), 2)
  paste0(mean_val, " \u00B1 ", se_val)
}

TableS3 = Ben_com %>%
  group_by(Region, Reef, Habitat) %>%
  summarise(
    across(c(Acro.arb:alg.ca), .fns = calc_mean_se)) %>%
  as.data.frame()

# Transpose and adjust column names
TableS3 = t(TableS3[,-(1:3)]) %>% as.data.frame()

TableS3 = rbind(
t(data.frame(
  Region = rep(c("GBR", "TS"), times = c(6,2)),
  Reef = rep(c("Heron", "Davies", "Lizard", "Masig"), each = 2),
  Habitat = rep(c("Back", "Fore"), 4)
)),
TableS3) %>% rownames_to_column(var = "MFG")
TableS3 = TableS3[-1,]
require(DT)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
datatable(TableS3, options = list(scrollX = TRUE))

#' 
#' #### Table S4: Benthic composition among GBR and TS sites
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
#| warning: false
# Compute mean ± SE in one step
Ben_R = Ben_com %>%
  group_by(Region) %>%
  summarise_at(vars(Acro.arb:Turf), calc_mean_se) %>%
  as.data.frame()

# Transpose and adjust column names
Ben_R = t(Ben_R[,-1]) %>% as.data.frame()
colnames(Ben_R) = c("GBR", "TS")
Ben_R$Benthos = row.names(Ben_R)
Ben_R = Ben_R %>% arrange(Benthos)

# Compare MFG between GBR and TS using Wilcoxon test
Ben_R_Wt = apply(Ben_com[,-c(1:6)], MARGIN = 2, FUN = function(x)(wilcox.test(x ~ Ben_com$Region))) %>% 
  lapply(., FUN = function(x)(x[["p.value"]])) %>% unlist() %>% as.data.frame()
Ben_R_Wt$Benthos = rownames(Ben_R_Wt)
colnames(Ben_R_Wt)[1] = "Significance"

TableS4 = left_join(Ben_R, Ben_R_Wt)
rownames(TableS4) = TableS4$Benthos
TableS4 = TableS4[,-3]
TableS4$Significance = ifelse(TableS4$Significance>0.05, "NS", 
                              ifelse(TableS4$Significance>0.01, "*",
                                     ifelse(TableS4$Significance>0.001, "**", "***")))

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
datatable(TableS4, options = list(scrollX = TRUE))

#' 
#' ### Drivers of corallivory within regions
#' 
#' #### Optimizing LMs among GBR sites
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry_GBR_t = CF_smry_GBR_s
CF_smry_GBR_t$Total_corallivory = as.vector(scale((CF_smry_GBR$Total_corallivory)^(1/3)))
CF_smry_GBR_t$Biomass_Piscivores = as.vector(scale((CF_smry_GBR$Biomass_Piscivores)^(1/4)))

#' 
#' #### GBR-LM1: total corallivory
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Total_corallivory ~  Feeding_intensity + Event_duration + 
     Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + 
     Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas+ Rugosity, CF_smry_GBR_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# LM model selection using AICc
model_selection = glmulti(Total_corallivory ~ Feeding_intensity + Event_duration + Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora +Porites.mas+ Rugosity,
                          data = CF_smry_GBR_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel
lm_total_corallivory_GBR = lm(Total_corallivory ~ Feeding_intensity + Event_duration + Richness_Corallivores + Acropora.tab, CF_smry_GBR_t)
summary(lm_total_corallivory_GBR)

# Test the normality of model residual
shapiro.test(resid(lm_total_corallivory_GBR))
# Q-Q plot
par(mfrow=c(1,2))
qqnorm(resid(lm_total_corallivory_GBR))
qqline(resid(lm_total_corallivory_GBR), col = "red")
hist(resid(lm_total_corallivory_GBR), breaks = 20, main = NULL)

#' 
#' #### GBR-LM2: feeding intensity
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Feeding_intensity ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora +Porites.mas + Rugosity, CF_smry_GBR_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Feeding_intensity ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas + Rugosity,
                          data = CF_smry_GBR_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel
lm_Feeding_intensity_GBR = lm(Feeding_intensity ~ Acropora.cory + Rugosity, CF_smry_GBR_t)
summary(lm_Feeding_intensity_GBR)
shapiro.test(resid(lm_Feeding_intensity_GBR))
par(mfrow=c(1,2))
qqnorm(resid(lm_Feeding_intensity_GBR))
qqline(resid(lm_Feeding_intensity_GBR), col = "red")
hist(resid(lm_Feeding_intensity_GBR), breaks = 20, main = NULL)

#' 
#' #### GBR-LM3: event duration
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Event_duration ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas + Rugosity, CF_smry_GBR_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Event_duration ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas + Rugosity,
                          data = CF_smry_GBR_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel
lm_Event_duration_GBR = lm(Event_duration ~ Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora, CF_smry_GBR_t)
summary(lm_Event_duration_GBR)
shapiro.test(resid(lm_Event_duration_GBR))
par(mfrow=c(1,2))
qqnorm(resid(lm_Event_duration_GBR))
qqline(resid(lm_Event_duration_GBR), col = "red")
hist(resid(lm_Event_duration_GBR), breaks = 20, main = NULL)

#' 
#' #### GBR-LM4: richness of corallivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Richness_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora +Porites.mas + Rugosity, CF_smry_GBR_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Richness_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas + Rugosity,
                          data = CF_smry_GBR_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Richness_Corallivores_GBR = lm(Richness_Corallivores ~ Acropora.cory + Acropora.arb + Montipora, CF_smry_GBR_t)
summary(lm_Richness_Corallivores_GBR)

shapiro.test(resid(lm_Richness_Corallivores_GBR))
par(mfrow=c(1,2))
qqnorm(resid(lm_Richness_Corallivores_GBR))
qqline(resid(lm_Richness_Corallivores_GBR), col = "red")
hist(resid(lm_Richness_Corallivores_GBR), breaks = 20, main = NULL)

#' 
#' #### GBR-LM5: biomass of corallivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Biomass_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas + Rugosity, CF_smry_GBR_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Biomass_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas+ Rugosity,
                          data = CF_smry_GBR_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Biomass_Corallivores_GBR = lm(Biomass_Corallivores ~ Acropora.tab + Acropora.arb, CF_smry_GBR_t)
summary(lm_Biomass_Corallivores_GBR)
shapiro.test(resid(lm_Biomass_Corallivores_GBR))
par(mfrow=c(1,2))
qqnorm(resid(lm_Biomass_Corallivores_GBR))
qqline(resid(lm_Biomass_Corallivores_GBR), col = "red")
hist(resid(lm_Biomass_Corallivores_GBR), breaks = 20, main = NULL)

#' 
#' #### GBR-LM6: biomass of piscivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Biomass_Piscivores ~ Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora + Porites.mas + Rugosity, CF_smry_GBR_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Biomass_Piscivores ~ Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora +Porites.mas+
                            Rugosity,
                          data = CF_smry_GBR_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Biomass_Piscivores_GBR = lm(Biomass_Piscivores ~ Acropora.tab + Acropora.arb + Montipora + Pocillopora, CF_smry_GBR_t)
summary(lm_Biomass_Piscivores_GBR)
shapiro.test(resid(lm_Biomass_Piscivores_GBR))

par(mfrow=c(1,2))
qqnorm(resid(lm_Biomass_Piscivores_GBR))
qqline(resid(lm_Biomass_Piscivores_GBR), col = "red")
hist(resid(lm_Biomass_Piscivores_GBR), breaks = 20, main = NULL)

#' 
#' #### GBR-LM7: rugosity
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Rugosity ~ Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora +Porites.mas,CF_smry_GBR_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Rugosity ~ Acropora.tab + Acropora.cory + Acropora.arb + 
                            Montipora + Pocillopora +Porites.mas,
                          data = CF_smry_GBR_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_rugosity_GBR = lm(Rugosity ~ Acropora.tab + Montipora + Porites.mas, CF_smry_GBR_t)
lm_rugosity_GBR %>% summary()
shapiro.test(resid(lm_rugosity_GBR))
par(mfrow=c(1,2))
qqnorm(resid(lm_rugosity_GBR))
qqline(resid(lm_rugosity_GBR), col = "red")
hist(resid(lm_rugosity_GBR), breaks = 20, main = NULL)

#' 
#' ### Optimizing LMs among TS sites
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
CF_smry_TS_t = CF_smry_TS_s
CF_smry_TS_t$Total_corallivory = as.vector(scale((CF_smry_TS$Total_corallivory)^(1/3)))
CF_smry_TS_t$Event_duration = as.vector(scale((CF_smry_TS$Event_duration)^(1/4)))

#' 
#' #### TS-LM1: total corallivory
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Total_corallivory ~  Feeding_intensity + Event_duration + 
     Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + 
     Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora +Porites.mas+ Rugosity, CF_smry_TS_t) %>% vif()
# Remove Pocillopora
lm(Total_corallivory ~  Feeding_intensity + Event_duration + 
     Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + 
     Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Porites.mas+
     Rugosity, CF_smry_TS_t) %>% vif()
# Remove tabular Acropora
lm(Total_corallivory ~  Feeding_intensity + Event_duration + 
     Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + 
     Acropora.cory + Acropora.arb + Montipora + Porites.mas+
     Rugosity, CF_smry_TS_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Total_corallivory ~  Feeding_intensity + Event_duration + Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.cory + Acropora.arb + Montipora + Porites.mas+ Rugosity, 
                          data = CF_smry_TS_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_total_corallivory_TS = lm(Total_corallivory ~ Montipora, CF_smry_TS_t)
lm_total_corallivory_TS %>% summary()
shapiro.test(resid(lm_total_corallivory_TS))
par(mfrow=c(1,2))
qqnorm(resid(lm_total_corallivory_TS))
qqline(resid(lm_total_corallivory_TS), col = "red")
hist(resid(lm_total_corallivory_TS), breaks = 20, main = NULL)

#' 
#' #### TS-LM2: feeding intensity
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Feeding_intensity ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + 
     Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+
     Rugosity, CF_smry_TS_t) %>% vif()
# remove Pocillopora
lm(Feeding_intensity ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + 
     Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Porites.mas+
     Rugosity, CF_smry_TS_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Feeding_intensity ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Porites.mas+ Rugosity, 
                          data = CF_smry_TS_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Feeding_intensity_TS = lm(Feeding_intensity ~ 1, CF_smry_TS_t)
summary(lm_Feeding_intensity_TS)
shapiro.test(resid(lm_Feeding_intensity_TS))
par(mfrow=c(1,2))
qqnorm(resid(lm_Feeding_intensity_TS))
qqline(resid(lm_Feeding_intensity_TS), col = "red")
hist(resid(lm_Feeding_intensity_TS), breaks = 20, main = NULL)

#' 
#' #### TS-LM3: event duration
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Event_duration ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+ Rugosity, CF_smry_TS_t) %>% vif()
# Remove Pocillopora
lm(Event_duration ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Porites.mas+
     Rugosity, CF_smry_TS_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Event_duration ~ Biomass_Corallivores + Richness_Corallivores + Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Porites.mas+ Rugosity, 
                          data = CF_smry_TS_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Event_duration_TS = lm(Event_duration ~ Biomass_Piscivores + Acropora.tab, CF_smry_TS_t)
lm_Event_duration_TS %>% summary()
shapiro.test(resid(lm_Event_duration_TS))
par(mfrow=c(1,2))
qqnorm(resid(lm_Event_duration_TS))
qqline(resid(lm_Event_duration_TS), col = "red")
hist(resid(lm_Event_duration_TS), breaks = 20, main = NULL)

#' 
#' #### TS-LM4: richness of corallivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Richness_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+ Rugosity, CF_smry_TS_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Richness_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+ Rugosity, 
                          data = CF_smry_TS_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Richness_Corallivores_TS = lm(Richness_Corallivores ~ Acropora.tab + Pocillopora + Porites.mas, CF_smry_TS_t)
lm_Richness_Corallivores_TS %>% summary()
shapiro.test(resid(lm_Richness_Corallivores_TS))
par(mfrow=c(1,2))
qqnorm(resid(lm_Richness_Corallivores_TS))
qqline(resid(lm_Richness_Corallivores_TS), col = "red")
hist(resid(lm_Richness_Corallivores_TS), breaks = 20, main = NULL)

#' 
#' #### TS-LM5: biomass of corallivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Biomass_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+ Rugosity, CF_smry_TS_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Biomass_Corallivores ~ Biomass_Piscivores + Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+ Rugosity, 
                          data = CF_smry_TS_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Biomass_Corallivores_TS = lm(Biomass_Corallivores ~ Acropora.arb + Montipora + Pocillopora, CF_smry_TS_t)
lm_Biomass_Corallivores_TS %>% summary()
shapiro.test(resid(lm_Biomass_Corallivores_TS))
par(mfrow=c(1,2))
qqnorm(resid(lm_Biomass_Corallivores_TS))
qqline(resid(lm_Biomass_Corallivores_TS), col = "red")
hist(resid(lm_Biomass_Corallivores_TS), breaks = 20, main = NULL)

#' 
#' #### TS-LM6: biomass of piscivores
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Biomass_Piscivores ~ Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+ Rugosity, CF_smry_TS_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Biomass_Piscivores ~ Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas+ Rugosity, 
                          data = CF_smry_TS_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_Biomass_Piscivores_TS = lm(Biomass_Piscivores ~ Acropora.tab + Montipora, CF_smry_TS_t)
lm_Biomass_Piscivores_TS %>% summary()
shapiro.test(resid(lm_Biomass_Piscivores_TS))
par(mfrow=c(1,2))
qqnorm(resid(lm_Biomass_Piscivores_TS))
qqline(resid(lm_Biomass_Piscivores_TS), col = "red")
hist(resid(lm_Biomass_Piscivores_TS), breaks = 20, main = NULL)

#' 
#' #### TS-LM7: rugosity
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
lm(Rugosity ~ Acropora.tab + Acropora.cory + Acropora.arb + Montipora + Pocillopora+ Porites.mas, CF_smry_TS_t) %>% vif()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
model_selection = glmulti(Rugosity ~ Acropora.tab + Acropora.cory + Acropora.arb + 
                            Montipora + Pocillopora+ Porites.mas, 
                          data = CF_smry_TS_t,
                          level = 1,               # only main effects
                          method = "h",            # heuristic search
                          crit = "aicc",           # use AICc!
                          fitfunction = lm)        # for lm models

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
summary(model_selection)$bestmodel

lm_rugosity_TS = lm(Rugosity ~ Porites.mas, CF_smry_TS_t)
lm_rugosity_TS %>% summary()
shapiro.test(resid(lm_rugosity_TS))
par(mfrow=c(1,2))
qqnorm(resid(lm_rugosity_TS))
qqline(resid(lm_rugosity_TS), col = "red")
hist(resid(lm_rugosity_TS), breaks = 20, main = NULL)

#' 
#' ### SEM among GBR sites
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
SEM_GBR = psem(
  lm_total_corallivory_GBR,
  lm_Feeding_intensity_GBR,
  lm_Event_duration_GBR,
  lm_Richness_Corallivores_GBR,
  lm_Biomass_Corallivores_GBR,
  lm_Biomass_Piscivores_GBR,
  lm_rugosity_GBR
)
smry_SEM_GBR = summary(SEM_GBR, .progressBar = FALSE)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
smry_SEM_GBR # AIC = 796.53

#' 
#' ### SEM among TS sites
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
SEM_TS = psem(
  lm_total_corallivory_TS,
  lm_Event_duration_TS,
  lm_Richness_Corallivores_TS,
  lm_Biomass_Corallivores_TS,
  lm_Biomass_Piscivores_TS,
  lm_rugosity_TS
)
smry_SEM_TS = summary(SEM_TS, .progressBar = FALSE)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
smry_SEM_TS # AIC = 227.98

#' 
#' ## **Plotting**
#' 
#' ### Figure 2
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false

Fig2a = ggplot(CF_smry, aes(x = Reef, y = Total_corallivory)) +
  geom_boxplot(aes(fill = factor(interaction(Region, Reef)))) +
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#998ec3"),
                    labels = c("GBR-Heron", "GBR-Davies", "GBR-Lizard", "TS-Masig"))+
  scale_x_discrete(labels = c("Heron", "Davies", "Lizard", "Masig"))+
  scale_y_continuous(limits = c(0, 350), breaks = seq(0, 300, 100))+
  labs(x = NULL, y = "Total bites per 30 mins", fill = "Region-Reef")

cld = data.frame(
  Habitat = unique(CF_smry$Habitat),
  Letters = c("a", "b")
)
Fig2b = ggplot(CF_smry, aes(x = Habitat, y = Total_corallivory, fill = Habitat))+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank())+
  labs(x = NULL, y = "Total bites per 30 mins")+
  scale_x_discrete(labels = c("Back", "Fore"))+
  scale_y_continuous(limits = c(0, 350), breaks = seq(0, 300, 100))+
  geom_boxplot()+
  geom_text(data = cld, aes(x = Habitat, y = 325, label = Letters, group = Habitat), 
            position = position_dodge(width = 0.75), size = 5)

Fig2c = ggplot(CF_smry, aes(x = Reef, y = Feeding_intensity)) +
  geom_boxplot(aes(fill = factor(interaction(Region, Reef)))) +
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#998ec3"),
                    labels = c("GBR-Heron", "GBR-Davies", "GBR-Lizard", "TS-Masig"))+
  scale_x_discrete(labels = c("Heron", "Davies", "Lizard", "Masig"))+
  labs(x = NULL, y = expression(bold("Feeding intensity"~" (bites"~min^{-1}*")")),
       fill = "Region-Reef")

Fig2d = ggplot(CF_smry, aes(x = Habitat, y = Feeding_intensity, fill = Habitat))+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank())+
  labs(x = NULL, y = expression(bold("Feeding intensity"~" (bites"~min^{-1}*")")),
       fill = "Region-Reef")+
  scale_x_discrete(labels = c("Back", "Fore"))+
  geom_boxplot()

Fig2e = ggplot(CF_smry, aes(x = Reef, y = Event_duration)) +
  geom_boxplot(aes(fill = factor(interaction(Region, Reef)))) +
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#998ec3"),
                    labels = c("GBR-Heron", "GBR-Davies", "GBR-Lizard", "TS-Masig"))+
  scale_x_discrete(labels = c("Heron", "Davies", "Lizard", "Masig"))+
  labs(x = NULL, y = "Event duration (min)", fill = "Region-Reef")

Fig2f = ggplot(CF_smry, aes(x = Habitat, y = Event_duration, fill = Habitat))+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank())+
  labs(x = NULL, y = "Event duration (min)")+
  scale_x_discrete(labels = c("Back", "Fore"))+
  geom_boxplot()

Fig2g = ggplot(CF_smry, aes(x = Reef, y = Richness_Corallivores)) +
  geom_boxplot(aes(fill = factor(interaction(Region, Reef)))) +
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#998ec3"),
                    labels = c("GBR-Heron", "GBR-Davies", "GBR-Lizard", "TS-Masig"))+
  scale_x_discrete(labels = c("Heron", "Davies", "Lizard", "Masig"))+
  scale_y_continuous(breaks = seq(0,10,2))+
  labs(x = NULL, y = "Richness of corallivores")

Fig2h = 
  ggplot(CF_smry, aes(x = Habitat, y = Richness_Corallivores, fill = Habitat))+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank())+
  labs(x = NULL, y = "Richness of corallivores")+
  scale_x_discrete(labels = c("Back", "Fore"))+
  scale_y_continuous(breaks = seq(0,10,2))+
  geom_boxplot()

cld = data.frame(
  Reef = unique(CF_smry$Reef),
  Letters = c("a", "ab", "ab", "b")
)
Fig2i = ggplot(CF_smry[CF_smry$Plot %in% c(1,3),], aes(x = Reef, y = Biomass_Corallivores)) +
  geom_boxplot(aes(fill = factor(interaction(Region, Reef)))) +
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#998ec3"),
                    labels = c("GBR-Heron", "GBR-Davies", "GBR-Lizard", "TS-Masig"))+
  scale_x_discrete(labels = c("Heron", "Davies", "Lizard", "Masig"))+
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 6, 2))+
  labs(x = NULL, y = expression(bold("Biomass of corallivores"~" (g"~m^{-2}*")")))+
  geom_text(data = cld, aes(x = Reef, y = 6.5, label = Letters, group = Reef), 
            position = position_dodge(width = 0.75), size = 5)

cld = data.frame(
  Habitat = unique(CF_smry$Habitat),
  Letters = c("a", "b")
)
Fig2j = ggplot(CF_smry[CF_smry$Plot %in% c(1,3),], aes(x = Habitat, y = Biomass_Corallivores, fill = Habitat))+
  geom_boxplot()+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank())+
  labs(x = NULL, y = expression(bold("Biomass of corallivores"~" (g"~m^{-2}*")")))+
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 6, 2))+
  scale_x_discrete(labels = c("Back", "Fore"))+
  geom_text(data = cld, aes(x = Habitat, y = 6.5, label = Letters, group = Habitat), 
            position = position_dodge(width = 0.75), size = 5)


Fig2_1 = plot_grid(Fig2a, Fig2b,
                   ncol = 2, nrow = 1, label_x = -0.02, label_y = 1.05, labels = c("(a)", "(b)"))
Fig2_1 = ggdraw() + 
  draw_plot(Fig2_1, x = 0.05, y = 0.05, width = 0.9, height = 0.9)
Fig2_2 = plot_grid(Fig2c, Fig2d,
                   ncol = 2, nrow = 1, label_x = -0.02, label_y = 1.2, labels = c("(c)", "(d)"))
Fig2_2 = ggdraw() + 
  draw_plot(Fig2_2, x = 0.05, y = 0.05, width = 0.9, height = 0.9)
Fig2_3 = plot_grid(Fig2e, Fig2f,
                   ncol = 2, nrow = 1, label_x = -0.02, label_y = 1.2, labels = c("(e)", "(f)"))
Fig2_3 = ggdraw() + 
  draw_plot(Fig2_3, x = 0.05, y = 0.05, width = 0.9, height = 0.9)

Fig2_4 = plot_grid(Fig2g, Fig2h,
                   ncol = 2, nrow = 1, label_x = -0.02, label_y = 1.2, labels = c("(g)", "(h)"))
Fig2_4 = ggdraw() + 
  draw_plot(Fig2_4, x = 0.05, y = 0.05, width = 0.9, height = 0.9)

Fig2_5 = plot_grid(Fig2i, Fig2j,
                   ncol = 2, nrow = 1, label_x = -0.02, label_y = 1.2, labels = c("(i)", "(j)"))
Fig2_5 = ggdraw() + 
  draw_plot(Fig2_5, x = 0.05, y = 0.05, width = 0.9, height = 0.9)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 13
#| fig-align: "center"
#| warning: false
plot_grid(Fig2_1, Fig2_2, Fig2_3, Fig2_4, Fig2_5, ncol = 1, nrow = 5)
ggsave("Plot/Fig2.png", width = 8, height = 13, units = "in", dpi = 600, bg = "white")

#' 
#' ### Figure 3
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
options(ggrepel.max.overlaps = Inf)
Fig3 = ggplot(df_nmds_bite_sp, aes(x = NMDS1, y = NMDS2))+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank())+
  geom_point(size = 4, position = position_jitter(h=0.15,w=0.15), aes(colour = factor(interaction(Region, Reef)), shape = Habitat))+
  labs(x = "NMDS1", y = "NMDS2", colour = "Region-Reef")+
  annotate(geom="text", x= 1.4, y=-1.2, label= paste("Stress =", round(nmds_bite_sp$stress,2)),
           size = 4, fontface = "bold")+
  scale_color_manual(values = c("#238443", "#78c679", "#c2e699", "#998ec3"),
                     labels = c("GBR-Heron", "GBR-Davies", "GBR-Lizard", "TS-Masig"))+
  scale_shape_discrete(labels = c('Back', 'Fore'))+
  # Add vectors for significant species
  geom_segment(data = sps_bite_sp, 
               aes(x = 0, y = 0, xend = NMDS1*0.9, yend = NMDS2*0.9), 
               arrow = arrow(length = unit(0.1, "cm")), 
               colour = "grey", linewidth = 0.7) +
    # Add labels at the end of vectors
  geom_text_repel(data = sps_bite_sp, aes(x = NMDS1, y = NMDS2, label = Valid_as), 
                  size = 4, colour = "#606060", fontface = "bold.italic", max.overlaps = 15)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 6
#| fig-align: "center"
#| warning: false
Fig3
ggsave("Plot/Fig3.png", width = 10, height = 6, units = "in", dpi = 600, bg = "white")

#' 
#' ### Figure 4
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
Fig4 = ggplot(df_nmds_ben, aes(x = NMDS1, y = NMDS2))+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank())+
  geom_point(size = 4, aes(colour = factor(interaction(Region, Reef)), shape = Habitat))+
  labs(x = "NMDS1", y = "NMDS2", color = "Region-Reef")+
  annotate(geom="text", x=0.4, y=-0.5, label= paste("Stress =", round(nmds_ben$stress,2)),
           size = 4, fontface = "bold")+
  scale_color_manual(values = c("#238443", "#78c679", "#c2e699", "#998ec3"),
                     labels = c("GBR-Heron", "GBR-Davies", "GBR-Lizard", "TS-Masig"))+
  scale_shape_discrete(labels = c('Back', 'Fore'))+
  # Add vectors for significant species
  geom_segment(data = sps_ben, 
               aes(x = 0, y = 0, xend = NMDS1*0.9, yend = NMDS2*0.9), 
               arrow = arrow(length = unit(0.1, "cm")), 
               colour = "grey", linewidth = 0.7) +
  
  # Add labels at the end of vectors
  geom_text_repel(data = sps_ben, 
                  aes(x = NMDS1, y = NMDS2, label = Benthos), 
                  colour = "#606060", size = 4, fontface = "bold")

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 6
#| fig-align: "center"
#| warning: false
Fig4
ggsave("Plot/Fig4.png", width = 10, height = 6, units = "in", dpi = 600, bg = "white")

#' 
#' ### Figure 5
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# Extract coefficient from SEM for plotting
var_list = data.frame(
  var = c("Total_corallivory", "Feeding_intensity", "Event_duration", "Richness_Corallivores", "Biomass_Corallivores", 
          "Biomass_Piscivores", "Rugosity", 
          "Acropora.tab", "Acropora.cory", "Acropora.arb", "Montipora", "Pocillopora",
          "Porites.mas"),
  lab = c("Total corallivory", "Feeding intensity", "Event duration", 
          "Corallivore richness", "Corallivore biomass",
          "Piscivore biomass", "Rugosity", "Tabular Acropora",
          "Corymbose Acropora", "Arborescent Acropora", "Montipora", "Pocillopora",
          "Massive Porites")
)
# Extract outputs of SEM_GBR
net_GBR  = coefs(SEM_GBR, standardize = "none") %>% as.data.frame()
net_GBR = net_GBR[,-8] %>% filter(P.Value<0.05) %>% dplyr::select(Predictor, Response, Estimate) %>% as.data.frame()

for (j in 1:nrow(net_GBR)) {
  net_GBR$Predictor[j] = var_list$lab[var_list$var == net_GBR$Predictor[j]]
  net_GBR$Response[j] = var_list$lab[var_list$var == net_GBR$Response[j]]
}
# Extract outputs of SEM_TS
net_TS  = coefs(SEM_TS, standardize = "none") %>% as.data.frame()
net_TS = net_TS[,-8] %>% filter(P.Value<0.05) %>% dplyr::select(Predictor, Response, Estimate) %>% as.data.frame()

for (j in 1:nrow(net_TS)) {
  net_TS$Predictor[j] = var_list$lab[var_list$var == net_TS$Predictor[j]]
  net_TS$Response[j] = var_list$lab[var_list$var == net_TS$Response[j]]
}

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 7
#| fig-height: 8
#| fig-align: "center"
#| warning: false
par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))  # Set layout (2 rows, 1 column)

custom_layout <- matrix(
  c(
    -7,1,
    7,1,
    10,0,
    -9,-2,
    -1.5,-2,
    -5,-1,
    -17,-2,
    12,-2,
    17,-2,
    6,-2,
    0,2,
    -10,0,
    5,-1
  ), 
  ncol = 2, byrow = TRUE
)
qgraph(net_GBR[,1:3], 
       layout = custom_layout,
       color = c("white", "white", "white","white","white"),
       groups = c(rep("Behavior", 2), "Corallivore metrics", rep("Coral", 2), "Predator & Structure", 
                  rep("Coral", 4), "Total bites", "Corallivore metrics", "Predator & Structure"),
       legend = F,
       shape = "rectangle",
       asize = 4,
       weighted = T,
       #curve = 0.4,
       #curveAll = T,
       vsize = c(nchar(unique(c(net_GBR$Predictor, net_GBR$Response)))+2),
       vsize2 = 3,
       esize = 8, 
       edge.color = ifelse(net_GBR$Estimate > 0, "#92c5de", ifelse(net_GBR$Estimate < 0, "#f4a582", NA)),
       fade=FALSE,
       edge.labels = TRUE,
       edge.label.cex = 0.8,
       edge.label.bg = T,
       edge.label.color = "black",
       edge.label.position = 0.77,
       label.cex = 1.2,
       label.scale.equal = T, 
       label.font = 2)
mtext("(a)", side = 3, line = 1, adj = 0, font = 2, cex = 1)  # Add label (a)

custom_layout <- matrix(
  c(
    0.6,-2,
    0.25,-1,
    -0.4, -2,
    1,-2,
    0.15,-2,
    -1,-2,
    0,2,
    0.25,1,
    0.5,0,
    -0.5,0
  ), 
  ncol = 2, byrow = TRUE
)
qgraph(net_TS[,1:3], 
       color = c("white", "white", "white","white","white"),
       groups = c("Coral", "Predator & Structure", rep("Coral",4), "Total corallivory",
                  "Behavior", rep("Corallivore metrics", 2)),
       layout = custom_layout,
       legend = F,
       shape = "rectangle",
       asize = 4,
       weighted = T,
       #curve = 0.7,
       #curveAll = T,
       vsize = c(nchar(unique(c(net_TS$Predictor, net_TS$Response)))+2),
       vsize2 = 3,
       esize = 8, 
       edge.color = ifelse(net_TS$Estimate > 0, "#92c5de", ifelse(net_TS$Estimate < 0, "#f4a582", NA)),
       fade=FALSE,
       edge.labels = TRUE,
       edge.label.cex = 0.8,
       edge.label.bg = T,
       edge.label.color = "black",
       edge.label.position = 0.77,
       label.cex = 1.2, 
       label.scale.equal = T,
       label.font = 2,
       repulsion = 1.2)
mtext("(b)", side = 3, line = 2, adj = 0, font = 2, cex = 1)  # Add label (b)

#' 
#' ### Marginal effect plots of LMs
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false
# Create summary table of mean and sd for each variable
var_s = rbind(
  data.frame(
    Mean = colMeans(CF_smry_GBR[, var_list$var], na.rm = TRUE),
    SD   = apply(CF_smry_GBR[, var_list$var], 2, sd, na.rm = TRUE)
  ),
  data.frame(
    Mean = colMeans(CF_smry_TS[, var_list$var], na.rm = TRUE),
    SD   = apply(CF_smry_TS[, var_list$var], 2, sd, na.rm = TRUE)
  )
)
rownames(var_s) = NULL
var_s$var = var_list$var
var_s$region = rep(c("GBR", "TS"), each = 13)

# These following variables were transformed
var_s[1,1] = mean((CF_smry_GBR$Total_corallivory)^(1/3), na.rm = TRUE)
var_s[1,2] = sd((CF_smry_GBR$Total_corallivory)^(1/3), na.rm = TRUE)
var_s[6,1] = mean((CF_smry_GBR$Biomass_Piscivores)^(1/4), na.rm = TRUE)
var_s[6,2] = sd((CF_smry_GBR$Biomass_Piscivores)^(1/4), na.rm = TRUE)
var_s[14,1] = mean((CF_smry_TS$Total_corallivory)^(1/3), na.rm = TRUE)
var_s[14,2] = sd((CF_smry_TS$Total_corallivory)^(1/3), na.rm = TRUE)
var_s[16, 1] = mean((CF_smry_TS$Event_duration)^(1/4), na.rm = TRUE)
var_s[16,2] = sd((CF_smry_TS$Event_duration)^(1/4), na.rm = TRUE)

# Aggregate all LMs across GBR and TS
allLM = list(lm_total_corallivory_GBR,
              lm_Feeding_intensity_GBR,
              lm_Event_duration_GBR,
              lm_Richness_Corallivores_GBR,
              lm_Biomass_Corallivores_GBR,
              lm_Biomass_Piscivores_GBR,
              lm_rugosity_GBR,
              lm_total_corallivory_TS,
              #lm_Feeding_intensity_TS,
              lm_Event_duration_TS,
              lm_Richness_Corallivores_TS,
              lm_Biomass_Corallivores_TS,
              lm_Biomass_Piscivores_TS,
              lm_rugosity_TS)

ResV = rep(var_list$var[1:7], 2)
ResV = ResV[-9] # no predictor for lm_Feeding_intensity_TS

# Initialize the master list
plot_ME = vector("list", length = length(allLM))

c("Total corallivory", "Feeding intensity", "Event duration", 
          "Corallivore richness", "Corallivore biomass",
          "Piscivore biomass", "Rugosity", "Tabular Acropora",
          "Corymbose Acropora", "Arborescent Acropora", "Montipora", "Pocillopora",
          "Massive Porites")

axis_title = list(
  Total_corallivory = expression(Total~corallivory~(bites~30~min^{-1})),
  Feeding_intensity = expression(Feeding~intensity~(bites~min^{-1})),
  Event_duration = expression(Event~duration~(min)),
  Richness_Corallivores = expression(Corallivore~richness),
  Biomass_Corallivores = expression(Corallivore~biomass~(g~m^{-2})),
  Biomass_Piscivores = expression(Piscivore~biomass~(kg~m^{-2})),
  Rugosity = expression(Rugosity),
  Acropora.tab = expression(tabular~italic(Acropora)~("%")),
  Acropora.cory = expression(corymbose~italic(Acropora)~("%")),
  Acropora.arb = expression(arborescent~italic(Acropora)~("%")),
  Montipora = expression(italic(Montipora)~("%")),
  Pocillopora = expression(italic(Pocillopora)~("%")),
  Porites.mas = expression(massive~italic(Porites)~("%"))
)

# Helper function to handle back-transformations by variable name
back_transform = function(x, var_name, stats_df, Region) {
  mu = stats_df[stats_df$var == var_name & stats_df$region == Region, "Mean"]
  sigma = stats_df[stats_df$var == var_name & stats_df$region == Region, "SD"]
  val = (x * sigma) + mu # Reverse the standardization
  
  # Reverse the power transformation based on the specific variable name
  if (var_name == "Total_corallivory") {
    return(val^3)
  } else if (var_name == "Biomass_Piscivores" && Region == "GBR") {
    return(val^4)
  }
  else if (var_name == "Event_duration" && Region == "TS") {
    return(val^4)
  }
  return(val)
}

# Main Loop
for (n in seq_along(allLM)) {
  mod = allLM[[n]]
  resp_var = ResV[n]
  mod_coefs = summary(mod)$coefficients
  pred_vars = rownames(mod_coefs)[-1] # Grab all predictors, excluding the Intercept
  
  # Partition data for either region
  if (n %in% 1:7) {
    df_raw = CF_smry_GBR
    Zone = "GBR"
  } else {
    df_raw = CF_smry_TS
    Zone = "TS"
  }
  
  # Create empty lists for storing marginal plots
  plot_ME[[n]] = vector("list", length = length(pred_vars))
  
    for (j in seq_along(pred_vars)) {
    pred_var = pred_vars[j]
    eff = mod_coefs[j + 1, 1]
    p_val = mod_coefs[j + 1, "Pr(>|t|)"]
    
    # Get marginal effects
    me = ggpredict(mod, terms = paste0(pred_var, "[n=20]"))
    
    # Apply back-transformations using the helper function
    me$x         = back_transform(me$x, pred_var, var_s, Zone)
    me$predicted = back_transform(me$predicted, resp_var, var_s, Zone)
    me$conf.low  = back_transform(me$conf.low, resp_var, var_s, Zone)
    me$conf.high = back_transform(me$conf.high, resp_var, var_s, Zone)
    
    # Map variable names for plotting
    raw_x = var_list$var[var_list$var == pred_var]
    raw_y = var_list$var[var_list$var == resp_var]
    lab_x = axis_title[[pred_var]]
    lab_y = axis_title[[resp_var]]
    

    # Build the BASE plot (Points, Theme, Labels)
    p = ggplot() +
      geom_point(data = df_raw, aes(x = .data[[raw_x]], y = .data[[raw_y]])) +
      theme_minimal() +
      theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
      labs(x = lab_x, y = lab_y)
    
    # Conditionally layer the ribbon and line ONLY if significant
    if (p_val < 0.05) {
      ribbon_color = ifelse(eff > 0, "#92c5de", "#f4a582")
      p = p + 
        geom_ribbon(data = me, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = ribbon_color) +
        geom_line(data = me, aes(x = x, y = predicted))
    }
    
    # Save the final plot to the list
    plot_ME[[n]][[j]] = p
  }
}

#' 
#' #### Figure S1
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 12
#| fig-height: 8
#| fig-align: "center"
#| warning: false
plot_grid(plotlist = unlist(plot_ME[1:7], recursive = FALSE), ncol = 6)
ggsave("Plot/FigS1.png", width = 12, height = 8, units = "in", dpi = 600, bg = "white")

#' 
#' #### Figure S2
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 6
#| fig-align: "center"
#| warning: false
plot_grid(plotlist = unlist(plot_ME[8:14], recursive = FALSE), ncol = 4)
ggsave("Plot/FigS2.png", width = 8, height = 6, units = "in", dpi = 600, bg = "white")

#' 
#' ### Biomass of corallivores \~ coral cover
#' 
#' #### Figure S3
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| output: false

# Total hard coral cover
CF_smry$HardCoralCover = Ben_MFG %>% ungroup() %>% select(
  Acro.arb, Acro.bran, Acro.cory, Acro.tab, Fire.coral, Foliose.lett, Foliose.plat,
  Freeliving, HC.other.bran, HC.other.mas, Montipora, Pocillopora, Porites.mas, Stylophora
) %>% rowSums()

# Total cover pf the major six coral prey
CF_smry = CF_smry %>% mutate(
  M6prey = rowSums(pick(`Acropora.arb`:`Pocillopora`), na.rm = TRUE)
)

CF_smry_site = CF_smry %>% group_by(Region, Reef, Site, TR) %>% summarise(across(where(is.numeric), mean))

require(ggpmisc)

FigS3a = ggplot(CF_smry_site, aes(x = HardCoralCover, y = Biomass_Corallivores))+
  theme_minimal()+
  geom_point()+
  geom_smooth(method = "lm", fill = "#92c5de", col = "black", linewidth = 0.65)+
  stat_poly_eq(formula = y ~ x,
               use_label(c("R2", "p")), small.p = TRUE)+
  labs(x = "Hard coral cover (%)", y = expression("Biomass of corallivores (g "~m^{-2}*")"))

FigS3b = 
ggplot(CF_smry_site, aes(x = M6prey, y = Biomass_Corallivores))+
  theme_minimal()+
  geom_point()+
  geom_smooth(method = "lm", fill = "#92c5de", col = "black", linewidth = 0.65)+
  stat_poly_eq(formula = y ~ x,
               use_label(c("R2", "p")), small.p = TRUE)+
  labs(x = "Cover of six major coral prey (%)", y = expression("Biomass of corallivores (g "~m^{-2}*")"))

FigS3 = plot_grid(FigS3a, FigS3b, labels = c("(a)", "(b)"))

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 4
#| fig-align: "center"
#| warning: false
FigS3
ggsave("Plot/FigS3.png", units = "in", width = 8, height = 4, dpi = 600, bg = "white")

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| eval: true
#| output: false

# Automatically extract code after rendering
knitr::purl(input = "CorallivoryDrivers_Hsu_et_al_2026.qmd", output = "CorallivoryDrivers_Hsu_et_al_2026.R", documentation = 2, quiet = TRUE)

