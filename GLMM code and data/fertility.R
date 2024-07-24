# Created by Nicky Faber, 2024.

##############################
## ---- Working directory ----
##############################

rm(list = ls())
setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/02 Evolutionary stability/HaploLethalFertilityDrive/GLMM code and data")

library(tidyverse)
library(multcompView)
library(MetBrewer)
library(showtext); showtext_auto(); showtext_opts(dpi=1000)
library(glmmTMB)
library(DHARMa)
library(car)

# Figure aesthetics
PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title=element_text(size=12, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12))

#########################################
########## Wrangle data #################
#########################################

resultz <- read.csv2(file = "fertility.csv")

resultz$cross <- factor(resultz$cross, 
                             levels = c("R", "RG"), 
                             labels = c("noDriveParent", "driveParent"))
resultz$driveParent <- factor(resultz$driveParent, 
                              levels = c("male", "female"), 
                              labels = c(intToUtf8(9794),intToUtf8(9792)))
resultz$line <- factor(resultz$line, 
                              levels = c("1A", "1D"), 
                              labels = c("1A", "1D"))
resultz$batch <- rownames(resultz)
resultz$batch <- factor(resultz$batch)
resultz$group <- paste(resultz$driveParent, resultz$line, resultz$cross, sep = "_")
resultz$group <- factor(resultz$group)

resultz <- resultz %>%
  mutate(eggToAdult = adultOffspring / eggOffspring, 
         drInh = drive / adultOffspring)

##################################
########## GLM ###################
##################################

m_full1 <- glmmTMB(eggToAdult ~ 1 + cross * driveParent * line + (1 | batch),
                  family = binomial, 
                  data = resultz,
                  weights = eggOffspring)
summary(m_full1)
plot(simulateResiduals(m_full1))
Anova(m_full1)

# Removed some insignificant interactions due to model convergence issues (overfitting)
m_full2 <- glmmTMB(eggToAdult ~ 1 + cross * driveParent + cross * line + (driveParent | batch),
                   family = binomial, 
                   data = resultz,
                   weights = eggOffspring)
summary(m_full2)
plot(simulateResiduals(m_full2))
Anova(m_full2)

# Removed some insignificant interactions due to model convergence issues (overfitting)
m_full3 <- glmmTMB(eggToAdult ~ 1 + cross * driveParent + line + (cross | batch),
                   family = binomial, 
                   data = resultz,
                   weights = eggOffspring)
summary(m_full3)
plot(simulateResiduals(m_full3))
Anova(m_full3)

glm.eggSurvival1 <- glmmTMB(eggToAdult ~ 1 + cross * driveParent + (1 | batch),
                   family = binomial, 
                   data = resultz,
                   weights = eggOffspring)
summary(glm.eggSurvival1)
plot(simulateResiduals(glm.eggSurvival1))
Anova(glm.eggSurvival1)

glm.eggSurvival2 <- glmmTMB(eggToAdult ~ 1 + cross * driveParent + cross * line + (driveParent | batch),
                  family = binomial, 
                  data = resultz,
                  weights = eggOffspring)
summary(glm.eggSurvival2)
plot(simulateResiduals(glm.eggSurvival2))
Anova(glm.eggSurvival2)

glm.eggSurvival3 <- glmmTMB(eggToAdult ~ 1 + cross * driveParent + (cross | batch),
                   family = binomial, 
                   data = resultz,
                   weights = eggOffspring)
summary(glm.eggSurvival3)
plot(simulateResiduals(glm.eggSurvival3))
Anova(glm.eggSurvival3)

anova(m_full1, m_full2, m_full3, glm.eggSurvival1, glm.eggSurvival2, glm.eggSurvival3)
# glm.eggSurvival2 is the final model (and it the same as m_full2)

##########################################
########## Group means and CIs ###########
##########################################

m_groupmeans <- glmmTMB(eggToAdult ~ group + (driveParent | batch) - 1,
                        family = binomial, 
                        data = resultz,
                        weights = eggOffspring)
summary(m_groupmeans)
Anova(m_groupmeans)

probabilities <- rownames_to_column(as_tibble(summary(m_groupmeans)$coefficients$cond, rownames = NA), "group") %>%
  mutate(group = str_sub(group, start = 6L)) %>%
  mutate(meanEggToAdult = exp(Estimate) / (1 + exp(Estimate)),
         sem = exp(Estimate + `Std. Error`) / (1 + exp(Estimate + `Std. Error`)) - meanEggToAdult)

resultzSummary <- resultz %>%
  group_by(driveParent, line, cross, group) %>%
  summarise() %>%
  right_join(probabilities, by='group')

#########################################
########## Plot plots ###################
#########################################

p1 <- ggplot(data = resultz, aes(x = driveParent, group = line)) +
  facet_grid(. ~ cross, labeller = label_both) +
  geom_point(aes(y = eggToAdult * 100, size = eggOffspring, colour = line), alpha = 0.75,
             position=position_jitterdodge(jitter.width = 0.5, dodge.width = 1, seed = 42)) +
  geom_text(data = resultzSummary, aes(y = c(60,60,50,50,60,60,50,50), label = sprintf("%.02f%% ± %.02f%%", meanEggToAdult * 100, sem * 100)), 
            position=position_dodge(width = 0.8), size = 3) +
  geom_errorbar(data = resultzSummary, aes(y = meanEggToAdult * 100, ymin = (meanEggToAdult - sem) * 100, ymax = (meanEggToAdult + sem) * 100), 
                linewidth = 0.5, width = 0.25, position=position_dodge(width = 1)) +
  ylim(c(-5,100)) +
  scale_size(name = "Number of eggs") +
  scale_color_manual(values = met.brewer("Hiroshige",2), name = "Drive line") +
  ylab("Egg to adult viability (%)") +
  xlab("Sex") +
  PaperTheme +
  guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
  theme(legend.box = "horizontal",
        strip.placement = "outside",
        axis.text.x = element_text(size = 20),
        strip.text.x = element_text(size = 12))
p1

ggsave(plot = p1, filename = "Fig10.pdf", height = 12, width = 20, unit = "cm", dpi = 1000)

##########################################
########## Group means and CIs ###########
##########################################

m_groupmeans <- glmmTMB(drInh ~ group + (cross | batch) - 1,
                        family = binomial, 
                        data = resultz,
                        weights = adultOffspring)
summary(m_groupmeans)
Anova(m_groupmeans)

probabilities <- rownames_to_column(as_tibble(summary(m_groupmeans)$coefficients$cond, rownames = NA), "group") %>%
  mutate(group = str_sub(group, start = 6L)) %>%
  mutate(meanInheritance = exp(Estimate) / (1 + exp(Estimate)),
         sem = exp(Estimate + `Std. Error`) / (1 + exp(Estimate + `Std. Error`)) - meanInheritance)

resultzSummary <- resultz %>%
  group_by(driveParent, line, cross, group) %>%
  summarise() %>%
  right_join(probabilities, by='group')
resultzSummary$sem <- replace_na(resultzSummary$sem, 0)

#########################################
########## Plot plots ###################
#########################################

p2 <- ggplot(data = resultz, aes(x = driveParent, group = line)) +
  facet_grid(. ~ cross, labeller = label_both) +
  geom_point(aes(y = drInh * 100, size = adultOffspring, colour = line), alpha = 0.75,
             position=position_jitterdodge(jitter.width = 0.5, dodge.width = 1, seed = 2)) +
  geom_text(data = resultzSummary, aes(y = c(20,40,30,50,20,40,30,50), label = sprintf("%.02f%% ± %.02f%%", meanInheritance * 100, sem * 100)), 
            position=position_dodge(width = 0.8), size = 3) +
  geom_errorbar(data = resultzSummary, aes(y = meanInheritance, ymin = (meanInheritance - sem) * 100, ymax = (meanInheritance + sem) * 100), 
                linewidth = 0.5, width = 0.25, position=position_dodge(width = 1)) +
  ylim(c(0,100)) +
  scale_size(name = "Number of offspring") +
  scale_color_manual(values = met.brewer("Hiroshige",2), name = "Drive line") +
  ylab("Drive inheritance (%)") +
  xlab("Sex") +
  PaperTheme +
  guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
  theme(legend.box = "horizontal",
        strip.placement = "outside",
        axis.text.x = element_text(size = 20),
        strip.text.x = element_text(size = 12))
p2

ggsave(plot = p2, filename = "FigS12.pdf", height = 12, width = 20, unit = "cm", dpi = 1000)






