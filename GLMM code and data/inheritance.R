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

resultz <- read.csv2(file = "inheritance.csv")

resultz$driveParent <- factor(resultz$driveParent, 
                              levels = c("male", "female"), 
                              labels = c(intToUtf8(9794),intToUtf8(9792)))
resultz$line <- factor(resultz$line, 
                              levels = c("1A", "1D"), 
                              labels = c("1A", "1D"))
resultz$cross <- factor(resultz$cross, 
                       levels = c("control", "drive"), 
                       labels = c("Control", "Drive"))
resultz$batch <- rownames(resultz)
resultz$batch <- factor(resultz$batch)
resultz$group <- paste(resultz$driveParent, resultz$line, resultz$cross, sep = "_")
resultz$group <- factor(resultz$group)

resultz <- resultz %>%
  mutate(drInh = drive / offspring)

##################################
########## GLM ###################
##################################

m_full1 <- glmmTMB(drInh ~ 1 + cross * driveParent * line + (1 | batch),
                   family = binomial, 
                   data = resultz,
                   weights = offspring)
summary(m_full1)
plot(simulateResiduals(m_full1))
Anova(m_full1)

m_full2 <- glmmTMB(drInh ~ 1 + cross * driveParent * line + (cross | batch),
                  family = binomial, 
                  data = resultz,
                  weights = offspring)
summary(m_full2)
plot(simulateResiduals(m_full2))
Anova(m_full2)

anova(m_full1, m_full2)

glm.drive1 <- glmmTMB(drInh ~ 1 + cross + driveParent + (1 | batch),
                     family = binomial, 
                     data = resultz,
                     weights = offspring)
summary(glm.drive1)
plot(simulateResiduals(glm.drive1))
Anova(glm.drive)

glm.drive2 <- glmmTMB(drInh ~ 1 + cross + driveParent + (cross | batch),
                 family = binomial, 
                 data = resultz,
                 weights = offspring)
summary(glm.drive2)
plot(simulateResiduals(glm.drive2))
Anova(glm.drive2)

anova(m_full1, m_full2, glm.drive1, glm.drive2)
# glm.drive2 is the final model

##########################################
########## Group means and CIs ###########
##########################################

m_groupmeans <- glmmTMB(drInh ~ group + (cross | batch) - 1,
                        family = binomial, 
                        data = resultz,
                        weights = offspring)
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

#########################################
########## Plot plots ###################
#########################################

p1 <- ggplot(data = resultz, aes(x = interaction(line, driveParent), group = cross)) +
  geom_point(aes(y = drInh * 100, size = offspring, colour = line, shape = cross), alpha = 0.75,
             position=position_jitterdodge(jitter.width = 0.5, dodge.width = 1, seed = 8)) +
  geom_text(data = resultzSummary, aes(y = c(40,95,35,100,40,95,35,90), label = sprintf("%.02f%% ± %.02f%%", meanInheritance * 100, sem * 100)), 
            position=position_dodge(width = 0.5), size = 3) +
  geom_errorbar(data = resultzSummary, aes(y = meanInheritance * 100, ymin = (meanInheritance - sem) * 100, ymax = (meanInheritance + sem) * 100), 
                linewidth = 0.5, width = 0.25, position=position_dodge(width = 1)) +
  scale_size(name = "Number of offspring") +
  scale_shape(name = "Cross") +
  scale_color_manual(values = c(met.brewer("Hiroshige",2)), name = "Drive line") +
  ylim(c(0,100)) +
  ylab("Drive inheritance (%)") +
  xlab("Sex") +
  PaperTheme +
  guides(colour = guide_legend(order = 2), size = guide_legend(order = 1), shape = guide_legend(order = 3)) +
  theme(legend.box = "horizontal",
        strip.placement = "outside",
        strip.text.x = element_text(size = 10),
        legend.position.inside = c(0.5, 0.15))
p1

ggsave(plot = p1, filename = "Fig8.pdf", height = 13, width = 12, unit = "cm", dpi = 600)



