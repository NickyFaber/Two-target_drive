# Created by Nicky Faber, 2022.
# Adapted from SLiM-Extras at https://github.com/MesserLab/SLiM-Extras.

##############################
## ---- Working directory ----
##############################

rm(list = ls())
setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/02 Evolutionary stability/HaploLethalFertilityDrive/fig3_r1_alleles")

library('foreach', quietly = TRUE)
library('doParallel', quietly = TRUE)
library('future', quietly = TRUE)
library('tidyverse', quietly = TRUE)
library('dplyr', quietly = TRUE)
library('ggplot2', quietly = TRUE)
library('MetBrewer', quietly = TRUE)
library('vcfR', quietly = TRUE)
library('patchwork', quietly = TRUE)
source("../Parse_SLiM_output.R")

# Figure aesthetics
PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title=element_text(size=14, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12))
options(scipen = 99999)

#############################
## ---- Parameter inputs ----
#############################

input <- list()

input$REPETITION <- 1:50
input$CAPACITY <- 100000
input$DRIVE_TYPE <- c("fertility",
                      "Xlinkedfertility",
                      "recessivelethal",
                      "fertility_fertility",
                      "haplolethal_fertility",
                      "haplolethal_recessivelethal",
                      "haplolethal_Xlinkedfertility",
                      "recessivelethal_fertility",
                      "recessivelethal_recessivelethal",
                      "recessivelethal_Xlinkedfertility")
input$DRIVE_CONVERSION_RATE <- 0.9
input$TOTAL_CUT_RATE <- 1
input$EMBRYO_RESISTANCE_CUT_RATE <- 0.1
input$SOMATIC_EXPRESSION_FITNESS_MULTIPLIER <- 0.9
input$R1_OCCURRENCE_RATE <- c(0,10^c(-6,-5,-4,-3,-2,-1))

inputs <- expand.grid(input)

#####################################
## ---- Run parallel simulations ---- 
#####################################

cl <- makeCluster(8, outfile="")
registerDoParallel(cl)

# Run SLiM in parallel:
raw_slim_output_matrix <- foreach(i=iter(1:nrow(inputs), by='row')) %dopar% {
  cat(sprintf("Working on iteration %i out of %i\n", i, nrow(inputs)))
  system(paste("slim",
               sprintf("-d CAPACITY=%s", inputs[i,2]),
               sprintf("-d \"DRIVE_TYPE='%s'\"", inputs[i,3]),
               sprintf('-d DRIVE_CONVERSION_RATE=%f', inputs[i,4]),
               sprintf('-d TOTAL_CUT_RATE=%f', inputs[i,5]),
               sprintf('-d EMBRYO_RESISTANCE_CUT_RATE=%f', inputs[i,6]),
               sprintf('-d SOMATIC_EXPRESSION_FITNESS_MULTIPLIER=%f', inputs[i,7]),
               sprintf('-d R1_OCCURRENCE_RATE=%f', inputs[i,8]),
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig3_r1_alleles.Rdata")

##############################################
########## Parse and calculate stuff #########
##############################################

model_output <- parse_slim_output(raw_slim_output_matrix, generations = 111, inputs)

model_output <- model_output %>%
  group_by(DRIVE_TYPE, repetition) %>%
  mutate(genetic_load = 1 - lead(popSize) * bonus_pop_factor / expectedPopSize)

#########################################
########## Save model ###################
#########################################

save(model_output, file = "Fig3_r1_alleles.Rdata")
# load(file = "Fig3_r1_alleles.Rdata")

#########################################
########## Plot plots ###################
#########################################

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility"))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue with distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue with distant-site')

model_output_gen110 <- main_drives[main_drives$generation==110,] %>%
  mutate(suppressed = case_when(popSize == 0 ~ 1,
                                popSize > 0 ~ 0)) %>%
  group_by(DRIVE_TYPE, R1_OCCURRENCE_RATE) %>%
  summarise(suppressionRate = sum(suppressed)/50)
model_output_gen110$R1_OCCURRENCE_RATE[model_output_gen110$R1_OCCURRENCE_RATE == 0] <- 10^-7

p1 <- ggplot(data = model_output_gen110) +
  geom_line(aes(x = log(R1_OCCURRENCE_RATE, 10), y = suppressionRate, colour = DRIVE_TYPE), na.rm = TRUE) +
  scale_color_manual(values = met.brewer("Hiroshige", 3), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  scale_x_continuous(breaks = seq(-7, -1, 1), 
                     labels = c("-Inf", seq(-6, -1, 1))) +
  xlab("Log(r1 formation rate)") +
  ylab("Complete suppression\nsuccess rate") +
  PaperTheme + theme(legend.position = "bottom")
p1

ggsave(plot = p1, filename = "Fig6.png", height = 10, width = 9, unit = "cm", dpi = 600)

p2 <- ggplot(data = main_drives) +
  facet_grid(. ~ R1_OCCURRENCE_RATE) +
  geom_line(aes(x = generation, y = popSize, group = interaction(DRIVE_TYPE, repetition), colour = DRIVE_TYPE), alpha = 0.4, linewidth = 0.4) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 3), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ylab("Population size") +
  ggtitle("Horizontal = Relative functional resistance rate") +
  PaperTheme + theme(legend.position = "bottom")
p2

ggsave(plot = p2, filename = "FigS6.png", height = 8, width = 15, unit = "cm", dpi = 600)

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility"))
main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression\ndrive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')
model_output_subset_drive <- pivot_longer(main_drives, cols = c(freq_wt, freq_dr, freq_r2, freq_r1),
                                          names_to = "allele", values_to = "freq")
model_output_subset_drive$allele <- factor(model_output_subset_drive$allele,
                                           levels = c("freq_wt", "freq_dr", "freq_r2", "freq_r1"),
                                           labels = c("Wild type","Drive","Non-functional (r2)","Functional resistance (r1)"))

p3 <- ggplot(data = model_output_subset_drive) +
  facet_grid(DRIVE_TYPE ~ R1_OCCURRENCE_RATE) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = c(rev(met.brewer("Hiroshige", 3)),"red"), name = "") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Horizontal = Relative functional resistance rate\nHoming site") +
  PaperTheme
p3

model_output_subset_offtarget <- pivot_longer(main_drives, cols = c("freq_offtarget_wt", "freq_offtarget_disrupted"),
                                              names_to = "allele", values_to = "freq")
model_output_subset_offtarget$allele <- factor(model_output_subset_offtarget$allele,
                                               labels = c("Wild type","Non-functional (r2)"))

p4 <- ggplot(data = model_output_subset_offtarget) +
  facet_grid(DRIVE_TYPE ~ R1_OCCURRENCE_RATE) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 2), name = "") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Distant site") +
  PaperTheme + theme(legend.position = "none")
p4

p <- (p3 / p4) + plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "FigS7.png", height = 30, width = 20, unit = "cm", dpi = 600)

model_output_labels <- model_output
model_output_labels$DRIVE_TYPE <- recode(model_output_labels$DRIVE_TYPE, 
                                         `Fertility` = "Standard female fertility\ndrive",
                                         `Recessivelethal` = "Both-sex viability drive",
                                         `X-linked-fertility` = 'X-linked female fertility\ndrive', 
                                         `Fertility Fertility` = "Standard female fertility\ndrive with distant-site\nfemale fertility",
                                         `Haplolethal Fertility` = 'Haplolethal rescue\ndrive with distant-site\nfemale fertility', 
                                         `Haplolethal Recessivelethal` = 'Haplolethal rescue\ndrive with distant-site\nviability', 
                                         `Haplolethal X-linked-fertility` = 'Haplolethal rescue\ndrive with distant-site\nX-linked female fertility', 
                                         `Recessivelethal Fertility` = 'Haplosufficient rescue\ndrive with distant-site\nfemale fertility',
                                         `Recessivelethal Recessivelethal` = 'Haplosufficient rescue\ndrive with distant-site\nviability', 
                                         `Recessivelethal X-linked-fertility` = 'Haplosufficient rescue\ndrive with distant-site\nX-linked female fertility')

model_output_gen110 <- model_output_labels[model_output_labels$generation==110,] %>%
  mutate(suppressed = case_when(popSize == 0 ~ 1,
                                popSize > 0 ~ 0)) %>%
  group_by(DRIVE_TYPE, R1_OCCURRENCE_RATE) %>%
  summarise(suppressionRate = sum(suppressed)/50)
model_output_gen110$R1_OCCURRENCE_RATE[model_output_gen110$R1_OCCURRENCE_RATE == 0] <- 10^-7

p5 <- ggplot(data = model_output_gen110) +
  facet_wrap(DRIVE_TYPE ~ ., nrow = 2) +
  geom_line(aes(x = log(R1_OCCURRENCE_RATE, 10), y = suppressionRate, colour = DRIVE_TYPE), na.rm = TRUE) +
  scale_color_manual(values = met.brewer("Hiroshige", 10), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  scale_x_continuous(breaks = seq(-7, -1, 1), 
                     labels = c("-Inf", seq(-6, -1, 1))) +
  xlab("Log(r1 formation rate)") +
  ylab("Complete suppression\nsuccess rate") +
  PaperTheme + theme(legend.position = "none")
p5

ggsave(plot = p5, filename = "FigS8.png", height = 10, width = 20, unit = "cm", dpi = 1000)

p6 <- ggplot(data = model_output_labels) +
  facet_grid(DRIVE_TYPE ~ R1_OCCURRENCE_RATE) +
  geom_line(aes(x = generation, y = popSize, group = interaction(DRIVE_TYPE, repetition), colour = DRIVE_TYPE), alpha = 0.4, linewidth = 0.4) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 10), name = "Drive type") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ylab("Population size") +
  ggtitle("Horizontal = Relative functional resistance rate") +
  PaperTheme + theme(legend.position = "none") + theme(strip.text.y = element_text(angle = 360))
p6

ggsave(plot = p6, filename = "FigS8.5.png", height = 18, width = 20, unit = "cm", dpi = 600)
