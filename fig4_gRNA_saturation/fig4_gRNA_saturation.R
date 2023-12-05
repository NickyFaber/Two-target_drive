# Created by Nicky Faber, 2022.
# Adapted from SLiM-Extras at https://github.com/MesserLab/SLiM-Extras.

##############################
## ---- Working directory ----
##############################

rm(list = ls())
setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/02 Evolutionary stability/HaploLethalFertilityDrive/fig4_gRNA_saturation")

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
        plot.title=element_text(size=12, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12))

#############################
## ---- Parameter inputs ----
#############################

input <- list()

input$REPETITION <- 1:10
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
input$GRNA_SATURATION_SIMULATED <- c("F","T")
input$TOTAL_CUT_RATE <- c(0.5,0.9,0.95,0.99)
input$GLOBAL_SATURATION_FACTOR <- seq(1,5.5,0.5)

inputs <- expand.grid(input)

inputs <- inputs[(inputs$GRNA_SATURATION_SIMULATED == "F" & inputs$GLOBAL_SATURATION_FACTOR == 5.5) | 
                   inputs$GRNA_SATURATION_SIMULATED == "T" & inputs$GLOBAL_SATURATION_FACTOR != 5.5,]

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
               sprintf("-d \"GRNA_SATURATION_SIMULATED=%s\"", inputs[i,4]),
               sprintf('-d TOTAL_CUT_RATE=%f', inputs[i,5]),
               sprintf("-d GLOBAL_SATURATION_FACTOR=%f", inputs[i,6]),
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig4_gRNA_saturation.Rdata")

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

save(model_output, file = "Fig4_gRNA_saturation.Rdata")
# load(file = "Fig4_gRNA_saturation.Rdata")

#########################################
########## Plot plots ###################
#########################################

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(GLOBAL_SATURATION_FACTOR %in% c(1:5,5.5))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue with distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue with distant-site')

model_output_plot <- select(main_drives, -GRNA_SATURATION_SIMULATED)
model_output_plot$GLOBAL_SATURATION_FACTOR <- factor(model_output_plot$GLOBAL_SATURATION_FACTOR,
                                                levels = seq(0.5, 5.5, 0.5),
                                                labels = c(seq(0.5, 5, 0.5), "None"))

p1 <- ggplot(data = model_output_plot) +
  facet_grid(TOTAL_CUT_RATE ~ GLOBAL_SATURATION_FACTOR) +
  geom_line(aes(x = generation, y = popSize, group = interaction(DRIVE_TYPE, repetition), colour = DRIVE_TYPE), alpha = 0.3, linewidth = 0.3) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 3), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ylab("Population size") +
  ggtitle("Horizontal = Saturation factor\nVertical = Total cut rate") +
  PaperTheme
p1

ggsave(plot = p1, filename = "FigS9.pdf", height = 15, width = 20, unit = "cm", dpi = 600)

model_output_subset <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(TOTAL_CUT_RATE == 0.9, GLOBAL_SATURATION_FACTOR %in% c(1:5,5.5))

model_output_subset$DRIVE_TYPE <- recode(model_output_subset$DRIVE_TYPE, `Fertility` = 'Standard suppression\ndrive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')
model_output_subset <- select(model_output_subset, -GRNA_SATURATION_SIMULATED)
model_output_subset$GLOBAL_SATURATION_FACTOR <- factor(model_output_subset$GLOBAL_SATURATION_FACTOR,
                                                     levels = seq(0.5, 5.5, 0.5),
                                                     labels = c(seq(0.5, 5, 0.5), "None"))
model_output_subset_drive <- pivot_longer(model_output_subset, cols = c(freq_wt, freq_dr, freq_r2),
                                          names_to = "allele", values_to = "freq")
model_output_subset_drive$allele <- factor(model_output_subset_drive$allele,
                                           levels = c("freq_wt", "freq_dr", "freq_r2"),
                                           labels = c("Wild type","Drive","Non-functional (r2)"))

p2 <- ggplot(data = model_output_subset_drive) +
  facet_grid(DRIVE_TYPE ~ GLOBAL_SATURATION_FACTOR) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, size = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 3), name = "Allele") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  PaperTheme +
  ggtitle("Horizontal = Saturation factor\nTotal cut rate = 0.9\nHoming site")
p2

model_output_subset_offtarget <- pivot_longer(model_output_subset, cols = c("freq_offtarget_wt", "freq_offtarget_disrupted"),
                                              names_to = "allele", values_to = "freq")

p3 <- ggplot(data = model_output_subset_offtarget) +
  facet_grid(DRIVE_TYPE ~ GLOBAL_SATURATION_FACTOR) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, size = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 2), name = "Allele") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  PaperTheme
p3

p <- (p2 / p3) + plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "FigS10.png", height = 30, width = 20, unit = "cm", dpi = 600)



