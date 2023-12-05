# Created by Nicky Faber, 2022.
# Adapted from SLiM-Extras at https://github.com/MesserLab/SLiM-Extras.

##############################
## ---- Working directory ----
##############################

rm(list = ls())
setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/02 Evolutionary stability/HaploLethalFertilityDrive/fig2_embryo_somatic_cutting")

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
input$EMBRYO_RESISTANCE_CUT_RATE <- seq(0,1,0.1)
input$SOMATIC_EXPRESSION_FITNESS_MULTIPLIER <- seq(0,1,0.1)
input$INTRO_FREQ <- 0.1

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
               sprintf('-d EMBRYO_RESISTANCE_CUT_RATE=%f', inputs[i,4]),
               sprintf('-d SOMATIC_EXPRESSION_FITNESS_MULTIPLIER=%f', inputs[i,5]),
               sprintf('-d INTRO_FREQ=%f', inputs[i,6]),
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig2_embryo_somatic_cutting.Rdata")

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

save(model_output, file = "Fig2_embryo_somatic_cutting.Rdata")
# load(file = "Fig2_embryo_somatic_cutting.Rdata")

#########################################
########## Plot plots ###################
#########################################

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(EMBRYO_RESISTANCE_CUT_RATE %in% seq(0,0.8,0.4) & SOMATIC_EXPRESSION_FITNESS_MULTIPLIER %in% seq(0.4,1,0.3))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue with distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue with distant-site')

p1 <- ggplot(data = main_drives) +
  facet_grid(EMBRYO_RESISTANCE_CUT_RATE ~ SOMATIC_EXPRESSION_FITNESS_MULTIPLIER) +
  geom_line(aes(x = generation, y = popSize, group = interaction(DRIVE_TYPE, repetition), colour = DRIVE_TYPE), alpha = 0.4, linewidth = 0.4) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 3), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  scale_y_continuous(breaks = seq(0,100000,25000)) +
  xlab("Generation") +
  ylab("Population size") +
  ggtitle("Horizontal = Somatic expression fertility effect\nVertical = Embryo cut rate") +
  PaperTheme
p1

ggsave(plot = p1, filename = "Fig4.pdf", height = 12, width = 10, unit = "cm", dpi = 1000)

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility"))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')

model_output_gen110 <- main_drives[main_drives$generation==110,] %>%
  mutate(suppressed = case_when(popSize == 0 ~ 1,
                                popSize > 0 ~ 0)) %>%
  group_by(DRIVE_TYPE, SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, EMBRYO_RESISTANCE_CUT_RATE) %>%
  summarise(suppressionRate = sum(suppressed)/10)

p2 <- ggplot(data = model_output_gen110) +
  facet_grid(. ~ DRIVE_TYPE) +
  geom_raster(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, y = EMBRYO_RESISTANCE_CUT_RATE, fill = suppressionRate)) +
  scale_fill_gradientn(colors=met.brewer("Greek"), limits = c(0,1), name = "Suppression success rate") +
  xlab("Somatic expression fitness multiplier") +
  ylab("Embryo cut rate") +
  PaperTheme
p2

model_output_gen100 <- main_drives[main_drives$generation %in% 100:109,] %>%
  group_by(DRIVE_TYPE, SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, EMBRYO_RESISTANCE_CUT_RATE) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genetic_load_mean = mean(genetic_load, na.rm = TRUE),
            genetic_load_sd = sd(genetic_load, na.rm = TRUE),
            genetic_load_n = sum(!is.na(genetic_load))) %>%
  mutate(genetic_load_se = genetic_load_sd / sqrt(genetic_load_n),
         genetic_load_lower_ci = genetic_load_mean - qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se,
         genetic_load_upper_ci = genetic_load_mean + qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se)

p3 <- ggplot(data = model_output_gen100) +
  facet_grid(. ~ DRIVE_TYPE) +
  geom_raster(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, y = EMBRYO_RESISTANCE_CUT_RATE, fill = genetic_load_mean)) +
  geom_point(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER - 0.025, y = EMBRYO_RESISTANCE_CUT_RATE, color=genetic_load_lower_ci), size = 1) +
  geom_point(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER + 0.025, y = EMBRYO_RESISTANCE_CUT_RATE, color=genetic_load_upper_ci), size = 1) +
  geom_text(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, y = EMBRYO_RESISTANCE_CUT_RATE, label=genetic_load_n), size = 1.5) +
  scale_fill_gradient(high = met.brewer("Hiroshige",3)[3], name = "Mean genetic load") +
  scale_colour_gradient(high = met.brewer("Hiroshige",3)[3], name = "95% confidence interval genetic load") +
  xlab("Somatic expression fitness multiplier") +
  ylab("Embryo cut rate") +
  PaperTheme + guides(color = "none")
p3

# load(file = "Fig2_embryo_somatic_cutting_GL_plot.Rdata")

p <- (p2 / p3) + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & theme(legend.position = 'bottom'); p

ggsave(plot = p, filename = "Fig5.pdf", height = 18, width = 18, unit = "cm", dpi = 1000)

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(EMBRYO_RESISTANCE_CUT_RATE == 0.4, SOMATIC_EXPRESSION_FITNESS_MULTIPLIER > 0.4)

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression\ndrive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')
model_output_subset_drive <- pivot_longer(main_drives, cols = c(freq_wt, freq_dr, freq_r2),
                                          names_to = "allele", values_to = "freq")
model_output_subset_drive$allele <- factor(model_output_subset_drive$allele,
                                           levels = c("freq_wt", "freq_dr", "freq_r2"),
                                           labels = c("Wild type","Drive","Non-functional (r2)"))

p4 <- ggplot(data = model_output_subset_drive) +
  facet_grid(DRIVE_TYPE ~ SOMATIC_EXPRESSION_FITNESS_MULTIPLIER) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = rev(met.brewer("Hiroshige", 3)), name = "Allele") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Horizontal = Somatic expression fitness multiplier\nEmbryo cut rate = 0.4\nHoming site") +
  PaperTheme
p4

model_output_subset_offtarget <- pivot_longer(main_drives, cols = c("freq_offtarget_wt", "freq_offtarget_disrupted"),
                                       names_to = "allele", values_to = "freq")
model_output_subset_offtarget$allele <- factor(model_output_subset_offtarget$allele,
                                               labels = c("Wild type","Non-functional (r2)"))

p5 <- ggplot(data = model_output_subset_offtarget) +
  facet_grid(DRIVE_TYPE ~ SOMATIC_EXPRESSION_FITNESS_MULTIPLIER) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 2), name = "Allele") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Distant site") +
  PaperTheme + theme(legend.position = "none")
p5

p <- (p4 / p5) + plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "FigS4.png", height = 30, width = 20, unit = "cm", dpi = 600)

model_output_labels <- model_output
model_output_labels$DRIVE_TYPE <- recode(model_output_labels$DRIVE_TYPE, 
                                         `Fertility` = "Standard female fertility drive",
                                         `Recessivelethal` = "Both-sex viability drive",
                                         `X-linked-fertility` = 'X-linked female fertility drive', 
                                         `Fertility Fertility` = "Standard female fertility drive\nwith distant-site\nfemale fertility",
                                         `Haplolethal Fertility` = 'Haplolethal rescue\ndrive with distant-site\nfemale fertility', 
                                         `Haplolethal Recessivelethal` = 'Haplolethal rescue\ndrive with distant-site\nviability', 
                                         `Haplolethal X-linked-fertility` = 'Haplolethal rescue\ndrive with distant-site\nX-linked female fertility', 
                                         `Recessivelethal Fertility` = 'Haplosufficient\nrescue drive with distant-site\nfemale fertility',
                                         `Recessivelethal Recessivelethal` = 'Haplosufficient\nrescue drive with distant-site\nviability', 
                                         `Recessivelethal X-linked-fertility` = 'Haplosufficient\nrescue drive with distant-site\nX-linked female fertility')

model_output_gen110 <- model_output_labels[model_output_labels$generation==110,] %>%
  mutate(suppressed = case_when(popSize == 0 ~ 1,
                                popSize > 0 ~ 0)) %>%
  group_by(DRIVE_TYPE, SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, EMBRYO_RESISTANCE_CUT_RATE) %>%
  summarise(suppressionRate = sum(suppressed)/10)

p6 <- ggplot(data = model_output_gen110) +
  facet_grid(DRIVE_TYPE ~ .) +
  geom_raster(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, y = EMBRYO_RESISTANCE_CUT_RATE, fill = suppressionRate)) +
  scale_fill_gradientn(colors=met.brewer("Greek"), limits = c(0,1), name = "Suppression success\nrate") + 
  guides(fill = guide_colourbar(title.position = "top")) +
  xlab("Somatic expression\nfitness multiplier") +
  ylab("Embryo cut rate") +
  PaperTheme + theme(strip.text.y = element_blank())
p6

model_output_gen100 <- model_output_labels[model_output_labels$generation %in% 100:109,] %>%
  group_by(DRIVE_TYPE, SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, EMBRYO_RESISTANCE_CUT_RATE) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genetic_load_mean = mean(genetic_load, na.rm = TRUE),
            genetic_load_sd = sd(genetic_load, na.rm = TRUE),
            genetic_load_n = sum(!is.na(genetic_load))) %>%
  mutate(genetic_load_se = genetic_load_sd / sqrt(genetic_load_n),
         genetic_load_lower_ci = genetic_load_mean - qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se,
         genetic_load_upper_ci = genetic_load_mean + qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se)

p7 <- ggplot(data = model_output_gen100) +
  facet_grid(DRIVE_TYPE ~ .) +
  geom_raster(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, y = EMBRYO_RESISTANCE_CUT_RATE, fill = genetic_load_mean)) +
  geom_point(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER - 0.025, y = EMBRYO_RESISTANCE_CUT_RATE, color=genetic_load_lower_ci), size = 1) +
  geom_point(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER + 0.025, y = EMBRYO_RESISTANCE_CUT_RATE, color=genetic_load_upper_ci), size = 1) +
  geom_text(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, y = EMBRYO_RESISTANCE_CUT_RATE, label=genetic_load_n), size = 1.5) +
  scale_fill_gradient(high = met.brewer("Hiroshige",3)[3], name = "Mean genetic load") +
  scale_colour_gradient(high = met.brewer("Hiroshige",3)[3], name = "95% confidence interval genetic load") +
  guides(fill = guide_colourbar(title.position = "top")) +
  xlab("Somatic expression\nfitness multiplier") +
  ylab("Embryo cut rate") +
  PaperTheme + guides(color = "none") + theme(strip.text.y = element_blank(), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank(),
                                              axis.title.y=element_blank())
p7

# load(file = "Fig2_embryo_somatic_cutting_GL_plot.Rdata")

model_output_last_gen <- model_output_labels %>%
  filter(popSize != 0) %>%
  group_by(repetition) %>%
  filter(generation == max(generation)) %>%
  ungroup() %>%
  group_by(DRIVE_TYPE, SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, EMBRYO_RESISTANCE_CUT_RATE) %>%
  summarise(generation_mean = mean(generation, na.rm = TRUE),
            generation_sd = sd(generation, na.rm = TRUE),
            generation_n = sum(!is.na(generation))) %>%
  mutate(generation_se = generation_sd / sqrt(generation_n),
         generation_lower_ci = generation_mean - qt(1 - (0.05 / 2), generation_n - 1) * generation_se,
         generation_upper_ci = generation_mean + qt(1 - (0.05 / 2), generation_n - 1) * generation_se)

p8 <- ggplot(data = model_output_last_gen) +
  facet_grid(DRIVE_TYPE ~ .) +
  geom_raster(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER, y = EMBRYO_RESISTANCE_CUT_RATE, fill = generation_mean)) +
  geom_point(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER - 0.025, y = EMBRYO_RESISTANCE_CUT_RATE, color=generation_lower_ci), size = 1) +
  geom_point(aes(x = SOMATIC_EXPRESSION_FITNESS_MULTIPLIER + 0.025, y = EMBRYO_RESISTANCE_CUT_RATE, color=generation_upper_ci), size = 1) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3"), limits = c(0,130), name = "Mean last generation") +
  scale_colour_gradientn(colors=met.brewer("VanGogh3"), limits = c(0,130), name = "95% confidence interval last generation") +
  guides(fill = guide_colourbar(title.position = "top")) +
  xlab("Somatic expression\nfitness multiplier") +
  ylab("Embryo cut rate") +
  PaperTheme + guides(color = "none") + theme(strip.text.y = element_text(angle = 360), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank(),
                                              axis.title.y=element_blank())
p8

p <- (p6 | p7 | p8) + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'keep'); p

ggsave(plot = p, filename = "FigS5.pdf", height = 30, width = 20, unit = "cm", dpi = 600)











