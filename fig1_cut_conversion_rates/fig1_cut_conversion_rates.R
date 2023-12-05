# Created by Nicky Faber, 2022.
# Adapted from SLiM-Extras at https://github.com/MesserLab/SLiM-Extras.

##############################
## ---- Working directory ----
##############################

rm(list = ls())
setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/02 Evolutionary stability/HaploLethalFertilityDrive/fig1_cut_conversion_rates")

library('foreach', quietly = TRUE)
library('doParallel', quietly = TRUE)
library('future', quietly = TRUE)
library('tidyverse', quietly = TRUE)
library('dplyr', quietly = TRUE)
library('ggplot2', quietly = TRUE)
library('ggh4x', quietly = TRUE)
library('MetBrewer', quietly = TRUE)
library('vcfR', quietly = TRUE)
library('patchwork', quietly = TRUE)
library('svglite', quietly = TRUE)
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

input$DRIVE_CONVERSION_RATE <- seq(0,1,0.1)
input$TOTAL_CUT_RATE <- seq(0,1,0.1)

inputs <- expand.grid(input)

inputs <- inputs[inputs$DRIVE_CONVERSION_RATE <= inputs$TOTAL_CUT_RATE,]

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
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig1_cut_conversion_rates.Rdata")

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

save(model_output, file = "Fig1_cut_conversion_rates.Rdata")
# load(file = "Fig1_cut_conversion_rates.Rdata")

#########################################
########## Plot plots ###################
#########################################

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(TOTAL_CUT_RATE %in% seq(0.6,1,0.2) & DRIVE_CONVERSION_RATE %in% c(0.2,0.4,0.6,0.8,1))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue with distant-site', 
                                  `Recessivelethal Fertility` = 'Haplosufficient rescue with distant-site')

design <- "
  ABC
  DEF
  GHI
  #JK
  ##L
"

p1 <- ggplot(data = main_drives) +
  facet_manual(DRIVE_CONVERSION_RATE ~ TOTAL_CUT_RATE, design = design, strip = strip_split(position = c("right", "top"))) +
  geom_line(aes(x = generation, y = popSize, group = interaction(DRIVE_TYPE, repetition), colour = DRIVE_TYPE), alpha = 0.4, linewidth = 0.4) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige",3), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ylab("Population size") +
  ggtitle("Horizontal = Total cut rate\nVertical = Drive conversion rate") +
  PaperTheme
p1

ggsave(plot = p1, filename = "Fig2.pdf", height = 18, width = 10, unit = "cm", dpi = 1000)

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility"))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')

model_output_gen110 <- main_drives[main_drives$generation==110,] %>%
  mutate(suppressed = case_when(popSize == 0 ~ 1,
                                popSize > 0 ~ 0)) %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, DRIVE_CONVERSION_RATE) %>%
  summarise(suppressionRate = sum(suppressed)/10)

p2 <- ggplot(data = model_output_gen110) +
  facet_grid(. ~ DRIVE_TYPE) +
  geom_raster(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, fill = suppressionRate)) +
  scale_fill_gradientn(colors=met.brewer("Greek"), limits = c(0,1), name = "Suppression success rate") +
  xlab("Total cut rate") +
  ylab("Drive conversion rate") +
  PaperTheme
p2

model_output_gen100 <- main_drives[main_drives$generation %in% 100:109,] %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, DRIVE_CONVERSION_RATE) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genetic_load_mean = mean(genetic_load, na.rm = TRUE),
            genetic_load_sd = sd(genetic_load, na.rm = TRUE),
            genetic_load_n = sum(!is.na(genetic_load))) %>%
  mutate(genetic_load_se = genetic_load_sd / sqrt(genetic_load_n),
         genetic_load_lower_ci = genetic_load_mean - abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se),
         genetic_load_upper_ci = genetic_load_mean + abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se))

p3 <- ggplot(data = model_output_gen100) +
  facet_grid(. ~ DRIVE_TYPE) +
  geom_raster(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, fill = genetic_load_mean)) +
  geom_point(aes(x = TOTAL_CUT_RATE - 0.025, y = DRIVE_CONVERSION_RATE, color=genetic_load_lower_ci), size = 1) +
  geom_point(aes(x = TOTAL_CUT_RATE + 0.025, y = DRIVE_CONVERSION_RATE, color=genetic_load_upper_ci), size = 1) +
  geom_text(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, label=genetic_load_n), size = 1.5) +
  scale_fill_gradient(high = met.brewer("Hiroshige",3)[3], name = "Mean genetic load") +
  scale_colour_gradient(high = met.brewer("Hiroshige",3)[3], name = "95% confidence interval genetic load") +
  xlab("Total cut rate") +
  ylab("Drive conversion rate") +
  PaperTheme + guides(color = "none")
p3

# load(file = "Fig1_cut_conversion_rates_GL_plot.Rdata")

p <- (p2 / p3) + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & theme(legend.position = 'bottom'); p

ggsave(plot = p, filename = "Fig3.pdf", height = 18, width = 15, unit = "cm", dpi = 1000)

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(DRIVE_CONVERSION_RATE == 0.6)

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression\ndrive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')
model_output_subset <- main_drives[main_drives$DRIVE_CONVERSION_RATE==0.6,]
model_output_subset_drive <- pivot_longer(model_output_subset, cols = c(freq_wt, freq_dr, freq_r2),
                                          names_to = "allele", values_to = "freq")
model_output_subset_drive$allele <- factor(model_output_subset_drive$allele,
                                           levels = c("freq_wt", "freq_dr", "freq_r2"),
                                           labels = c("Wild type","Drive","Non-functional (r2)"))

p4 <- ggplot(data = model_output_subset_drive) +
  facet_grid(DRIVE_TYPE ~ TOTAL_CUT_RATE) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = rev(met.brewer("Hiroshige", 3)), name = "Allele") +
  guides(colour = guide_legend(override.aes = list(alpha = 1), title.position = "left")) +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Horizontal = Total cut rate\nDrive conversion rate = 0.6\nHoming site") +
  PaperTheme + theme(legend.position = "bottom")
p4

model_output_subset_offtarget <- pivot_longer(model_output_subset, cols = c("freq_offtarget_wt", "freq_offtarget_disrupted"),
                                              names_to = "allele", values_to = "freq")
model_output_subset_offtarget$allele <- factor(model_output_subset_offtarget$allele,
                                           labels = c("Wild type","Non-functional (r2)"))

p5 <- ggplot(data = model_output_subset_offtarget) +
  facet_grid(DRIVE_TYPE ~ TOTAL_CUT_RATE) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 2), name = "Allele") +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Distant site") +
  PaperTheme + theme(legend.position = "none")
p5

p <- (p4 / p5) + plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "FigS1.png", height = 30, width = 20, unit = "cm", dpi = 600)

model_output_subset$DRIVE_TYPE <- recode(model_output_subset$DRIVE_TYPE, `Standard supression\ndrive` = 'Standard supression drive')
model_output_subset$TOTAL_CUT_RATE <- factor(model_output_subset$TOTAL_CUT_RATE)

p6 <- ggplot(data = model_output_subset) +
  facet_grid(. ~ DRIVE_TYPE) +
  geom_line(aes(x = generation, y = popSize, group = interaction(TOTAL_CUT_RATE, repetition), colour = TOTAL_CUT_RATE), alpha = 0.4, linewidth = 0.4) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 5), name = "Total cut rate") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ylab("Population size") +
  ggtitle("Drive conversion = 0.6") +
  PaperTheme + theme(legend.position = "right")
p6

ggsave(plot = p6, filename = "FigS2.png", height = 8, width = 20, unit = "cm", dpi = 600)



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
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, DRIVE_CONVERSION_RATE) %>%
  summarise(suppressionRate = sum(suppressed)/10)

p7 <- ggplot(data = model_output_gen110) +
  facet_grid(DRIVE_TYPE ~ ., labeller = ) +
  geom_raster(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, fill = suppressionRate)) +
  scale_fill_gradientn(colors=met.brewer("Greek"), limits = c(0,1), name = "Suppression success\nrate")+ 
  guides(fill = guide_colourbar(title.position = "top")) +
  xlab("Total cut rate") +
  ylab("Drive conversion rate") +
  PaperTheme + theme(strip.text.y = element_blank())
p7

model_output_gen100 <- model_output_labels[model_output_labels$generation %in% 100:109,] %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, DRIVE_CONVERSION_RATE) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genetic_load_mean = mean(genetic_load, na.rm = TRUE),
            genetic_load_sd = sd(genetic_load, na.rm = TRUE),
            genetic_load_n = sum(!is.na(genetic_load))) %>%
  mutate(genetic_load_se = genetic_load_sd / sqrt(genetic_load_n),
         genetic_load_lower_ci = genetic_load_mean - abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se),
         genetic_load_upper_ci = genetic_load_mean + abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se))

p8 <- ggplot(data = model_output_gen100) +
  facet_grid(DRIVE_TYPE ~ .) +
  geom_raster(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, fill = genetic_load_mean)) +
  geom_point(aes(x = TOTAL_CUT_RATE - 0.025, y = DRIVE_CONVERSION_RATE, color=genetic_load_lower_ci), size = 1) +
  geom_point(aes(x = TOTAL_CUT_RATE + 0.025, y = DRIVE_CONVERSION_RATE, color=genetic_load_upper_ci), size = 1) +
  geom_text(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, label=genetic_load_n), size = 1.5) +
  scale_fill_gradient(high = met.brewer("Hiroshige",3)[3], name = "Mean genetic load") +
  scale_colour_gradient(high = met.brewer("Hiroshige",3)[3], name = "95% confidence interval genetic load") +
  guides(fill = guide_colourbar(title.position = "top")) +
  xlab("Total cut rate") +
  ylab("Drive conversion rate") +
  PaperTheme + guides(color = "none") + theme(strip.text.y = element_blank(), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank(),
                                              axis.title.y=element_blank())
p8

# load(file = "Fig1_cut_conversion_rates_GL_plot.Rdata")

model_output_last_gen <- model_output_labels %>%
  filter(popSize != 0) %>%
  group_by(repetition) %>%
  filter(generation == max(generation)) %>%
  ungroup() %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, DRIVE_CONVERSION_RATE) %>%
  summarise(generation_mean = mean(generation, na.rm = TRUE),
            generation_sd = sd(generation, na.rm = TRUE),
            generation_n = sum(!is.na(generation))) %>%
  mutate(generation_se = generation_sd / sqrt(generation_n),
         generation_lower_ci = generation_mean - abs(qt(1 - (0.05 / 2), generation_n - 1) * generation_se),
         generation_upper_ci = generation_mean + abs(qt(1 - (0.05 / 2), generation_n - 1) * generation_se))

p9 <- ggplot(data = model_output_last_gen) +
  facet_grid(DRIVE_TYPE ~ .) +
  geom_raster(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, fill = generation_mean)) +
  geom_point(aes(x = TOTAL_CUT_RATE - 0.025, y = DRIVE_CONVERSION_RATE, color=generation_lower_ci), size = 1) +
  geom_point(aes(x = TOTAL_CUT_RATE + 0.025, y = DRIVE_CONVERSION_RATE, color=generation_upper_ci), size = 1) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3"), limits = c(0,130), name = "Mean last generation") +
  scale_colour_gradientn(colors=met.brewer("VanGogh3"), limits = c(0,130), name = "95% confidence interval last generation") +
  guides(fill = guide_colourbar(title.position = "top")) +
  xlab("Total cut rate") +
  ylab("Drive conversion rate") +
  PaperTheme + guides(color = "none") + theme(strip.text.y = element_text(angle = 360), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank(),
                                              axis.title.y=element_blank())
p9

p <- (p7 | p8 | p9) + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'keep'); p

ggsave(plot = p, filename = "FigS3.pdf", height = 30, width = 20, unit = "cm", dpi = 1000)



main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(TOTAL_CUT_RATE > 0.5, DRIVE_CONVERSION_RATE == 0.3)

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression\ndrive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')
model_output_subset_drive <- pivot_longer(main_drives, cols = c(freq_wt, freq_dr, freq_r2),
                                          names_to = "allele", values_to = "freq")
model_output_subset_drive$allele <- factor(model_output_subset_drive$allele,
                                           levels = c("freq_wt", "freq_dr", "freq_r2"),
                                           labels = c("Wild type","Drive","Non-functional (r2)"))

p10 <- ggplot(data = model_output_subset_drive) +
  facet_grid(DRIVE_TYPE ~ TOTAL_CUT_RATE) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = rev(met.brewer("Hiroshige", 3)), name = "Allele") +
  guides(colour = guide_legend(override.aes = list(alpha = 1), title.position = "left")) +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Horizontal = Total cut rate\nDrive conversion rate = 0.3\nHoming site") +
  PaperTheme + theme(legend.position = "bottom")
p10

model_output_subset_offtarget <- pivot_longer(main_drives, cols = c("freq_offtarget_wt", "freq_offtarget_disrupted"),
                                              names_to = "allele", values_to = "freq")
model_output_subset_offtarget$allele <- factor(model_output_subset_offtarget$allele,
                                               labels = c("Wild type","Non-functional (r2)"))

p11 <- ggplot(data = model_output_subset_offtarget) +
  facet_grid(DRIVE_TYPE ~ TOTAL_CUT_RATE) +
  geom_line(aes(x = generation, y = freq, colour = allele, group = interaction(allele, repetition)), alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(aes(xintercept = 10), linetype = "dotted") +
  scale_color_manual(values = met.brewer("Hiroshige", 2), name = "Allele") +
  ylim(c(0,1)) +
  xlab("Generation") +
  ylab("Frequency") +
  ggtitle("Distant site") +
  PaperTheme + theme(legend.position = "none")
p11

p <- (p10 / p11) + plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "FigS2.5.png", height = 30, width = 20, unit = "cm", dpi = 600)

