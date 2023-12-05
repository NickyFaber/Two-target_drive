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
input$GENERATIONS <- 161
input$CAPACITY <- 100000
input$INTRO_FREQ <- 0.5
input$LOW_DENSITY_GROWTH_RATE <- 100
input$MAX_OFFSPRING <- 500
input$GL_RUN <- "F"
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
               sprintf("-d GENERATIONS=%s", inputs[i,2]),
               sprintf("-d CAPACITY=%s", inputs[i,3]),
               sprintf("-d INTRO_FREQ=%s", inputs[i,4]),
               sprintf("-d LOW_DENSITY_GROWTH_RATE=%s", inputs[i,5]),
               sprintf("-d MAX_OFFSPRING=%s", inputs[i,6]),
               sprintf("-d GL_RUN=%s", inputs[i,7]),
               sprintf("-d \"DRIVE_TYPE='%s'\"", inputs[i,8]),
               sprintf("-d \"GRNA_SATURATION_SIMULATED=%s\"", inputs[i,9]),
               sprintf('-d TOTAL_CUT_RATE=%f', inputs[i,10]),
               sprintf("-d GLOBAL_SATURATION_FACTOR=%f", inputs[i,11]),
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig4_gRNA_saturation_GL_equilibria.Rdata")

##############################################
########## Parse and calculate stuff #########
##############################################

model_output <- parse_slim_output(raw_slim_output_matrix, generations = 161, inputs)

model_output <- model_output %>%
  group_by(DRIVE_TYPE, repetition) %>%
  mutate(genetic_load = 1 - lead(popSize) * bonus_pop_factor / expectedPopSize)

#########################################
########## Save model ###################
#########################################

save(model_output, file = "Fig4_gRNA_saturation_GL_equilibria.Rdata")
# load(file = "Fig4_gRNA_saturation_GL_equilibria.Rdata")

###################################################################
## ---- Find suppressed populations and rerun with GL settings ---- 
###################################################################

run_with_GL <- model_output[model_output$generation == 160,] %>%
  group_by(DRIVE_TYPE, GRNA_SATURATION_SIMULATED, TOTAL_CUT_RATE, GLOBAL_SATURATION_FACTOR) %>%
  summarise(nonSuppressed = sum(popSize > 0)) %>%
  filter(nonSuppressed < 10) %>%
  select(DRIVE_TYPE:GLOBAL_SATURATION_FACTOR)

inputs <- run_with_GL %>%
  mutate(REPETITION = 10,
         GENERATIONS = 161,
         CAPACITY = 100000,
         INTRO_FREQ = 0.5,
         LOW_DENSITY_GROWTH_RATE = 100,
         MAX_OFFSPRING = 500,
         GL_RUN = "T") %>%
  uncount(REPETITION) %>%
  mutate(REPETITION = 1) %>%
  relocate(DRIVE_TYPE:GLOBAL_SATURATION_FACTOR, .after = last_col()) %>%
  relocate(REPETITION, .before = GENERATIONS)

inputs$REPETITION <- rep(1:10, times = nrow(inputs)/10)
inputs$DRIVE_TYPE <- as.character(inputs$DRIVE_TYPE)
inputs$DRIVE_TYPE <- case_when(
  inputs$DRIVE_TYPE == "Fertility" ~ "fertility",
  inputs$DRIVE_TYPE == "X-linked-fertility" ~ "Xlinkedfertility",
  inputs$DRIVE_TYPE == "Recessivelethal" ~ "recessivelethal",
  inputs$DRIVE_TYPE == "Fertility Fertility" ~ "fertility_fertility",
  inputs$DRIVE_TYPE == "Haplolethal Fertility" ~ "haplolethal_fertility",
  inputs$DRIVE_TYPE == "Haplolethal Recessivelethal" ~ "haplolethal_recessivelethal",
  inputs$DRIVE_TYPE == "Haplolethal X-linked-fertility" ~ "haplolethal_Xlinkedfertility",
  inputs$DRIVE_TYPE == "Recessivelethal Fertility" ~ "recessivelethal_fertility",
  inputs$DRIVE_TYPE == "Recessivelethal Recessivelethal" ~ "recessivelethal_recessivelethal",
  inputs$DRIVE_TYPE == "Recessivelethal X-linked-fertility" ~ "recessivelethal_Xlinkedfertility",
)

#####################################
## ---- Run parallel simulations ---- 
#####################################

cl <- makeCluster(8, outfile="")
registerDoParallel(cl)

# Run SLiM in parallel:
raw_slim_output_matrix <- foreach(i=iter(1:nrow(inputs), by='row')) %dopar% {
  cat(sprintf("Working on iteration %i out of %i\n", i, nrow(inputs)))
  system(paste("slim",
               sprintf("-d GENERATIONS=%s", inputs[i,2]),
               sprintf("-d CAPACITY=%s", inputs[i,3]),
               sprintf("-d INTRO_FREQ=%s", inputs[i,4]),
               sprintf("-d LOW_DENSITY_GROWTH_RATE=%s", inputs[i,5]),
               sprintf("-d MAX_OFFSPRING=%s", inputs[i,6]),
               sprintf("-d GL_RUN=%s", inputs[i,7]),
               sprintf("-d \"DRIVE_TYPE='%s'\"", inputs[i,8]),
               sprintf("-d \"GRNA_SATURATION_SIMULATED=%s\"", inputs[i,9]),
               sprintf('-d TOTAL_CUT_RATE=%f', inputs[i,10]),
               sprintf("-d GLOBAL_SATURATION_FACTOR=%f", inputs[i,11]),
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig4_gRNA_saturation_GL.Rdata")

##############################################
########## Parse and calculate stuff #########
##############################################

model_output <- parse_slim_output(raw_slim_output_matrix, generations = 161, inputs)

model_output <- model_output %>%
  group_by(DRIVE_TYPE, repetition) %>%
  mutate(genetic_load = 1 - lead(popSize) * bonus_pop_factor / expectedPopSize)

#########################################
########## Save model ###################
#########################################

save(raw_slim_output_matrix, file = "Fig4_gRNA_saturation_GL.Rdata")
# load(file = "Fig4_gRNA_saturation_GL.Rdata")

#########################################
########## Plot plots ###################
#########################################

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility")) %>%
  filter(TOTAL_CUT_RATE %in% c(0.5,0.9,0.95,0.99))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard supression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue with distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue with distant-site')

model_output_gen109 <- main_drives[main_drives$generation %in% 150:159,] %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, GLOBAL_SATURATION_FACTOR) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genLoad_mean = mean(genetic_load, na.rm = TRUE),
            genLoad_sd = sd(genetic_load, na.rm = TRUE),
            genLoad_n = sum(!is.na(genetic_load))) %>%
  mutate(genLoad_se = genLoad_sd / sqrt(genLoad_n),
         genLoad_lower_ci = genLoad_mean - qt(1 - (0.05 / 2), genLoad_n - 1) * genLoad_se,
         genLoad_upper_ci = genLoad_mean + qt(1 - (0.05 / 2), genLoad_n - 1) * genLoad_se)

model_output_gen109_noSat <- filter(model_output_gen109, GLOBAL_SATURATION_FACTOR == 5.5)
model_output_gen109_SatOnly <- filter(model_output_gen109, GLOBAL_SATURATION_FACTOR < 5.5)

p1 <- ggplot(data = model_output_gen109_SatOnly, aes(x = GLOBAL_SATURATION_FACTOR, y = genLoad_mean, group = DRIVE_TYPE)) +
  facet_grid(. ~ TOTAL_CUT_RATE) +
  geom_line(position=position_dodge(width=0.2), aes(colour = DRIVE_TYPE)) +
  geom_errorbar(position=position_dodge(width=0.2), linewidth = 0.125, width = 0.1, colour = "black", aes(ymin = genLoad_lower_ci, ymax = genLoad_upper_ci)) +
  geom_point(data = model_output_gen109_noSat, aes(x = GLOBAL_SATURATION_FACTOR, y = genLoad_mean, colour = DRIVE_TYPE), position=position_dodge(width=0.2), na.rm = TRUE) +
  geom_errorbar(data = model_output_gen109_noSat, aes(x = GLOBAL_SATURATION_FACTOR, y = genLoad_mean, ymin = genLoad_lower_ci, ymax = genLoad_upper_ci), 
                size = 0.125, width = 0.1, colour = "black", position=position_dodge(width=0.2)) +
  scale_color_manual(values = met.brewer("Cross", 3), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 3, title.position = "top", title.hjust = 0.5, override.aes = list(alpha = 1))) +
  xlab("gRNA saturation factor") +
  scale_x_continuous(breaks = c(1:5,5.5), 
                     labels = c(1:5, "None")) +
  ylab("Genetic load") + ylim(c(0,1)) +
  ggtitle("Horizontal = Total cut rate") +
  PaperTheme + theme(axis.text.x = element_text(angle = 45, hjust=1))
p1

ggsave(plot = p1, filename = "Fig7.png", height = 9, width = 18, unit = "cm", dpi = 600)

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

model_output_gen109 <- model_output_labels[model_output_labels$generation %in% 150:159,] %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, GLOBAL_SATURATION_FACTOR) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genLoad_mean = mean(genetic_load, na.rm = TRUE),
            genLoad_sd = sd(genetic_load, na.rm = TRUE),
            genLoad_n = sum(!is.na(genetic_load))) %>%
  mutate(genLoad_se = genLoad_sd / sqrt(genLoad_n),
         genLoad_lower_ci = genLoad_mean - qt(1 - (0.05 / 2), genLoad_n - 1) * genLoad_se,
         genLoad_upper_ci = genLoad_mean + qt(1 - (0.05 / 2), genLoad_n - 1) * genLoad_se)

model_output_gen109_noSat <- filter(model_output_gen109, GLOBAL_SATURATION_FACTOR == 5.5)
model_output_gen109_SatOnly <- filter(model_output_gen109, GLOBAL_SATURATION_FACTOR < 5.5)

p2 <- ggplot(data = model_output_gen109_SatOnly, aes(x = GLOBAL_SATURATION_FACTOR, y = genLoad_mean, group = DRIVE_TYPE)) +
  facet_grid(DRIVE_TYPE ~ TOTAL_CUT_RATE) +
  geom_line(position=position_dodge(width=0.2), aes(colour = DRIVE_TYPE)) +
  geom_errorbar(position=position_dodge(width=0.2), linewidth = 0.125, width = 0.1, colour = "black", aes(ymin = genLoad_lower_ci, ymax = genLoad_upper_ci)) +
  geom_point(data = model_output_gen109_noSat, aes(x = GLOBAL_SATURATION_FACTOR, y = genLoad_mean, colour = DRIVE_TYPE), position=position_dodge(width=0.2), na.rm = TRUE) +
  geom_errorbar(data = model_output_gen109_noSat, aes(x = GLOBAL_SATURATION_FACTOR, y = genLoad_mean, ymin = genLoad_lower_ci, ymax = genLoad_upper_ci), 
                size = 0.125, width = 0.1, colour = "black", position=position_dodge(width=0.2)) +
  scale_color_manual(values = met.brewer("Hiroshige", 10), name = "Female fertility suppression drive type") +
  guides(colour = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, override.aes = list(alpha = 1))) +
  xlab("gRNA saturation factor") +
  scale_x_continuous(breaks = c(1:5,5.5), 
                     labels = c(1:5, "None")) +
  ylab("Genetic load") + ylim(c(0,1)) +
  ggtitle("Horizontal = Total cut rate") +
  PaperTheme + guides(color = "none") + theme(strip.text.y = element_text(angle = 360),
                                              axis.text.x = element_text(angle = 45, hjust=1))
p2

ggsave(plot = p2, filename = "FigS11.png", height = 30, width = 20, unit = "cm", dpi = 600)






