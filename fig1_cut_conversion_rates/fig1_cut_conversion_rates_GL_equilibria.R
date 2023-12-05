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
               sprintf("-d GENERATIONS=%s", inputs[i,2]),
               sprintf("-d CAPACITY=%s", inputs[i,3]),
               sprintf("-d INTRO_FREQ=%s", inputs[i,4]),
               sprintf("-d LOW_DENSITY_GROWTH_RATE=%s", inputs[i,5]),
               sprintf("-d MAX_OFFSPRING=%s", inputs[i,6]),
               sprintf("-d GL_RUN=%s", inputs[i,7]),
               sprintf("-d \"DRIVE_TYPE='%s'\"", inputs[i,8]),
               sprintf('-d DRIVE_CONVERSION_RATE=%f', inputs[i,9]),
               sprintf('-d TOTAL_CUT_RATE=%f', inputs[i,10]),
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig1_cut_conversion_rates_GL_equilibria.Rdata")

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

save(model_output, file = "Fig1_cut_conversion_rates_GL_equilibria.Rdata")
# load(file = "Fig1_cut_conversion_rates_GL_equilibria.Rdata")

###################################################################
## ---- Find suppressed populations and rerun with GL settings ---- 
###################################################################

run_with_GL <- model_output[model_output$generation == 160,] %>%
  group_by(DRIVE_TYPE, DRIVE_CONVERSION_RATE, TOTAL_CUT_RATE) %>%
  summarise(nonSuppressed = sum(popSize > 0)) %>%
  filter(nonSuppressed < 10) %>%
  select(DRIVE_TYPE:TOTAL_CUT_RATE)

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
  relocate(DRIVE_TYPE:TOTAL_CUT_RATE, .after = last_col()) %>%
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
               sprintf('-d DRIVE_CONVERSION_RATE=%f', inputs[i,9]),
               sprintf('-d TOTAL_CUT_RATE=%f', inputs[i,10]),
               "../haplolethal_suppression_drive.slim"), intern=TRUE)
}

stopCluster(cl)

save(raw_slim_output_matrix, file = "Fig1_cut_conversion_rates_GL.Rdata")

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

save(model_output, file = "Fig1_cut_conversion_rates_GL.Rdata")
# load(file = "Fig1_cut_conversion_rates_GL.Rdata")

#########################################
########## Plot plots ###################
#########################################

load(file = "Fig1_cut_conversion_rates_GL_equilibria.Rdata")

model_output$rerun_with_GL <- ifelse(is.na(match(paste0(model_output$DRIVE_TYPE, model_output$TOTAL_CUT_RATE, model_output$DRIVE_CONVERSION_RATE), 
                                         paste0(run_with_GL$DRIVE_TYPE, run_with_GL$TOTAL_CUT_RATE, run_with_GL$DRIVE_CONVERSION_RATE))),"No", "Yes")
model_output_equilibria <- model_output %>%
  filter(rerun_with_GL == "No") %>%
  select(-rerun_with_GL)

load(file = "Fig1_cut_conversion_rates_GL.Rdata")

model_output <- bind_rows(model_output, model_output_equilibria)

main_drives <- filter(model_output, DRIVE_TYPE %in% c("Fertility",
                                                      "Haplolethal Fertility",
                                                      "Recessivelethal Fertility"))

main_drives$DRIVE_TYPE <- recode(main_drives$DRIVE_TYPE, `Fertility` = 'Standard suppression drive', 
                                 `Haplolethal Fertility` = 'Haplolethal rescue\nwith distant-site', 
                                 `Recessivelethal Fertility` = 'Haplosufficient rescue\nwith distant-site')

model_output_gen150 <- main_drives[main_drives$generation %in% 150:159,] %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, DRIVE_CONVERSION_RATE) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genetic_load_mean = mean(genetic_load, na.rm = TRUE),
            genetic_load_sd = sd(genetic_load, na.rm = TRUE),
            genetic_load_n = sum(!is.na(genetic_load))) %>%
  mutate(genetic_load_se = genetic_load_sd / sqrt(genetic_load_n),
         genetic_load_lower_ci = genetic_load_mean - abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se),
         genetic_load_upper_ci = genetic_load_mean + abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se))

p3 <- ggplot(data = model_output_gen150) +
  facet_grid(. ~ DRIVE_TYPE) +
  geom_raster(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, fill = genetic_load_mean)) +
  scale_fill_gradient(high = met.brewer("Hiroshige",3)[3], name = "Mean genetic load") +
  scale_colour_gradient(high = met.brewer("Hiroshige",3)[3], name = "95% confidence interval genetic load") +
  xlab("Total cut rate") +
  ylab("Drive conversion rate") +
  PaperTheme + guides(color = "none")
p3

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

model_output_gen150 <- model_output_labels[model_output_labels$generation %in% 150:159,] %>%
  group_by(DRIVE_TYPE, TOTAL_CUT_RATE, DRIVE_CONVERSION_RATE) %>%
  mutate_all(~replace(., is.nan(.), NA)) %>%
  mutate(genetic_load = replace_na(genetic_load, 1)) %>%
  summarise(genetic_load_mean = mean(genetic_load, na.rm = TRUE),
            genetic_load_sd = sd(genetic_load, na.rm = TRUE),
            genetic_load_n = sum(!is.na(genetic_load))) %>%
  mutate(genetic_load_se = genetic_load_sd / sqrt(genetic_load_n),
         genetic_load_lower_ci = genetic_load_mean - abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se),
         genetic_load_upper_ci = genetic_load_mean + abs(qt(1 - (0.05 / 2), genetic_load_n - 1) * genetic_load_se))

p8 <- ggplot(data = model_output_gen150) +
  facet_grid(DRIVE_TYPE ~ .) +
  geom_raster(aes(x = TOTAL_CUT_RATE, y = DRIVE_CONVERSION_RATE, fill = genetic_load_mean)) +
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

save(p3, p8, file = "Fig1_cut_conversion_rates_GL_plot.Rdata")
