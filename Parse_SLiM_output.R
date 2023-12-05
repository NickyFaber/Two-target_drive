parse_slim_output <- function (raw_slim_output_matrix, generations, inputs) {
  # Elements in raw_slim_output_matrix contain all of the text that SLiM outputs.
  # It could be parsed and put into a matrix, or just flattened into a list, as below:
  parsed_output = c()
  for (col in raw_slim_output_matrix)
    for (row in col)
      for (line in row)
        if (startsWith(line, "OUT:"))
          # Strip off the "OUT:" from the beginning of the line.
          parsed_output = c(parsed_output, substring(line, 6))
  
  model_output <- as_tibble(do.call(rbind, strsplit(parsed_output, split = ",")), .name_repair = "minimal")
  model_output <- cbind(sapply(model_output[,1],as.numeric),
                        model_output[,2:3],
                        sapply(model_output[,4:25], as.numeric))
  colnames(model_output) <- c("generation", 
                              "DRIVE_TYPE", "GRNA_SATURATION_SIMULATED", 
                              "DRIVE_CONVERSION_RATE", "TOTAL_CUT_RATE", 
                              "EMBRYO_RESISTANCE_CUT_RATE", "DRIVE_FITNESS_VALUE", 
                              "R1_OCCURRENCE_RATE", "OFF_TARGET_CUT_RATE",
                              "OFF_TARGET_EMBRYO", "NUM_GRNAS", "INTRO_FREQ", 
                              "OFF_TARGET_FITNESS_MULTIPLIER", "SOMATIC_EXPRESSION_FITNESS_MULTIPLIER", 
                              "GLOBAL_SATURATION_FACTOR",
                              "popSize", "expectedPopSize", "bonus_pop_factor",
                              "freq_wt", "freq_dr", "freq_has_drive", "freq_r1", "freq_r2", 
                              "freq_offtarget_wt", "freq_offtarget_disrupted")

  model_output$repetition <- rep(1:nrow(inputs), each=generations)
  model_output$DRIVE_TYPE <- factor(model_output$DRIVE_TYPE, 
                                    levels = c("fertility","Xlinkedfertility","recessivelethal","fertility_fertility","neutral_fertility","haplolethal_fertility","haplolethal_recessivelethal","haplolethal_Xlinkedfertility","recessivelethal_fertility","recessivelethal_recessivelethal","recessivelethal_Xlinkedfertility"), 
                                    labels = c("Fertility","X-linked-fertility","Recessivelethal","Fertility Fertility","Neutral Fertility","Haplolethal Fertility","Haplolethal Recessivelethal","Haplolethal X-linked-fertility","Recessivelethal Fertility","Recessivelethal Recessivelethal","Recessivelethal X-linked-fertility"))

  return(model_output)
}



