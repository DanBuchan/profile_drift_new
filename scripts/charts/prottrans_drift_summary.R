library(ggplot2)
library(reshape2)
library(tidyr)
library(tidyverse)

prottrans_drift_summary <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/best_hits/prottrans_closest//drift_summary.csv", header=T)
target_lookup <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/generation_or_af_targets/alphafold_targets.csv")
colnames(target_lookup) <- c('target_family', 'type')

prottrans_merged_data <- merge(target_lookup, prottrans_drift_summary, by="target_family")

masked25_data <- subset(prottrans_merged_data, mask == "masked25")
masked50_data <- subset(prottrans_merged_data, mask == "masked50")
masked75_data <- subset(prottrans_merged_data, mask == "masked75")

masked25_data["percentage_correct"] <- masked25_data["correct_count"]/masked25_data["tot_seqs"]
masked25_data["percentage_drift"] <- masked25_data["drift_count"]/masked25_data["tot_seqs"]
masked25_data["percentage_nobest"] <- masked25_data["nobest"]/masked25_data["tot_seqs"]

percentages25 <- aggregate(masked25_data$percentage_correct, list(masked25_data$type), mean)
colnames(percentages25) <- c('type', 'mean_correct')
percentages25$mean_drift <- tapply(masked25_data$percentage_drift, list(masked25_data$type), mean)
percentages25$mean_nobest <-tapply(masked25_data$percentage_nobest, list(masked25_data$type), mean)

###

masked50_data["percentage_correct"] <- masked50_data["correct_count"]/masked50_data["tot_seqs"]
masked50_data["percentage_drift"] <- masked50_data["drift_count"]/masked50_data["tot_seqs"]
masked50_data["percentage_nobest"] <- masked50_data["nobest"]/masked50_data["tot_seqs"]

percentages50 <- aggregate(masked50_data$percentage_correct, list(masked50_data$type), mean)
colnames(percentages50) <- c('type', 'mean_correct')
percentages50$mean_drift <- tapply(masked50_data$percentage_drift, list(masked50_data$type), mean)
percentages50$mean_nobest <-tapply(masked50_data$percentage_nobest, list(masked50_data$type), mean)

###

masked75_data["percentage_correct"] <- masked75_data["correct_count"]/masked75_data["tot_seqs"]
masked75_data["percentage_drift"] <- masked75_data["drift_count"]/masked75_data["tot_seqs"]
masked75_data["percentage_nobest"] <- masked75_data["nobest"]/masked75_data["tot_seqs"]

percentages75 <- aggregate(masked75_data$percentage_correct, list(masked75_data$type), mean)
colnames(percentages75) <- c('type', 'mean_correct')
percentages75$mean_drift <- tapply(masked75_data$percentage_drift, list(masked75_data$type), mean)
percentages75$mean_nobest <-tapply(masked75_data$percentage_nobest, list(masked75_data$type), mean)

write.csv(percentages25, "/home/dbuchan/Projects/profile_drift/results_data/drift/best_hits/prottrans25_drift_percentages.csv")
write.csv(percentages50, "/home/dbuchan/Projects/profile_drift/results_data/drift/best_hits/prottrans50_drift_percentages.csv")
write.csv(percentages75, "/home/dbuchan/Projects/profile_drift/results_data/drift/best_hits/prottrans75_drift_percentages.csv")




