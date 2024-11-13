library(ggplot2)
library(reshape2)
library(tidyr)
library(tidyverse)

hhm_drift_summary <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/best_hits/hmm_closest/drift_summary.csv", header=T)
target_lookup <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/generation_or_af_targets/alphafold_targets.csv", header=T)
colnames(target_lookup) <- c('target_family', 'type')

merged_data <- merge(target_lookup,hhm_drift_summary, by="target_family")

# classes <- unique(merged_data["type"])
merged_data["percentage_correct"] <- merged_data["correct_count"]/merged_data["tot_seqs"]
merged_data["percentage_drift"] <- merged_data["drift_count"]/merged_data["tot_seqs"]
merged_data["percentage_nobest"] <- merged_data["nobest"]/merged_data["tot_seqs"]

percentages_hmm <- aggregate(merged_data$percentage_correct, list(merged_data$type), mean)
colnames(percentages_hmm) <- c('type', 'mean_correct')
percentages_hmm$mean_drift <- tapply(merged_data$percentage_drift, list(merged_data$type), mean)
percentages_hmm$mean_nobest <-tapply(merged_data$percentage_nobest, list(merged_data$type), mean)

write.csv(percentages_hmm, "/home/dbuchan/Projects/profile_drift/results_data/drift/best_hits/hmm_drift_percentages.csv")
