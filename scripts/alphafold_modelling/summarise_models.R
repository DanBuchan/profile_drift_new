library(ggplot2)
library(reshape2)
library(tidyr)
library(tidyverse)

plDDT_summary <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/alphafold_models/plddt_summary.csv", header=T)

merizo_summary <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/alphafold_models/merizo_hits.csv", header=T)
