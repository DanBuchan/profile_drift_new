library(ggplot2)
library(reshape2)
library(tidyr)
library(tidyverse)

drops <- c("extra")

# PF12878
contaminants_grew <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/17447_A0A1Y3DLV4.1_1002-1171_blast_summary.csv")
colnames(contaminants_grew) <- c("iteration","query_family","hit_family","count","extra")
contaminants_grew <- contaminants_grew[ , !(names(contaminants_grew) %in% drops)]
ggplot(contaminants_grew, aes(x=iteration, y=count)) + geom_line(aes(linetype = hit_family)) + xlab("Iteration") + ylab("Family membership Count") + scale_linetype_manual("Pfam Family", values = c(1, 2))

# PF00050
contaminants_complex <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/11711_P00995.2_32-79_blast_summary.csv")
colnames(contaminants_complex) <- c("iteration","query_family","hit_family","count","extra")
contaminants_complex <- contaminants_complex[ , !(names(contaminants_complex) %in% drops)]
ggplot(contaminants_complex, aes(x=iteration, y=count)) + geom_line(aes(linetype = hit_family)) + xlab("Iteration") + ylab("Family membership Count") + scale_linetype_manual("Pfam Family", values = c(1,2,3,4,5,6,7,8))

# PF17135
contaminants_purified <-read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/16478_E1RIT8.1_5-120_blast_summary.csv")
colnames(contaminants_purified) <- c("iteration","query_family","hit_family","count","extra")
contaminants_purified <- contaminants_purified[ , !(names(contaminants_purified) %in% drops)]
ggplot(contaminants_purified, aes(x=iteration, y=count)) + geom_line(aes(linetype = hit_family)) + xlab("Iteration") + ylab("Family membership Count") + scale_linetype_manual("Pfam Family", values = c(1, 2, 3))


# PF08711,insig_drift
insig_drift <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/9637_A0A5N6Y7X4.1_1850-1947_blast_summary.csv")
colnames(insig_drift) <- c("iteration","query_family","hit_family","count","extra")
insig_drift <- insig_drift[ , !(names(insig_drift) %in% drops)]
ggplot(insig_drift, aes(x=iteration, y=count)) + geom_line(aes(linetype = hit_family)) + xlab("Iteration") + ylab("Family membership Count") + scale_linetype_manual("Pfam Family", values = c(1, 2))

# PF10269
non_drift <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/18819_A0A7K6MAT4.1_14-124_blast_summary.csv")
colnames(non_drift) <- c("iteration","query_family","hit_family","count","extra")
non_drift <- non_drift[ , !(names(non_drift) %in% drops)]
ggplot(non_drift, aes(x=iteration, y=count)) + geom_line(aes(linetype = hit_family)) + xlab("Iteration") + ylab("Family membership Count") + scale_linetype_manual("Pfam Family", values = c(1))

# PF02063
query_purified <- read.csv("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/12444_A0A553NMI6.1_4-203_blast_summary.csv")
colnames(query_purified) <- c("iteration","query_family","hit_family","count","extra")
query_purified <- query_purified[ , !(names(query_purified) %in% drops)]
ggplot(query_purified, aes(x=iteration, y=count)) + geom_line(aes(linetype = hit_family)) + xlab("Iteration") + ylab("Family membership Count") + scale_linetype_manual("Pfam Family", values = c(1,2,3,4,5,6,7,8,9,10,11))
