rm(list = ls())
source("src/misc.R")
library(tidyverse)

load("data/eSets/trainingSetNames.RData")
load("data/eSets/setNames.RData")
summaryTable <- read.table("data/eSets/datasetSummary.txt", 
                           row.names=1, header=T, sep="\t")

table1 <- summaryTable[, 1:7]
table1$latestage.perc <- format(round(table1$latestage.perc, 2), nsmall=2)
table1[ is.na(table1) ] <- "-"
table1[ table1 == "    NA"  ] <- "-"
table1[ table1 == "   NA" ] <- "-"

table1 <- as.matrix(table1)
table1 <- rbind(c("Training Sets", rep("", ncol(table1) -1)), 
                 table1[trainingSetNames, ], 
                 c("Validation Sets", rep("", ncol(table1) -1)), 
                 table1[setdiff(setNames, trainingSetNames), ])

colnames(table1) <- c("Dataset", 
                      "Accession ID", 
                      "Platform", 
                      "# Normal / Tumor Samples", 
                      "% Late Stage Tumors", 
                      "Staging System", 
                      "Availability of Metastasis Info")
write.table(table1, "manuscript/table1.txt", quote=F, row.names=F, sep="\t")
