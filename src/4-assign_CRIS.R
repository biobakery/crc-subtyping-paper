rm(list = ls())
library(CRISclassifier)
library(tidyverse)
source("src/misc.R")

load("data/eSets/setNames.RData")
for(set in setNames) load(paste0("data/eSets/new/",
                                 set, ".RData"))
data("features")

dir.create("results/CRIS/classify/tmp/",
           recursive = T,
           showWarnings = F)
seed <- 1
for(set in setNames) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type == "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf %>%
    apply(1, function(x) x - mean(x)) %>%
    t
  genes_common <- intersect(rownames(exprs_tmp), features$`Gene ID`)
  exprs_tmp <- exprs_tmp[genes_common, ]
  
  exprs_tmp %>%
    data.frame(symbol = rownames(exprs_tmp), .) %>%
    write.table(file = paste0("results/CRIS/classify/tmp/",
                              set, ".txt"),
                quote = F,
                sep = "\t",
                row.names = F)
  
  cris_classifier(paste0("results/CRIS/classify/tmp/",
                         set, ".txt"),
                  paste0("results/CRIS/classify/tmp/",
                         set, "_classify.txt"),
                  nresmpl = 1000,
                  rnd.seed = seed)
  CRIS_results <- read_tsv(paste0("results/CRIS/classify/tmp/",
                                  set, 
                                  "_classify.txt_prediction_result.xls")) %>% 
    data.frame(row.names = 1,
               check.names = F)
  colnames(CRIS_results) <- c("CRIS_label", paste0("CRIS_", colnames(CRIS_results)[-1]))
  CRIS_results <- CRIS_results %>% 
    rownames_to_column("samples") %>% 
    mutate(CRIS_label = CRIS_label %>% 
             recode(.missing = "unlabeled") %>% 
             factor(levels = c("CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E", "unlabeled"))) %>% 
    column_to_rownames("samples")
  if(!(all(rownames(pdata_tmp) == rownames(CRIS_results))))
    stop("Row names don't match!")
  pdata_tmp <- merge(pdata_tmp, CRIS_results, by = "row.names") %>% 
    data.frame(row.names = 1, check.names = F)
  pdata <- pdata %>% left_join(pdata_tmp)
  rownames(pdata) <- sampleNames(eSet)
  pData(eSet) <- pdata
  assign(set, eSet)
  save(list = set, file = paste0(
    "data/eSets/new/",
    set,
    ".RData"
  ))
}
