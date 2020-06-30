rm(list=ls())
library(tidyverse)
library(CMSclassifier)
library(org.Hs.eg.db)
source("src/misc.R")

load('data/eSets/setNames.RData')
for (set in setNames) {
  load(paste0('data/eSets/new/', 
              set, '.RData'))
}

# overall CMS labels
df_CMS_orig <- read_tsv("data/crcsc/mergedPhenotype/cms_labels_public_all.txt") %>% 
  mutate(sample = sample %>% 
           gsub("-", ".", ., fixed = T),
         dataset = dataset %>% 
           gsub("gse", "GSE", ., fixed = T) %>% 
           paste0("_eset"))
df_CMS_orig[df_CMS_orig$dataset == "GSE33113_eset", ]$sample <- 
  sampleNames(GSE33113_eset)[match(pData(GSE33113_eset)$alt_sample_name,
                                   df_CMS_orig[df_CMS_orig$dataset == "GSE33113_eset", ]$sample)] %>% 
  na.omit

data("model")
df_mapping_RF <- select(org.Hs.eg.db, 
                        finalModel$importance %>% rownames, 
                        columns = "SYMBOL")
data("centroids")
df_mapping_SSP <- select(org.Hs.eg.db, 
                         rownames(centroids), 
                         columns = "SYMBOL")
for (set in setNames) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type == "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf %>% 
    apply(1, function(x) x - mean(x)) %>% 
    t
  
  # classify RF
  exprs_tmp_RF <- exprs_tmp[intersect(rownames(exprs_tmp), df_mapping_RF$SYMBOL),]
  rownames(exprs_tmp_RF) <- intersect(rownames(exprs_tmp), df_mapping_RF$SYMBOL) %>% 
    match(df_mapping_RF$SYMBOL) %>% 
    `[`(df_mapping_RF$ENTREZID, .)
  exprs_tmp_RF <- exprs_tmp_RF %>% fillInZero(df_mapping_RF$ENTREZID)
  lbl_RF <- classifyCMS.RF(exprs_tmp_RF)$RF.predictedCMS %>% 
    replace(., is.na(.), "unlabeled") %>% 
    factor(levels = c("CMS1", "CMS2", "CMS3", "CMS4", "unlabeled"))
  
  # classify SSP
  exprs_tmp_SSP <- exprs_tmp[intersect(rownames(exprs_tmp), df_mapping_SSP$SYMBOL),]
  rownames(exprs_tmp_SSP) <- intersect(rownames(exprs_tmp), df_mapping_SSP$SYMBOL) %>% 
    match(df_mapping_SSP$SYMBOL) %>% 
    `[`(df_mapping_SSP$ENTREZID, .)
  exprs_tmp_SSP <- exprs_tmp_SSP %>% fillInZero(df_mapping_SSP$ENTREZID)
  lbl_SSP <- classifyCMS.SSP(exprs_tmp_SSP)$SSP.predictedCMS %>% 
    replace(., is.na(.), "unlabeled") %>% 
    factor(levels = c("CMS1", "CMS2", "CMS3", "CMS4", "unlabeled"))
  
  # classification from CRCSC
  lbl_CRCSC <- sapply(pdata_tmp %>% rownames, function(sample) {
    ifelse(sample %in% df_CMS_orig$sample, 
           df_CMS_orig$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples[sample == df_CMS_orig$sample], 
           NA)
  }) %>% 
    replace(. == "NOLBL", "unlabeled") %>% 
    factor(levels = c("CMS1", "CMS2", "CMS3", "CMS4", "unlabeled"))
  
  pdata_tmp$cms_label_orig <- lbl_CRCSC
  pdata_tmp$cms_label_RF <- lbl_RF
  pdata_tmp$cms_label_SSP <- lbl_SSP
  pdata <- pdata %>% left_join(pdata_tmp)
  rownames(pdata) <- sampleNames(eSet)
  pData(eSet) <- pdata
  assign(set, eSet)
  save(list = set,
       file = paste0("data/eSets/new/",
                     set,
                     ".RData"))
}
