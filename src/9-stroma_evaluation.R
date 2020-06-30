rm(list = ls())
library(tidyverse)
library(cowplot)

load("data/eSets/new/TCGA.COAD_eset.RData")
load("data/eSets/new/TCGA.RNASeqV2_eset.RData")

TCGA_purity <- read_tsv("data/purity/TCGA_pancan12.sample_info.txt")
TCGA_purity <- TCGA_purity %>% 
  filter(disease %in% c("COAD", "READ")) %>% 
  mutate(sample_name = tcga_id %>% 
           gsub("-", ".", ., fixed = T) %>% 
           substr(1, 12))

COAD <- pData(TCGA.COAD_eset) %>% 
  mutate(sample_name = rownames(pData(TCGA.COAD_eset)),
         study = "COAD") %>% 
  left_join(TCGA_purity, by = "sample_name")
RNASeqV2 <- pData(TCGA.RNASeqV2_eset) %>% 
  mutate(sample_name = rownames(pData(TCGA.RNASeqV2_eset)),
         study = "RNASeqV2") %>% 
  left_join(TCGA_purity, by = "sample_name")

rbind(COAD, RNASeqV2) %>% 
  lm(abs_purity ~ PCSS1, .) %>% 
  summary

rbind(COAD, RNASeqV2) %>% 
  lm(abs_purity ~ PCSS2, .) %>% 
  summary
