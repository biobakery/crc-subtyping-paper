rm(list=ls())
library(tidyverse)
source("src/misc.R")

load('data/eSets/setNames.RData')
for (set in setNames) {
  load(set %>% paste0('data/eSets/new/', ., '.RData'))
}

# fill in survival info from Isella
data_Isella <- read_tsv("data/CRIS/Isella_Supp_Tab_16.txt", col_names = F) %>%
  mutate(X2 = gsub("-", ".", X2, fixed = T)) %>%
  as.data.frame
rownames(data_Isella) <- data_Isella$X2
unique(data_Isella$X1)

for(set in c("GSE17536_eset", "TCGA.COAD_eset", "TCGA.RNASeqV2_eset")) {
  eset <- get(set)
  pdata <- pData(eset)
  common_samples <- intersect(sampleNames(eset), data_Isella$X2)
  pdata[common_samples, ]$days_to_recurrence_or_death <-
    data_Isella[common_samples, ]$X3 * 365
  pdata[common_samples, ]$dfs_status <-
    c("living_norecurrence", "deceased_or_recurrence")[
      data_Isella[common_samples, ]$X4 + 1
      ]
  pData(eset) <- pdata
  assign(x = set, value = eset)
  save(list = set, file = paste0("data/eSets/new/", set, ".RData"))
}

# fill in survival from CRCSC
df_CMS_clinical <- read_tsv("data/crcsc/mergedPhenotype/clinical_molecular_public_all.txt") %>%
  mutate(sample = sample %>%
           gsub("-", ".", ., fixed = T),
         dataset = dataset %>%
           gsub("gse", "GSE", ., fixed = T) %>%
           paste0("_eset")) %>%
  filter(!is.na(rfsStat))
df_CMS_clinical[df_CMS_clinical$dataset == "GSE33113_eset", ]$sample <-
  sampleNames(GSE33113_eset)[match(pData(GSE33113_eset)$alt_sample_name,
                                   df_CMS_clinical[df_CMS_clinical$dataset == "GSE33113_eset", ]$sample)] %>%
  na.omit
df_CMS_clinical <- as.data.frame(df_CMS_clinical)
rownames(df_CMS_clinical) <- df_CMS_clinical$sample
for(set in c("GSE33113_eset")) {
  eset <- get(set)
  pdata <- pData(eset)
  common_samples <- intersect(sampleNames(eset), df_CMS_clinical$sample)
  pdata[common_samples, ]$dfs_status <-
    c("living_norecurrence", "deceased_or_recurrence")[
      df_CMS_clinical[common_samples, ]$rfsStat + 1
      ]
  pData(eset) <- pdata
  assign(x = set, value = eset)
  save(list = set, file = paste0("data/eSets/new/", set, ".RData"))
}
