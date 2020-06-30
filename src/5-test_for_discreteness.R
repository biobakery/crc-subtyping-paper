rm(list = ls())
source("src/misc.R")
library(tidyverse)
library(CMSclassifier)
library(CRISclassifier)
library(org.Hs.eg.db)
library(cluster)

load('data/eSets/setNames.RData')
for (set in setNames) {
  load(set %>% paste0('data/eSets/new/', ., '.RData'))
}

data("centroids")
df_mapping_SSP <- select(org.Hs.eg.db, 
                         rownames(centroids), 
                         columns = "SYMBOL")
data("features")

# silhouette width between pre-defined subtypes ---------------------------
df_silWidth <- lapply(setNames, function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf %>% 
    apply(1, function(x) x - mean(x)) %>% 
    t
  
  exprs_tmp_SSP <- exprs_tmp[intersect(rownames(exprs_tmp), df_mapping_SSP$SYMBOL),]
  exprs_tmp_CRIS <- exprs_tmp[intersect(rownames(exprs_tmp), features$`Gene ID`), ]
  
  distances_SSP <- list(
    "Euclidean" = dist(exprs_tmp_SSP %>% t, 'eucl') %>% as.matrix,
    "Manhattan" = dist(exprs_tmp_SSP %>% t, "manh") %>% as.matrix,
    "1 - Pearson" = 1 - cor(exprs_tmp_SSP, method = "pearson"),
    "1 - Spearman" = 1 - cor(exprs_tmp_SSP, method = "spearman")
  )
  distances_CRIS <- list(
    "Euclidean" = dist(exprs_tmp_CRIS %>% t, 'eucl') %>% as.matrix,
    "Manhattan" = dist(exprs_tmp_CRIS %>% t, "manh") %>% as.matrix,
    "1 - Pearson" = 1 - cor(exprs_tmp_CRIS, method = "pearson"),
    "1 - Spearman" = 1 - cor(exprs_tmp_CRIS, method = "spearman")
  )
  ind_lbl_SSP <- pdata_tmp$cms_label_SSP != "unlabeled"
  ind_lbl_CRIS <- pdata_tmp$CRIS_label != "unlabeled"
  
  lapply(c("Euclidean",
           "Manhattan",
           "1 - Pearson",
           "1 - Spearman"), 
         function(metric) {
    dist_tmp_SSP <- distances_SSP[[metric]]
    dist_tmp_CRIS <- distances_CRIS[[metric]]
    rbind(
      data.frame(sw = silhouette(pdata_tmp$cms_label_SSP[ind_lbl_SSP] %>% 
                                   factor %>% 
                                   as.numeric,
                                 dist = dist_tmp_SSP[ind_lbl_SSP, ind_lbl_SSP] %>% 
                                   as.dist)[, "sil_width"],
                 samples = "labeled",
                 subtype = "CMS"),
      data.frame(sw = silhouette(pdata_tmp$cms_label_SSP %>% 
                                   factor %>% 
                                   as.numeric,
                                 dist = dist_tmp_SSP %>% 
                                   as.dist)[, "sil_width"],
                 samples = "all",
                 subtype = "CMS"),
      data.frame(sw = silhouette(pdata_tmp$CRIS_label[ind_lbl_CRIS] %>% 
                                   factor %>% 
                                   as.numeric,
                                 dist = dist_tmp_CRIS[ind_lbl_CRIS, ind_lbl_CRIS] %>% 
                                   as.dist)[, "sil_width"],
                 samples = "labeled",
                 subtype = "CRIS"),
      data.frame(sw = silhouette(pdata_tmp$CRIS_label %>% 
                                   factor %>% 
                                   as.numeric,
                                 dist = dist_tmp_CRIS %>% 
                                   as.dist)[, "sil_width"],
                 samples = "all",
                 subtype = "CRIS")
    ) %>% 
      data.frame(metric = metric)
  }) %>% 
    Reduce("rbind", .) %>% 
    data.frame(study = set) %>% 
    return
}) %>% Reduce("rbind", .)
dir.create("results/discrete_evaluation/", recursive = T, showWarnings = F)
save(df_silWidth, file = "results/discrete_evaluation/all_silWidth.RData")



# additional evaluation for CRIS datasets ---------------------------------
dir.create("results/discrete_evaluation/additional/", 
           recursive = T, showWarnings = F)
load("data/CRIS/GSE59857_eset.RData")
load("data/CRIS/GSE76402_eset.RData")
df_silWidth_additional <- lapply(c("GSE59857_eset", "GSE76402_eset"), function(set) {
  eSet_tmp <- get(set)
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf %>% 
    apply(1, function(x) x - mean(x)) %>% 
    t
  
  exprs_tmp_CRIS <- exprs_tmp[intersect(rownames(exprs_tmp), features$`Gene ID`), ]
  
  distances_CRIS <- list(
    "Euclidean" = dist(exprs_tmp_CRIS %>% t, 'eucl') %>% as.matrix,
    "Manhattan" = dist(exprs_tmp_CRIS %>% t, "manh") %>% as.matrix,
    "1 - Pearson" = 1 - cor(exprs_tmp_CRIS, method = "pearson"),
    "1 - Spearman" = 1 - cor(exprs_tmp_CRIS, method = "spearman")
  )
  ind_lbl_CRIS <- pdata_tmp$CRIS_label != "unlabeled"
  
  lapply(c("Euclidean",
           "Manhattan",
           "1 - Pearson",
           "1 - Spearman"), 
         function(metric) {
           dist_tmp_CRIS <- distances_CRIS[[metric]]
           rbind(
             data.frame(sw = silhouette(pdata_tmp$CRIS_label[ind_lbl_CRIS] %>% 
                                          factor %>% 
                                          as.numeric,
                                        dist = dist_tmp_CRIS[ind_lbl_CRIS, ind_lbl_CRIS] %>% 
                                          as.dist)[, "sil_width"],
                        samples = "labeled",
                        subtype = "CRIS"),
             data.frame(sw = silhouette(pdata_tmp$CRIS_label %>% 
                                          factor %>% 
                                          as.numeric,
                                        dist = dist_tmp_CRIS %>% 
                                          as.dist)[, "sil_width"],
                        samples = "all",
                        subtype = "CRIS")
           ) %>% 
             data.frame(metric = metric)
         }) %>% 
    Reduce("rbind", .) %>% 
    data.frame(study = set) %>% 
    return
}) %>% Reduce("rbind", .)
save(df_silWidth_additional, file = "results/discrete_evaluation/silWidth_additional.RData")
