rm(list=ls())
library(Biobase)
library(tidyverse)
library(cluster)
library(fpc)
library(NMF)
library(ConsensusClusterPlus)
source("src/misc.R" )

load("data/eSets/setNames.RData")
for (set in setNames) {
  load(set %>% paste0("data/eSets/new/", ., ".RData"))
}
genes_common <- lapply(setNames, function(set) set %>% get %>% featureNames) %>% 
  Reduce("intersect", .)
dir.create("results/clustering/", recursive = T, showWarnings = F)
var_perc <- sapply(setNames, function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)[genes_common, ]
  pdata_tmp <- pData(eSet_tmp)
  
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf %>% 
    apply(1, function(x) x - mean(x)) %>% 
    t
  var_all <- apply(exprs_tmp, 1, var, na.rm = T)
  sum(var_all[order(var_all, decreasing = T)[1:3000]]) / sum(var_all)
})
save(var_perc, file = "results/clustering/var_perc.RData")

# clustering with most variable genes
K_max <- 8
for(set in setNames) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)[genes_common, ]
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf
  pdata_tmp <- pData(eSet_tmp)
  
  var_all <- apply(exprs_tmp, 1, var, na.rm = T)
  exprs_tmp <- exprs_tmp[order(var_all, decreasing = T)[1:3000], ]
  exprs_tmp_centered <- exprs_tmp %>% 
    apply(1, function(x) x - mean(x)) %>% 
    t
  
  # Gap statistic with pam clustering
  df_toreturnGap <- data.frame(clusGap(t(exprs_tmp_centered), 
                                       FUN = pam, 
                                       K.max = K_max, 
                                       B = 100)$Tab[-1, 3:4])
  colnames(df_toreturnGap) <- c("statistics", "se")
  df_toreturnGap <- df_toreturnGap %>% 
    mutate(K = 2:K_max,
           metric = "Gap Statistic",
           method = "k-medoid")
  
  
  # Prediction strength with pam clustering
  predstr_tmp <- prediction.strength(t(exprs_tmp_centered), 
                                     Gmin = 2, 
                                     Gmax = K_max, 
                                     M = 100, 
                                     clustermethod=claraCBI)
  df_toreturnPredStr <- data.frame(statistics = predstr_tmp$mean.pred[2:K_max], 
                                   se = sapply(predstr_tmp$predcorr[2:K_max], sd),
                                   K = 2:K_max, 
                                   metric = "Prediction Strength",
                                   method = "k-medoid")
  
  # Silhouette Width with different clustering structures
  dist_eucl <- dist(t(exprs_tmp_centered), method = "eucl")
  dist_pearson <- (1 - cor(exprs_tmp, method = "pearson")) %>% as.dist
  # pam
  l_clst_pam <- lapply(
    2:K_max, 
    function(k) pam(x = dist_eucl, k = k)$clustering
  )
  sw_pam <- sapply(
    l_clst_pam, 
    function(clst_tmp) silhouette(clst_tmp, dist = dist_eucl)[, "sil_width"]
  )
  df_toreturnSWpam <- data.frame(statistics = apply(sw_pam, 2, mean),
                                 se = apply(sw_pam, 2, sd),
                                 K = 2:K_max,
                                 metric = "Silhouette Width",
                                 method = "k-medoid")     
  # NMF
  if(mean(exprs_tmp < 0, na.rm = T) < 0.01) {
    exprs_tmp_nmf <- exprs_tmp
    exprs_tmp_nmf[exprs_tmp_nmf < 0] <- 0
    l_clst_nmf <- lapply(
      2, 
      function(k) nmf(exprs_tmp_nmf, rank = k, nrun = 30,
                      .options = list(p = 6, parallel.required = TRUE)) %>% 
        predict(what = "samples") %>% 
        as.numeric
    )
    sw_nmf <- sapply(
      l_clst_nmf, 
      function(clst_tmp) silhouette(clst_tmp, dist = dist_eucl)[, "sil_width"]
    )
    df_toreturnSWnmf <- data.frame(statistics = apply(sw_nmf, 2, mean),
                                   se = apply(sw_nmf, 2, sd),
                                   K = 2:K_max,
                                   metric = "Silhouette Width",
                                   method = "NMF") 
  } else {
    df_toreturnSWnmf <- data.frame(statistics = rep(NA, K_max - 1),
                                   se = rep(NA, K_max - 1),
                                   K = 2:K_max,
                                   metric = "Silhouette Width",
                                   method = "NMF") 
  }
  # consensus hc
  conshc_fit <- ConsensusClusterPlus(d = dist_pearson, 
                                     maxK = K_max, reps = 100, pItem = 0.8, pFeature = 1,
                                     clusterAlg = "hc", innerLinkage = "ward.D")
  l_clst_conshc <- lapply(
    2:K_max,
    function(k) conshc_fit[[k]]$consensusClass
  )
  sw_conshc <- sapply(
    l_clst_conshc, 
    function(clst_tmp) silhouette(clst_tmp, dist = dist_pearson)[, "sil_width"]
  )
  df_toreturnSWconshc <- data.frame(statistics = apply(sw_conshc, 2, mean),
                                 se = apply(sw_conshc, 2, sd),
                                 K = 2:K_max,
                                 metric = "Silhouette Width",
                                 method = "Consensus Hierarchical")
  assign(x = paste0("df_", set), 
         value = rbind(df_toreturnGap, 
                       df_toreturnPredStr, 
                       df_toreturnSWpam,
                       df_toreturnSWnmf, 
                       df_toreturnSWconshc)
         )
  save(list = paste0("df_", set),
       file = paste0("results/clustering/df_", set, ".RData"))
}

# clustering with CRIS genes
library(CRISclassifier)
data("features")
dir.create("results/clustering/CRIS/", recursive = T, showWarnings = F)
load("data/CRIS/GSE59857_eset.RData")
load("data/CRIS/GSE76402_eset.RData")
K_max <- 8
for(set in c(setNames, "GSE59857_eset", "GSE76402_eset")) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf
  
  exprs_tmp <- exprs_tmp[intersect(rownames(exprs_tmp), features$`Gene ID`), ]
  exprs_tmp_centered <- exprs_tmp %>%
    apply(1, function(x) x - mean(x)) %>%
    t

  # Gap statistic with pam clustering
  df_toreturnGap <- data.frame(clusGap(t(exprs_tmp_centered),
                                       FUN = pam,
                                       K.max = K_max,
                                       B = 100)$Tab[-1, 3:4])
  colnames(df_toreturnGap) <- c("statistics", "se")
  df_toreturnGap <- df_toreturnGap %>%
    mutate(K = 2:K_max,
           metric = "Gap Statistic",
           method = "k-medoid")


  # Prediction strength with pam clustering
  predstr_tmp <- prediction.strength(t(exprs_tmp_centered),
                                     Gmin = 2,
                                     Gmax = K_max,
                                     M = 100,
                                     clustermethod=claraCBI)
  df_toreturnPredStr <- data.frame(statistics = predstr_tmp$mean.pred[2:K_max],
                                   se = sapply(predstr_tmp$predcorr[2:K_max], sd),
                                   K = 2:K_max,
                                   metric = "Prediction Strength",
                                   method = "k-medoid")

  # Silhouette Width with different clustering structures
  dist_eucl <- dist(t(exprs_tmp_centered), method = "eucl")
  dist_pearson <- (1 - cor(exprs_tmp, method = "pearson")) %>% as.dist
  # pam
  l_clst_pam <- lapply(
    2:K_max,
    function(k) pam(x = dist_eucl, k = k)$clustering
  )
  sw_pam <- sapply(
    l_clst_pam,
    function(clst_tmp) silhouette(clst_tmp, dist = dist_eucl)[, "sil_width"]
  )
  df_toreturnSWpam <- data.frame(statistics = apply(sw_pam, 2, mean),
                                 se = apply(sw_pam, 2, sd),
                                 K = 2:K_max,
                                 metric = "Silhouette Width",
                                 method = "k-medoid")
  # NMF
  if(mean(exprs_tmp < 0, na.rm = T) < 0.01) {
    exprs_tmp_nmf <- exprs_tmp
    exprs_tmp_nmf[exprs_tmp_nmf < 0] <- 0
    l_clst_nmf <- lapply(
      2,
      function(k) nmf(exprs_tmp_nmf, rank = k, nrun = 30,
                      .options = list(p = 6, parallel.required = TRUE)) %>%
        predict(what = "samples") %>%
        as.numeric
    )
    sw_nmf <- sapply(
      l_clst_nmf,
      function(clst_tmp) silhouette(clst_tmp, dist = dist_eucl)[, "sil_width"]
    )
    df_toreturnSWnmf <- data.frame(statistics = apply(sw_nmf, 2, mean),
                                   se = apply(sw_nmf, 2, sd),
                                   K = 2:K_max,
                                   metric = "Silhouette Width",
                                   method = "NMF")
  } else {
    df_toreturnSWnmf <- data.frame(statistics = rep(NA, K_max - 1),
                                   se = rep(NA, K_max - 1),
                                   K = 2:K_max,
                                   metric = "Silhouette Width",
                                   method = "NMF")
  }
  # consensus hc
  conshc_fit <- ConsensusClusterPlus(d = dist_pearson,
                                     maxK = K_max, reps = 100, pItem = 0.8, pFeature = 1,
                                     clusterAlg = "hc", innerLinkage = "ward.D")
  l_clst_conshc <- lapply(
    2:K_max,
    function(k) conshc_fit[[k]]$consensusClass
  )
  sw_conshc <- sapply(
    l_clst_conshc,
    function(clst_tmp) silhouette(clst_tmp, dist = dist_pearson)[, "sil_width"]
  )
  df_toreturnSWconshc <- data.frame(statistics = apply(sw_conshc, 2, mean),
                                    se = apply(sw_conshc, 2, sd),
                                    K = 2:K_max,
                                    metric = "Silhouette Width",
                                    method = "Consensus Hierarchical")
  assign(x = paste0("df_", set),
         value = rbind(df_toreturnGap,
                       df_toreturnPredStr,
                       df_toreturnSWpam,
                       df_toreturnSWnmf,
                       df_toreturnSWconshc)
  )
  save(list = paste0("df_", set),
       file = paste0("results/clustering/CRIS/df_", set, ".RData"))
}
