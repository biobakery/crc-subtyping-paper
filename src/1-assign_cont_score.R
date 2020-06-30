rm(list = ls())
library(tidyverse)
source("src/misc.R")

load("data/eSets/setNames.RData")
for(set in setNames) {
  load(paste0("data/eSets/", set, ".RData"))
}
avg_loadings <- read.csv("data/loadings/avg_loadings.csv",
                         row.names = 1) %>% 
  as.matrix

for (set in setNames ) {
  print( set )
  eSet <- get(set)
  exprs <- exprs(eSet) %>% rmNaInf
  exprs <- apply(exprs, 1, function(x) x - mean(x)) %>% t
  genesCommon <- intersect(rownames(exprs), rownames(avg_loadings))
  pcss <- t(exprs[genesCommon, ]) %*% avg_loadings[genesCommon, ]
  pcss <- (t(pcss) / apply(pcss, 2, sd)) %>% t
  colnames(pcss) <- paste0( 'PCSS', 1:4 ) 
  
  pdata <- pData(eSet)
  pdata[, paste0("PCSS", 1:4)] <- pcss
  pData(eSet) <- pdata
  
  assign(x = set, value = eSet)
  save(list = set, file = paste0("data/eSets/new/",
                                 set, ".RData"))
}
