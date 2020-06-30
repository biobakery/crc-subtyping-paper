rm(list = ls())
source("src/misc.R")
library(metafor)
library(logistf)
library(tidyverse)

load("data/eSets/setNames.RData")
for (set in setNames) {
  load(paste0("data/eSets/new/", set, ".RData"))
}

# CRIS subtypes with GSE76402 and GSE59857 ---------------------------------------------
load("data/CRIS/GSE76402_eset.RData")
load("data/CRIS/GSE59857_eset.RData")
eSets <- list()
for (set in c(setNames, "GSE76402_eset", "GSE59857_eset")) {
  eset.tmp <- get(set)
  pdata.tmp <- pData(eset.tmp)
  pdata.tmp <- pdata.tmp %>% mutate(CRISA = (CRIS_label == "CRIS-A"),
                                    CRISB = (CRIS_label == "CRIS-B"),
                                    CRISC = (CRIS_label == "CRIS-C"),
                                    CRISD = (CRIS_label == "CRIS-D"),
                                    CRISE = (CRIS_label == "CRIS-E"))
  pData(eset.tmp) <- pdata.tmp
  eSets[[set]] <- eset.tmp
}

output.dir <- "results/clinical/"
dir.create(output.dir, showWarnings = F, recursive = T)

vars <- c("CRISB")
vars.contrasts <- c("TRUE")
names(vars.contrasts) <- vars

var <- vars[1]


effect.size <- matrix(NA, 4, 0)
rownames(effect.size) <- c("PCSS1.uni", "PCSS2.uni", "PCSS1.biv", "PCSS2.biv")
effect.se <- effect.size

for(set in c(setNames, "GSE76402_eset", "GSE59857_eset")) {
  p.values <- matrix(NA, 5, 3)
  rownames(p.values) <- c("n", "PCSS1.uni", "PCSS2.uni", "PCSS1.biv", "PCSS2.biv")
  colnames(p.values) <- c("effect_size", "sd", "p_value")
  
  pData.tmp <- pData(eSets[[set]])[pData(eSets[[set]])[, "sample_type"]=="tumor", ]
  outcome <- (pData.tmp[, var] == vars.contrasts[var])
  
  count.table <- table(outcome)
  if(length(count.table) <= 1) next # all outcome values are the same, not able to evaluate
  else {
    p.values["n", 1] <- sum(!is.na(pData.tmp[, var]))
    p.values["PCSS1.uni", ] <- tryCatch(summary(glm((outcome == "TRUE") ~ PCSS1, 
                                                    data=pData.tmp, family="binomial"))$coef[2, c(1, 2, 4)], 
                                        warning = function(w) {
                                          logistf.fit.tmp <- logistf((outcome == "TRUE") ~ PCSS1, data=pData.tmp)
                                          c(logistf.fit.tmp$coef[2], 
                                             sqrt(logistf.fit.tmp$var[2, 2]), 
                                             logistf.fit.tmp$prob[2]) %>% return
                                        })
    p.values["PCSS2.uni", ] <- tryCatch(summary(glm((outcome == "TRUE") ~ PCSS2, 
                                                    data=pData.tmp, family="binomial"))$coef[2, c(1, 2, 4)], 
                                        warning = function(w) {
                                          logistf.fit.tmp <- logistf((outcome == "TRUE") ~ PCSS2, data=pData.tmp)
                                          c(logistf.fit.tmp$coef[2], 
                                             sqrt(logistf.fit.tmp$var[2, 2]), 
                                             logistf.fit.tmp$prob[2]) %>% return
                                        })
    p.values[c("PCSS1.biv", "PCSS2.biv"), ] <- tryCatch(summary(glm((outcome == "TRUE") ~ PCSS1 + PCSS2, 
                                                                    data=pData.tmp, family="binomial"))$coef[2:3, c(1, 2, 4)], 
                                                        warning = function(w) {
                                                          logistf.fit.tmp <- logistf((outcome == "TRUE") ~ PCSS1 + PCSS2, data=pData.tmp)
                                                          c(logistf.fit.tmp$coef[2:3], 
                                                             c(sqrt(logistf.fit.tmp$var[2, 2]), 
                                                                sqrt(logistf.fit.tmp$var[3, 3])), 
                                                             logistf.fit.tmp$prob[2:3])
                                                        })
  }
  
  effect.size <- cbind(effect.size, p.values[2:5, 1])
  colnames(effect.size)[ncol(effect.size)] <- set
  effect.se <- cbind(effect.se, p.values[2:5, 2])
  colnames(effect.se)[ncol(effect.se)] <- set
}

# meta-analysis
dimnames(results.rma) <- list(c("PCSS1.uni", "PCSS2.uni", "PCSS1.biv", "PCSS2.biv"), 
                                c("Effect", "se", "pValue", "I2"))
pdf("manuscript/suppFig_CRISB.pdf", width = 12, height = 6)
par(mfrow=c(1, 2))
rma.fit <- rma(yi=effect.size[1, ], sei=effect.se[1, ], method="FE")
forest(rma.fit, 
       slab=colnames(effect.size) %>% gsub("_eset", "", ., fixed = T), 
       xlab="CRIS-B log OR",
       main = "PCSS1")
rma.fit <- rma(yi=effect.size[2, ], sei=effect.se[2, ], method="FE")
forest(rma.fit, 
       slab=colnames(effect.size) %>% gsub("_eset", "", ., fixed = T), 
       xlab="CRIS-B log OR",
       main = "PCSS2")
dev.off()
