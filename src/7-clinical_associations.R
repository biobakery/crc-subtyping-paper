rm(list=ls())
source("src/misc.R")
library(Biobase)
library(tidyverse)
library(survival)
library(metafor)
library(logistf)

load("data/eSets/setNames.RData")
for(set in setNames) load(paste0("data/eSets/new/",
                                 set, ".RData"))

eSets <- list()
for (set in setNames){
  eset.tmp <- get(set)
  pdata.tmp <- pData(eset.tmp)
  pdata.tmp <- pdata.tmp %>% mutate(CMS1 = (cms_label_SSP == 'CMS1'),
                                    CMS2 = (cms_label_SSP == 'CMS2'),
                                    CMS3 = (cms_label_SSP == 'CMS3'),
                                    CMS4 = (cms_label_SSP == 'CMS4'))
  pData(eset.tmp) <- pdata.tmp
  eSets[[set]] <- eset.tmp
}

output.dir <- 'results/clinical/'
dir.create(output.dir, showWarnings = F, recursive = T)

# binary outcomes ---------------------------------------------------------
vars <- c("msi", "summarylocation", "summarygrade", "summarystage", 
          "CMS1", "CMS2", "CMS3", "CMS4")
vars.contrasts <- c("MSI", "r", "high", "late", 
                    "TRUE", "TRUE", "TRUE", "TRUE")
names(vars.contrasts) <- vars

for (var in vars) {
  dir.create(paste0(output.dir, var, '/'), showWarnings=F)
  
  effect.size <- matrix(NA, 4, 0)
  rownames(effect.size) <- c('PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv')
  effect.se <- effect.size
  
  for(set in setNames) {
    p.values <- matrix(NA, 5, 3)
    rownames(p.values) <- c('n', 'PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv')
    colnames(p.values) <- c('effect_size', 'sd', 'p_value')
    
    pData.tmp <- pData(eSets[[set]])[pData( eSets[[set]] )[, 'sample_type']=='tumor', ]
    outcome <- (pData.tmp[, var] == vars.contrasts[var])
    
    count.table <- table(outcome)
    if(length(count.table) <= 1) next # all outcome values are the same, not able to evaluate
    else {
      p.values['n', 1] <- sum( !is.na( pData.tmp[, var] ) )
      p.values['PCSS1.uni', ] <- tryCatch(summary(glm((outcome == 'TRUE') ~ PCSS1, 
                                                      data=pData.tmp, family='binomial'))$coef[2, c(1, 2, 4)], 
                                          warning = function(w) {
                                            logistf.fit.tmp <- logistf((outcome == 'TRUE') ~ PCSS1, data=pData.tmp)
                                            c( logistf.fit.tmp$coef[2], 
                                               sqrt( logistf.fit.tmp$var[2, 2] ), 
                                               logistf.fit.tmp$prob[2] ) %>% return
                                          })
      p.values['PCSS2.uni', ] <- tryCatch(summary(glm((outcome == 'TRUE') ~ PCSS2, 
                                                      data=pData.tmp, family='binomial'))$coef[2, c(1, 2, 4)], 
                                          warning = function(w) {
                                            logistf.fit.tmp <- logistf((outcome == 'TRUE') ~ PCSS2, data=pData.tmp)
                                            c( logistf.fit.tmp$coef[2], 
                                               sqrt( logistf.fit.tmp$var[2, 2] ), 
                                               logistf.fit.tmp$prob[2] ) %>% return
                                          })
      p.values[c('PCSS1.biv', 'PCSS2.biv'), ] <- tryCatch(summary(glm((outcome == 'TRUE') ~ PCSS1 + PCSS2, 
                                                                      data=pData.tmp, family='binomial' ) )$coef[2:3, c(1, 2, 4)], 
                                                          warning = function(w) {
                                                            logistf.fit.tmp <- logistf((outcome == 'TRUE') ~ PCSS1 + PCSS2, data=pData.tmp)
                                                            c( logistf.fit.tmp$coef[2:3], 
                                                               c( sqrt( logistf.fit.tmp$var[2, 2] ), 
                                                                  sqrt( logistf.fit.tmp$var[3, 3] ) ), 
                                                               logistf.fit.tmp$prob[2:3] )
                                                          })
    }
    
    effect.size <- cbind( effect.size, p.values[2:5, 1] )
    colnames( effect.size )[ncol( effect.size )] <- set
    effect.se <- cbind( effect.se, p.values[2:5, 2] )
    colnames( effect.se )[ncol( effect.se )] <- set
    write.csv( p.values, file=paste0( output.dir, var, '/', var, '_', set, '.csv' ) )
  }
  
  # meta-analysis
  results.rma <- matrix( NA, 4, 4 )
  dimnames( results.rma ) <- list(c('PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv'), 
                                  c('Effect', 'se', 'pValue', 'I2'))
  pdf(paste0( output.dir, var, '/meta_', var, '.pdf' ), 4*2, 4*2)
  par(mfrow=c(2, 2))
  for (i in 1:4) {
    rma.fit <- rma(yi=effect.size[i, ], sei=effect.se[i, ], method='FE')
    results.rma[i, ] <- c(rma.fit$b, rma.fit$se, rma.fit$pval, rma.fit$I2)
    forest(rma.fit, slab=colnames(effect.size), xlab=rownames(results.rma)[i])
  }
  dev.off()
  write.csv(results.rma, file=paste0( output.dir, var, '/meta_', var, '.csv'))
  
  # meta-analysis with REML
  pdf(paste0(output.dir, var, '/meta_', var, '_RE.pdf'), 4*2, 4*2)
  par(mfrow=c(2, 2))
  for (i in 1:4) {
    rma.fit <- rma(yi=effect.size[i, ], sei=effect.se[i, ], method='REML')
    results.rma[i, ] <- c(rma.fit$b, rma.fit$se, rma.fit$pval, rma.fit$I2)
    forest(rma.fit, slab=colnames( effect.size ), xlab=rownames(results.rma)[i])
  }
  dev.off()
  write.csv(results.rma, file=paste0(output.dir, var, '/meta_', var, '_RE.csv'))
}



# disease free survival ---------------------------------------------------
library(survival)
var <- 'dfs'
sets.avail <- c("GSE12945_eset",
                "GSE14333_eset",
                "GSE39582_eset",
                "GSE33113_eset",
                "GSE17536_eset",
                "TCGA.COAD_eset",
                "TCGA.RNASeqV2_eset")
dir.create(paste0( output.dir, var, '/' ), showWarnings=F)

effect.size <- matrix( NA, 4, 0 )
rownames( effect.size ) <- c('PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv')
effect.se <- effect.size

for(set in sets.avail) {
  p.values <- matrix(NA, 5, 3)
  rownames(p.values) <- c('n', 'PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv')
  colnames(p.values) <- c('effect_size', 'sd', 'p_value')
  
  pData.tmp <- pData( eSets[[set]] )[pData( eSets[[set]] )[, 'sample_type']=='tumor', ]
  p.values['n', 1] <- sum( !is.na( pData.tmp[, 'dfs_status'] ) )
  surv.outcome <- Surv(pData.tmp$days_to_recurrence_or_death, 
                       pData.tmp$dfs_status=='deceased_or_recurrence')
  coxph.fit <- coxph(surv.outcome ~ PCSS1, data = pData.tmp)
  p.values['PCSS1.uni', ] <- summary(coxph.fit)$coef[c(1, 3, 5)]
  coxph.fit <- coxph(surv.outcome ~ PCSS2, data = pData.tmp)
  p.values['PCSS2.uni', ] <- summary(coxph.fit)$coef[c(1, 3, 5)]
  coxph.fit <- coxph(surv.outcome ~ PCSS1 + PCSS2, data = pData.tmp)
  p.values[c('PCSS1.biv', 'PCSS2.biv'), ] <- summary(coxph.fit)$coef[, c(1, 3, 5)]
  
  effect.size <- cbind( effect.size, p.values[2:5, 1] )
  colnames( effect.size )[ncol( effect.size )] <- set
  effect.se <- cbind( effect.se, p.values[2:5, 2] )
  colnames( effect.se )[ncol( effect.se )] <- set
  write.csv( p.values, file=paste0( output.dir, var, '/', var, '_', set, '.csv' ) )
}


# meta-analysis
results.rma <- matrix( NA, 4, 4 )
dimnames( results.rma ) <- list(c('PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv'), 
                                c('Effect', 'se', 'pValue', 'I2'))
pdf(paste0( output.dir, var, '/meta_', var, '.pdf' ), 4*2, 4*2)
par(mfrow=c(2, 2))
for (i in 1:4) {
  rma.fit <- rma(yi=effect.size[i, ], sei=effect.se[i, ], method='FE')
  results.rma[i, ] <- c(rma.fit$b, rma.fit$se, rma.fit$pval, rma.fit$I2)
  forest(rma.fit, slab=colnames(effect.size), xlab=rownames(results.rma)[i])
}
dev.off()
write.csv(results.rma, file=paste0( output.dir, var, '/meta_', var, '.csv'))

# meta-analysis with REML
pdf(paste0(output.dir, var, '/meta_', var, '_RE.pdf'), 4*2, 4*2)
par(mfrow=c(2, 2))
for (i in 1:4) {
  rma.fit <- rma(yi=effect.size[i, ], sei=effect.se[i, ], method='REML')
  results.rma[i, ] <- c(rma.fit$b, rma.fit$se, rma.fit$pval, rma.fit$I2)
  forest(rma.fit, slab=colnames( effect.size ), xlab=rownames(results.rma)[i])
}
dev.off()
write.csv(results.rma, file=paste0(output.dir, var, '/meta_', var, '_RE.csv'))


# CRIS subtypes with GSE76402 and GSE59857 ---------------------------------------------
load("data/CRIS/GSE76402_eset.RData")
load("data/CRIS/GSE59857_eset.RData")
eSets <- list()
for (set in c(setNames, "GSE76402_eset", "GSE59857_eset")) {
  eset.tmp <- get(set)
  pdata.tmp <- pData(eset.tmp)
  pdata.tmp <- pdata.tmp %>% mutate(CRISA = (CRIS_label == 'CRIS-A'),
                                    CRISB = (CRIS_label == 'CRIS-B'),
                                    CRISC = (CRIS_label == 'CRIS-C'),
                                    CRISD = (CRIS_label == 'CRIS-D'),
                                    CRISE = (CRIS_label == 'CRIS-E'))
  pData(eset.tmp) <- pdata.tmp
  eSets[[set]] <- eset.tmp
}

output.dir <- 'results/clinical/'
dir.create(output.dir, showWarnings = F, recursive = T)

vars <- c("CRISA", "CRISB", "CRISC", "CRISD", "CRISE")
vars.contrasts <- c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE")
names(vars.contrasts) <- vars

for (var in vars) {
  dir.create(paste0(output.dir, var, '/'), showWarnings=F)
  
  effect.size <- matrix(NA, 4, 0)
  rownames(effect.size) <- c('PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv')
  effect.se <- effect.size
  
  for(set in c(setNames, "GSE76402_eset", "GSE59857_eset")) {
    p.values <- matrix(NA, 5, 3)
    rownames(p.values) <- c('n', 'PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv')
    colnames(p.values) <- c('effect_size', 'sd', 'p_value')
    
    pData.tmp <- pData(eSets[[set]])[pData( eSets[[set]] )[, 'sample_type']=='tumor', ]
    outcome <- (pData.tmp[, var] == vars.contrasts[var])
    
    count.table <- table(outcome)
    if(length(count.table) <= 1) next # all outcome values are the same, not able to evaluate
    else {
      p.values['n', 1] <- sum( !is.na( pData.tmp[, var] ) )
      p.values['PCSS1.uni', ] <- tryCatch(summary(glm((outcome == 'TRUE') ~ PCSS1, 
                                                      data=pData.tmp, family='binomial'))$coef[2, c(1, 2, 4)], 
                                          warning = function(w) {
                                            logistf.fit.tmp <- logistf((outcome == 'TRUE') ~ PCSS1, data=pData.tmp)
                                            c( logistf.fit.tmp$coef[2], 
                                               sqrt( logistf.fit.tmp$var[2, 2] ), 
                                               logistf.fit.tmp$prob[2] ) %>% return
                                          })
      p.values['PCSS2.uni', ] <- tryCatch(summary(glm((outcome == 'TRUE') ~ PCSS2, 
                                                      data=pData.tmp, family='binomial'))$coef[2, c(1, 2, 4)], 
                                          warning = function(w) {
                                            logistf.fit.tmp <- logistf((outcome == 'TRUE') ~ PCSS2, data=pData.tmp)
                                            c( logistf.fit.tmp$coef[2], 
                                               sqrt( logistf.fit.tmp$var[2, 2] ), 
                                               logistf.fit.tmp$prob[2] ) %>% return
                                          })
      p.values[c('PCSS1.biv', 'PCSS2.biv'), ] <- tryCatch(summary(glm((outcome == 'TRUE') ~ PCSS1 + PCSS2, 
                                                                      data=pData.tmp, family='binomial' ) )$coef[2:3, c(1, 2, 4)], 
                                                          warning = function(w) {
                                                            logistf.fit.tmp <- logistf((outcome == 'TRUE') ~ PCSS1 + PCSS2, data=pData.tmp)
                                                            c( logistf.fit.tmp$coef[2:3], 
                                                               c( sqrt( logistf.fit.tmp$var[2, 2] ), 
                                                                  sqrt( logistf.fit.tmp$var[3, 3] ) ), 
                                                               logistf.fit.tmp$prob[2:3] )
                                                          })
    }
    
    effect.size <- cbind( effect.size, p.values[2:5, 1] )
    colnames( effect.size )[ncol( effect.size )] <- set
    effect.se <- cbind( effect.se, p.values[2:5, 2] )
    colnames( effect.se )[ncol( effect.se )] <- set
    write.csv( p.values, file=paste0( output.dir, var, '/', var, '_', set, '.csv' ) )
  }
  
  # meta-analysis
  results.rma <- matrix( NA, 4, 4 )
  dimnames( results.rma ) <- list(c('PCSS1.uni', 'PCSS2.uni', 'PCSS1.biv', 'PCSS2.biv'), 
                                  c('Effect', 'se', 'pValue', 'I2'))
  pdf(paste0( output.dir, var, '/meta_', var, '.pdf' ), 4*2, 4*2)
  par(mfrow=c(2, 2))
  for (i in 1:4) {
    rma.fit <- rma(yi=effect.size[i, ], sei=effect.se[i, ], method='FE')
    results.rma[i, ] <- c(rma.fit$b, rma.fit$se, rma.fit$pval, rma.fit$I2)
    forest(rma.fit, slab=colnames(effect.size), xlab=rownames(results.rma)[i])
  }
  dev.off()
  write.csv(results.rma, file=paste0( output.dir, var, '/meta_', var, '.csv'))
  
  # meta-analysis with REML
  pdf(paste0(output.dir, var, '/meta_', var, '_RE.pdf'), 4*2, 4*2)
  par(mfrow=c(2, 2))
  for (i in 1:4) {
    rma.fit <- rma(yi=effect.size[i, ], sei=effect.se[i, ], method='REML')
    results.rma[i, ] <- c(rma.fit$b, rma.fit$se, rma.fit$pval, rma.fit$I2)
    forest(rma.fit, slab=colnames( effect.size ), xlab=rownames(results.rma)[i])
  }
  dev.off()
  write.csv(results.rma, file=paste0(output.dir, var, '/meta_', var, '_RE.csv'))
}
