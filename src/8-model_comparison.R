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
  eSets[[set]] <- eset.tmp
}

vars <- c("msi", "summarylocation", "summarygrade", "summarystage")
vars.contrasts <- c("MSI", "r", "high", "late")
names(vars.contrasts) <- vars

# vs. CMS subtypes --------------
output.dir <- "results/model_comparison/CMS/"
dir.create(output.dir, showWarnings = F, recursive = T)

# binary outcomes
for (var in vars) {
  results.table <- matrix(NA, length(setNames), 5)
  dimnames(results.table) <- list(setNames, 
                                  c("AIC.full", "AIC.disc", "LR.disc",
                                    "AIC.cont", "LR.cont"))
  for( set in setNames ) {
    pData.tmp <- pData(eSets[[set]]) %>% 
      subset(sample_type == "tumor")
    outcome <- (pData.tmp[, var] == vars.contrasts[var])
    
    count.table <- table(outcome)
    if(length(count.table) <= 1) next # all outcome values are the same, not able to evaluate
    else {
      model.all <- tryCatch(glm((outcome == "TRUE") ~ factor(cms_label_SSP) + PCSS1 + PCSS2,
                                data = pData.tmp, family = "binomial", 
                                control = list(maxit = 1000)),
                            warning = function(w) NULL)
      model.disc <- tryCatch(glm((outcome == "TRUE") ~ factor(cms_label_SSP),
                                 data = pData.tmp, family = "binomial", 
                                 control = list(maxit = 1000)),
                             warning = function(w) NULL)
      model.cont <- tryCatch(glm((outcome == "TRUE") ~ PCSS1 + PCSS2,
                                 data = pData.tmp, family = "binomial", 
                                 control = list(maxit = 1000)),
                             warning = function(w) NULL)
      if(any(c(is.null(model.all), is.null(model.disc), is.null(model.cont)))) next
      else {
        results.table[set, c("AIC.full", "AIC.disc", "AIC.cont")] <- AIC(model.all, model.disc, model.cont)[, 2]
        results.table[set, c("LR.disc", "LR.cont")] <- c(anova(model.all, model.disc, test="LRT")[2, 5],
                                                         anova(model.all, model.cont, test="LRT")[2, 5])
      }
    }
  }
  write.csv(results.table, file = paste0(output.dir, 
                                         var, ".csv"))
}


# DFS
var <- 'dfs'
sets.avail <- c("GSE12945_eset",
                "GSE14333_eset",
                "GSE39582_eset",
                "GSE33113_eset",
                "GSE17536_eset",
                "TCGA.COAD_eset",
                "TCGA.RNASeqV2_eset")
dir.create( paste0( output.dir, var, '/' ), showWarnings=F )

results.table <- lapply(sets.avail, function(set) {
  df.tmp <- get(set) %>% pData %>% 
    filter(sample_type == "tumor",
           !is.na(dfs_status),
           !is.na(days_to_recurrence_or_death)) %>% 
    filter(cms_label_SSP == "CMS4") %>% 
    mutate(dfs_status = as.character(dfs_status)) %>% 
    mutate(dfs_status = ifelse(days_to_recurrence_or_death > 365*5,
                               "living_norecurrence", 
                               dfs_status),
           days_to_recurrence_or_death = ifelse(days_to_recurrence_or_death > 365*5,
                                                365*5, 
                                                days_to_recurrence_or_death)
    )
  df.tmp <- df.tmp %>% 
    mutate(Strata = 
             ifelse(PCSS1 > quantile(PCSS1, 0.75) |
                      PCSS2 > quantile(PCSS1, 0.75),
                    "High PCSS1/PCSS2",
                    "Others") %>% 
             factor(levels = c("Others", "High PCSS1/PCSS2")))
  coxfit <- coxph(Surv(days_to_recurrence_or_death, 
                       dfs_status=='deceased_or_recurrence') ~ Strata,
                  data = df.tmp)
  data.frame(beta = summary(coxfit)$coef[1], sd = summary(coxfit)$coef[3], study = set)
}) %>% 
  Reduce("rbind", .) %>% 
  filter(sd < 100)
fit_rma <- rma(results.table$beta, sei = results.table$sd,  method = "FE")
results.table <- rbind(data.frame(beta = fit_rma$beta,
                                  sd = fit_rma$se,
                                  study = "RE Model"),
                       results.table)
write.csv(results.table, file = paste0(output.dir, 
                                       var, ".csv"),
          row.names = F)

# vs. CRIS subtypes --------------
output.dir <- "results/model_comparison/CRIS/"
dir.create(output.dir, showWarnings = F, recursive = T)

# binary outcomes
for (var in vars) {
  results.table <- matrix(NA, length(setNames), 5)
  dimnames(results.table) <- list(setNames, 
                                  c("AIC.full", "AIC.disc", "LR.disc",
                                    "AIC.cont", "LR.cont"))
  for( set in setNames ) {
    pData.tmp <- pData(eSets[[set]]) %>% 
      subset(sample_type == "tumor")
    outcome <- (pData.tmp[, var] == vars.contrasts[var])
    
    count.table <- table(outcome)
    if(length(count.table) <= 1) next # all outcome values are the same, not able to evaluate
    else {
      model.all <- tryCatch(glm((outcome == "TRUE") ~ factor(CRIS_label) + PCSS1 + PCSS2,
                                data = pData.tmp, family = "binomial", 
                                control = list(maxit = 1000)),
                            warning = function(w) NULL)
      model.disc <- tryCatch(glm((outcome == "TRUE") ~ factor(CRIS_label),
                                 data = pData.tmp, family = "binomial", 
                                 control = list(maxit = 1000)),
                             warning = function(w) NULL)
      model.cont <- tryCatch(glm((outcome == "TRUE") ~ PCSS1 + PCSS2,
                                 data = pData.tmp, family = "binomial", 
                                 control = list(maxit = 1000)),
                             warning = function(w) NULL)
      if(any(c(is.null(model.all), is.null(model.disc), is.null(model.cont)))) next
      else {
        results.table[set, c("AIC.full", "AIC.disc", "AIC.cont")] <- AIC(model.all, model.disc, model.cont)[, 2]
        results.table[set, c("LR.disc", "LR.cont")] <- c(anova(model.all, model.disc, test="LRT")[2, 5],
                                                         anova(model.all, model.cont, test="LRT")[2, 5])
      }
    }
  }
  write.csv(results.table, file = paste0(output.dir, 
                                         var, ".csv"))
}

# DFS
var <- 'dfs'
sets.avail <- c("GSE12945_eset",
                "GSE14333_eset",
                "GSE39582_eset",
                "GSE33113_eset",
                "GSE17536_eset",
                "TCGA.COAD_eset",
                "TCGA.RNASeqV2_eset")
dir.create( paste0( output.dir, var, '/' ), showWarnings=F )
results.table <- lapply(sets.avail, function(set) {
  df.tmp <- get(set) %>% pData %>% 
    filter(sample_type == "tumor",
           !is.na(dfs_status),
           !is.na(days_to_recurrence_or_death)) %>% 
    filter(CRIS_label == "CRIS-B") %>% 
    mutate(dfs_status = as.character(dfs_status)) %>% 
    mutate(dfs_status = ifelse(days_to_recurrence_or_death > 365*5,
                               "living_norecurrence", 
                               dfs_status),
           days_to_recurrence_or_death = ifelse(days_to_recurrence_or_death > 365*5,
                                                365*5, 
                                                days_to_recurrence_or_death)
    )
  df.tmp <- df.tmp %>% 
    mutate(Strata = 
             ifelse(PCSS1 > quantile(PCSS1, 0.5) & PCSS2 < quantile(PCSS2, 0.5),
                    "High PCSS1, low PCSS2",
                    "Others") %>% 
             factor(levels = c("Others", "High PCSS1, low PCSS2")))
  coxfit <- coxph(Surv(days_to_recurrence_or_death, 
                       dfs_status=='deceased_or_recurrence') ~ Strata,
                  data = df.tmp)
  return(data.frame(beta = summary(coxfit)$coef[1], sd = summary(coxfit)$coef[3], study = set))
}) %>% 
  Reduce("rbind", .) %>% 
  filter(sd < 100)
fit_rma <- rma(results.table$beta, sei = results.table$sd,  method = "FE")
results.table <- rbind(data.frame(beta = fit_rma$beta,
                                  sd = fit_rma$se,
                                  study = "RE Model"),
                       results.table)
write.csv(results.table, file = paste0(output.dir, 
                                       var, ".csv"),
          row.names = F)