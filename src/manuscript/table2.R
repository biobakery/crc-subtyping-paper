library(tidyverse)
vars <- c("CMS1", "CMS2", "CMS3", "CMS4", "msi", "summarylocation",
          "summarystage", "summarygrade", "dfs")
vars %>% lapply(function(var) {
  df_tmp <- read_csv(paste0("results/clinical/", var, "/meta_", var, ".csv"))
  df_tmp <- df_tmp[1:2, ] %>% 
    mutate(Variable = var %>% 
             recode(
               "CMS1" = "CMS1 subtype", 
               "CMS2" = "CMS2 subtype", 
               "CMS3" = "CMS3 subtype", 
               "CMS4" = "CMS4 subtype", 
               "msi" = "MSI", 
               "summarylocation" = "Right location",
               "summarystage" = "Late stage", 
               "summarygrade" = "High grade", 
               "dfs" = "Disease recurrence or death"
             ) %>% c(""),
           `Continuous score` = c("PCSS1", "PCSS2"),
           `Effect size` = Effect,
           `p-value` = pValue) %>% 
      dplyr::select(Variable,
             `Continuous score`,
             `Effect size`,
             `p-value`
             )
}) %>% 
  Reduce("rbind", .) %>% 
  write_csv("manuscript/table2.csv")

vars %>% lapply(function(var) {
  sets <- list.files(paste0("results/clinical/", var, "/"), 
                     pattern = paste0("^", var, "_")) %>% 
    gsub(paste0(var, "_"), "", ., fixed = T) %>% 
    gsub("_eset.csv", "", ., fixed = T)
  df_I2 <- read_csv(paste0("results/clinical/", var, "/meta_", var, ".csv"))
  sets %>% lapply(function(set){
    df_tmp <- read_csv(paste0("results/clinical/", var, "/", var, "_",
                              set, "_eset.csv"))
    tibble(
      Variable = ifelse(set == sets[1], var, ""),
      Study = set,
      n = df_tmp$effect_size[1],
      `PCSS1 effect size` = df_tmp$effect_size[2],
      `PCSS1 standard error` = df_tmp$sd[2],
      `PCSS2 effect size` = df_tmp$effect_size[3],
      `PCSS2 standard error` = df_tmp$sd[3],
      `PCSS1 I2` = ifelse(set == sets[1], df_I2$I2[1] %>% round(digits = 4), ""),
      `PCSS2 I2` = ifelse(set == sets[1], df_I2$I2[2] %>% round(digits = 4), "")
    )
  }) %>% 
    Reduce("rbind", .)
}) %>% 
  Reduce("rbind", .) %>% 
  mutate(Variable = Variable %>% 
           recode(
             "CMS1" = "CMS1 subtype", 
             "CMS2" = "CMS2 subtype", 
             "CMS3" = "CMS3 subtype", 
             "CMS4" = "CMS4 subtype", 
             "msi" = "MSI", 
             "summarylocation" = "Right location",
             "summarystage" = "Late stage", 
             "summarygrade" = "High grade", 
             "dfs" = "Disease recurrence or death"
           )) %>% 
  write_csv("manuscript/table2_sheet2.csv")

c("msi", "summarylocation", "summarystage", "summarygrade") %>% 
  lapply(function(var) {
    read_csv(paste0("results/model_comparison/CMS/", var, ".csv")) %>% 
      filter(!is.na(`AIC.full`)) %>% 
      mutate(Variable = ifelse(X1 == X1[1], var, ""),
             Study = X1,
             `p value (subtype)` = LR.disc,
             `p value (continuous scores)` = LR.cont,
             `AIC (full model)` = AIC.full,
             `AIC (subtype)` = AIC.disc,
             `AIC (continuous scores)` = AIC.cont) %>% 
      dplyr::select(Variable,
             Study,
             `p value (subtype)`,
             `p value (continuous scores)`,
             `AIC (full model)`,
             `AIC (subtype)`,
             `AIC (continuous scores)`
             )
    }) %>% 
  Reduce("rbind", .) %>% 
  write_csv("manuscript/table2_sheet3.csv")