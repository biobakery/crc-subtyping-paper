rm(list = ls())
source("src/misc.R")
library(metafor)
library(tidyverse)

vars <- c("CMS1", "CMS2", "CMS3", "CMS4", "CMS5",
          "msi", "summarygrade", "summarylocation", "summarystage",
          "dfs")
df_loo <- lapply(vars, function(var) {
  sets <- list.files(paste0("results/clinical/", var),
                     patter = "_eset.csv") %>% 
    gsub(paste0(var, "_"), "", ., fixed = T) %>% 
    gsub(".csv", "", ., fixed = T)
  if("NHS.HPFS_eset" %in% sets) {
    # data frame for summary statistics
    df_ss <- lapply(sets, function(set) {
      df_result <- paste0("results/clinical/", var, "/",
                          var, "_", set, ".csv") %>% 
        read.csv
      data.frame(n = df_result$effect_size[1],
                 effect_size = c(df_result$effect_size[2],
                                 df_result$effect_size[3]),
                 sd = c(df_result$sd[2],
                        df_result$sd[3]),
                 score = c("PCSS1", "PCSS2"),
                 set = set)
    }) %>% 
      Reduce("rbind", .)
    df_originalI2 <- df_ss %>% 
      group_by(score) %>% 
      do(rma.uni(yi = effect_size,
                 sei = sd,
                 method = "FE",
                 data = .)$I2 %>% 
           tibble(`Original I2` = .))
    df_I2 <- sets %>% 
      lapply(function(set_loo) {
        df_ss %>% 
          filter(set != set_loo) %>% 
          group_by(score) %>% 
          do(rma.uni(yi = effect_size,
                     sei = sd,
                     method = "FE",
                     data = .)$I2 %>% 
               data.frame(I2 = .,
                          set_exclude = set_loo)) %>% 
          ungroup
      }) %>% 
      Reduce("rbind", .) %>% 
      mutate(Variable = var)
    df_I2 <- df_I2 %>% 
      left_join(df_originalI2, 
                by = "score")
    return(df_I2)
  }
  return(tibble(score = NA,
                I2 = NA,
                set_exclude = NA,
                Variable = NA,
                `Original I2` = NA))
}) %>% Reduce("rbind", .)
figure <- df_loo %>% 
  filter(!is.na(score)) %>% 
  mutate(`Excluded Study` =
           set_exclude %>% 
           recode_factor(
             "NHS.HPFS_eset" = "NHS/HPPFS",
             .default = "others",
             .ordered = TRUE
           )) %>% 
  ggplot(aes(x = `Excluded Study`,
             y = I2)) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_grid(score ~ Variable) +
  geom_hline(aes(yintercept = `Original I2`),
             linetype = "dashed",
             color = "red") +
  theme_bw()
ggsave("response_letter/I2_comparison.pdf", 
       figure, 
       width = 12, height = 4)
