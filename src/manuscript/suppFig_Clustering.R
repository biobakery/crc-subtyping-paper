rm(list = ls())
source("src/misc.R")
library(tidyverse)

load('data/eSets/setNames.RData')
for(set in setNames) {
  load(paste0("results/clustering/df_", set, ".RData"))
}
df_all <- setNames %>% 
  lapply(function(set) paste0("df_", set) %>% get %>% data.frame(study = set)) %>% 
  Reduce("rbind", .)
# p <- df_all %>% 
#   mutate(study = study %>% gsub("_eset", "", ., fixed = T),
#          metric = factor(metric, levels = c("Prediction Strength",
#                                             "Gap Statistic",
#                                             "Silhouette Width")),
#          method = factor(method, levels = c("k-medoid", "NMF", "Consensus Hierarchical"))) %>% 
#   ggplot(aes(x = K, y = statistics, shape = method)) +
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_line(aes(linetype = method), position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = statistics - se, ymax = statistics + se), width = 0.5, position = position_dodge(width = 0.5)) +
#   facet_grid(metric ~ study, scale = "free_y") + 
#   theme_bw()
ggsave(p, file = "manuscript/SuppFigure_Clustering.pdf", width = 26, height = 8)

for(set in c(setNames, "GSE59857_eset", "GSE76402_eset")) {
  load(paste0("results/clustering/CRIS/df_", set, ".RData"))
}
df_all_CRIS <- c(setNames, "GSE59857_eset", "GSE76402_eset") %>% 
  lapply(function(set) paste0("df_", set) %>% get %>% data.frame(study = set)) %>% 
  Reduce("rbind", .) %>% 
  mutate(`Genes used` = "CRIS")
p <- df_all %>% mutate(`Genes used` = "Top 3000") %>% 
  rbind(df_all_CRIS) %>% 
  mutate(study = study %>% 
           gsub("_eset", "", ., fixed = T) %>% 
           factor(levels = c(setNames, "GSE76402_eset", "GSE59857_eset") %>% 
                    gsub("_eset", "", ., fixed = T)) ,
         metric = factor(metric, levels = c("Prediction Strength",
                                            "Gap Statistic",
                                            "Silhouette Width")),
         method = factor(method, levels = c("k-medoid", "NMF", "Consensus Hierarchical"))) %>% 
  ggplot(aes(x = K, y = statistics, shape = method, color = `Genes used`)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = statistics - se, ymax = statistics + se), width = 0.5, position = position_dodge(width = 0.5)) +
  facet_grid(metric ~ study, scale = "free_y") + 
  theme_bw()
ggsave(p, file = "manuscript/SuppFigure_Clustering.pdf", width = 30, height = 8)
