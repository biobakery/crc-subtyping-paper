rm(list = ls())
source("src/misc.R")
library(tidyverse)

load('data/eSets/setNames.RData')
load("results/discrete_evaluation/all_silWidth.RData")
load("results/discrete_evaluation/silWidth_additional.RData")
df_plot <- rbind(df_silWidth, df_silWidth_additional)

# panel A -----------------------------------------------------------------
pA <- df_plot %>% 
  mutate(study = study %>% 
           gsub("_eset", "", ., fixed = T) %>% 
           factor(levels = c(
             "GSE76402",
             "GSE59857",
             gsub("_eset", "", setNames, fixed = T)
           ))) %>% 
  filter(samples == "all",
         metric == "Euclidean",
         subtype == "CRIS") %>% 
  ggplot(aes(study, sw)) + 
  geom_boxplot() + 
  stat_summary(fun.y=mean, geom="point", shape=18, 
               size=5, show.legend = FALSE,
               position = position_dodge(width = 0.75))  + 
  theme_bw() + 
  labs(x = "Study", y = "Silhouette Width") + 
  theme(axis.text.x = element_text(angle = 30, hjust=0.9)) +
  geom_hline(aes(yintercept=0.25), color='grey', linetype='dashed') +
  annotate('text', x=16.85, y=0.3, color='grey', label ='Weak and could be artificial') + 
  annotate('text', x=17.05, y=-0.3, color='grey', label ='No substantial structure')

# panel B -----------------------------------------------------------------
pB <- df_plot %>% 
  filter(samples == "all") %>% 
  group_by(study, metric, subtype) %>% 
  summarise(`Average Silhouette Width` = mean(sw)) %>% 
  ggplot(aes(metric, `Average Silhouette Width`)) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  theme_bw() + 
  labs(x = "Dissimilary Metric", y = "Average Silhouette Width") + 
  geom_hline(aes(yintercept=0.25), linetype='dashed') +
  geom_text(data = data.frame(x = 0.5, y = 0.27, label = "Weak and could be artificial", 
                              subtype = c("CMS", "CRIS")), 
            aes(x = x, y = y, label = label, alpha = subtype), hjust = 0) + 
  geom_text(data = data.frame(x = 0.5, y = 0.23, label = "No substantial structure", 
                              subtype = c("CMS", "CRIS")), 
            aes(x = x, y = y, label = label, alpha = subtype), hjust = 0) +
  scale_alpha_manual(values = c("CMS" = 1, "CRIS" = 0), guide = FALSE) +
  facet_wrap(~subtype, nrow = 1)


# plotting ----------------------------------------------------------------
library(cowplot)
p <- plot_grid(pA, pB, labels = c("A", "B"), ncol = 1)
ggsave(p, file = "manuscript/SuppFigure_SilWidth.pdf", width = 8, height = 8)
