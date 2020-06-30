rm(list=ls())
source("src/misc.R")
library(tidyverse)

load('data/eSets/setNames.RData')

vars <- c('msi', 'summarygrade', 'summarystage', 'summarylocation')
df.toplot <- lapply(vars, function(var) {
  df.results <- read.csv(paste0('results/model_comparison/CRIS/', var, '.csv'))
  # df.results <- df.results[df.results$PCSS1.sd < 100, ]
  return(data.frame(df.results, Variable = var))
}) %>% Reduce('rbind', .) %>% 
  mutate(
    Variable = Variable %>% 
      recode("msi" = 'MSI status', 
             "summarygrade" = 'grade', 
             "summarystage" = 'stage', 
             "summarylocation" = 'tumor location')
  )

df.toplot <- df.toplot %>% 
  gather(key = Model,
         value = log10_pvalue,
         LR.cont, LR.disc) 
df.toplot$Model[df.toplot$Model == 'LR.cont'] <- 'score'
df.toplot$Model[df.toplot$Model == 'LR.disc'] <- 'subtype'
df.toplot$Model <- factor(df.toplot$Model, levels = c('subtype', 'score'))
p <- df.toplot %>%
  ggplot(aes( x= Model, y = -log10(log10_pvalue) )) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2)) +
  facet_wrap(~Variable, nrow = 1, scales = 'free') +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed') +
  ylab('-log-10 p value')

ggsave(p, file = "manuscript/SuppFigure_modelComp.pdf", width = 6, height = 4)
