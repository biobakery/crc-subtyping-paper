rm(list=ls())
source("src/misc.R")
library(Biobase)
library(tidyverse)
library(survival)
library(metafor)
library(logistf)
library(tidyr)

load("data/eSets/setNames.RData")
for (set in setNames) {
  load(paste0("data/eSets/new/", set, ".RData"))
}


# dist’n of PCSS scores ---------------------------------------------------
df.results <- setNames %>% lapply(function(set) {
  df.meta <- set %>% get %>% pData
  df.toreturn <- lapply(c('PCSS1', 'PCSS2'), function(score) {
    score.tmp <- df.meta[, score]
    score.tmp <- score.tmp - mean(score.tmp, na.rm = T)
    df.toreturn <- data.frame(mean = tapply(score.tmp, df.meta$cms_label_RF, 
                                            function(x) c(mean(x, na.rm=T))),
                              sd = tapply(score.tmp, df.meta$cms_label_RF, 
                                          function(x) c(sd(x, na.rm=T)))) 
    colnames(df.toreturn) <- paste0(c('mean_', 'sd_'), score)
    return(df.toreturn)
  } ) %>% Reduce('cbind', .)
  df.toreturn$subtype <- rownames(df.toreturn)
  return(data.frame(df.toreturn, set = set))
}) %>% Reduce('rbind', .)

colors <- gg_color_hue(4)
colors.toplot <- c(colors, 'grey')
names(colors.toplot) <- c(paste0('CMS', 1:4), 'unlabeled')

pA <- df.results %>% 
  filter(!is.na(sd_PCSS1),
         subtype != 'not labeled',
         mean_PCSS1 < 2.5) %>% 
  ggplot(aes( x = mean_PCSS1, y = mean_PCSS2, color = subtype)) +
  geom_point() +
  geom_errorbar(aes( x = mean_PCSS1, ymin = mean_PCSS2 - sd_PCSS2, ymax = mean_PCSS2 + sd_PCSS2)) +
  geom_errorbarh(aes( y = mean_PCSS2, xmin = mean_PCSS1 - sd_PCSS1, xmax = mean_PCSS1 + sd_PCSS1)) +
  scale_color_manual(values = colors.toplot, name = '') +
  # scale_alpha_manual(values = c('TRUE' = 0.5, 'FALSE' = 1), guide = FALSE) +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(vjust=1.5,face='italic'),
        axis.title.y=element_text(vjust=1.5,face='italic'),
        legend.position="top",
        panel.background=element_rect(fill = 'white', color = 'black'),
        # panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  coord_cartesian(xlim = c(-2, 2.5), ylim = c(-2.8, 2)) +
  xlab('PCSS1') + ylab('PCSS2') 
ggsave(pA, file = "manuscript/Figure1_toyFigure4A.pdf",
       width = 3.75,
       height = 3.75)

vars <- c('msi', 'summarygrade', 'summarystage', 'summarylocation')
df.toplot <- lapply(vars, function(var) {
  df.results <- read.csv(paste0('results/model_comparison/CMS/', var, '.csv'))
  # df.results <- df.results[df.results$PCSS1.sd < 100, ]
  return(data.frame(df.results, Variable = var))
}) %>% Reduce('rbind', .) %>% 
  mutate(
    Variable = Variable %>% 
      recode("msi" = 'MSI status', 
             "summarygrade" = 'grade', 
             "summarystage" = 'stage', 
             "summarylocation" = 'location', 
             "dfs" = 'DFS')
  )
df.toplot <- df.toplot %>% 
  gather(key = Model,
         value = log10_pvalue,
         LR.cont, LR.disc) 
df.toplot$Model[df.toplot$Model == 'LR.cont'] <- 'score'
df.toplot$Model[df.toplot$Model == 'LR.disc'] <- 'subtype'
df.toplot$Model[df.toplot$Model == 'score' & df.toplot$Variable != 'MSI status'] <- ''
df.toplot$Model[df.toplot$Model == 'subtype' & df.toplot$Variable != 'MSI status'] <- ' '
df.toplot$Model <- factor(df.toplot$Model, levels = c('subtype', 'score', ' ', ''))
df.toplot$Variable <- factor(df.toplot$Variable, levels = c('MSI status', 'grade', 'stage', 'location', 'DFS'))
pB <- df.toplot %>%
  ggplot(aes( x= Model, y = -log(log10_pvalue) )) +
  geom_boxplot() +
  facet_wrap(~Variable, nrow = 1, scales = 'free') +
  theme_bw() +
  ylab('-log-10 p value') +
  theme(axis.text.y=element_blank(),
        axis.text.x=element_text(size = 8, face='italic'),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(vjust=1.5,face='italic'),
        legend.position="none",
        panel.background=element_rect(fill = 'white', color = 'black'),
        # panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
ggsave(pB, file = "manuscript/Figure1_toyFigure4B.pdf",
       width = 4.8,
       height = 2.8)