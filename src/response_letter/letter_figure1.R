rm(list = ls())
source("src/misc.R")
library(tidyverse)
library(CMSclassifier)
library(org.Hs.eg.db)

load('data/eSets/setNames.RData')
for (set in setNames) {
  load(set %>% paste0('data/eSets/new/', ., '.RData'))
}

# Panel A -----------------------------------------------------------------
# PCSS1 vs PCSS2 for GSE17536
set.toplot <- 'GSE17536_eset'
eSet <- get(set.toplot)
pdata <- pData(eSet)
eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
exprs_tmp <- exprs(eSet_tmp)
pdata_tmp <- pData(eSet_tmp)

pdata_tmp <- pdata_tmp %>% 
  mutate(cms_label_SSP = cms_label_SSP %>% 
           recode("unlabeled" = "not labeled"))
colors <- gg_color_hue(4)
colors.toplot <- c(colors, 'grey')
names(colors.toplot) <- c(paste0('CMS', 1:4), 'not labeled')
pA <- pdata_tmp %>% 
  ggplot(aes(x = PCSS1, y = PCSS2, color = cms_label_SSP)) +
  geom_point(aes(alpha = (cms_label_SSP == 'not labeled'))) +
  scale_color_manual(values = colors.toplot, name = 'CMS subtype') +
  scale_alpha_manual(values = c('TRUE' = 0.5, 'FALSE' = 1), guide = F) +
  theme_bw() +
  theme(legend.direction = 'horizontal',
        legend.position = c(0, 1),
        legend.justification=c(0,1),
        legend.background = element_rect(colour='black')) +
  ggtitle('GSE17536, PCSS1 vs. PCSS2')


# panel B -----------------------------------------------------------------
load("data/eSets/trainingSetNames.RData")
df_silWidth <- lapply(setNames, function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  distance <- dist(pdata_tmp[, c("PCSS1", "PCSS2")], 'eucl') %>% as.matrix
  ind_lbl <- pdata_tmp$cms_label_SSP != "unlabeled"
  
  rbind(data.frame(sw=silhouette((pdata_tmp$cms_label_SSP[ind_lbl]) %>% 
                                   factor %>% 
                                   as.numeric, 
                                 dist=distance[ind_lbl, ind_lbl] %>% 
                                   as.dist )[, 'sil_width'], 
                   study=set,
                   samples='labeled'),
        data.frame(sw=silhouette((pdata_tmp$cms_label_SSP) %>% 
                                   factor %>% 
                                   as.numeric, 
                                 dist=distance %>% 
                                   as.dist )[, 'sil_width'], 
                   study=set,
                   samples='all')
  ) %>% return
}) %>% Reduce("rbind", .)
df_silWidth <- df_silWidth %>% 
  mutate(study = study %>% 
           gsub("_eset", "", ., fixed = T),
         samples = samples %>% 
           factor(levels = c("labeled", "all")))
df_avgWidth_ranked <- df_silWidth %>% 
  filter(samples == "all") %>% 
  group_by(study) %>% 
  summarise(avg_width = mean(sw, na.rm = T),
            sd_width = sd(sw, na.rm = T)) %>% 
  arrange(desc(avg_width))
df_avgWidth_ranked <- df_avgWidth_ranked %>% 
  mutate(study = study %>% 
           factor(levels = df_avgWidth_ranked$study))
# boxplots
pB <- df_avgWidth_ranked %>% 
  ggplot(aes(study, avg_width)) + 
  geom_point(shape=18, 
             size=5) +
  geom_errorbar(aes(ymin = avg_width - sd_width,
                    ymax = avg_width + sd_width),
                alpha = 0.5) +
  theme_bw() + 
  labs(x = "Study", y = "Average Silhouette Width") + 
  theme(axis.text.x = element_text(angle = 30, hjust=0.9)) +
  geom_hline(aes(yintercept=0.25), color='black', linetype='dashed') +
  annotate('text', x=16, y=0.2, color='black', label ='No substantial structure') 

# arranging plots ---------------------------------------------------------
library(cowplot)
plot <- ggdraw() +
  draw_plot(pA, 0, 0, 0.4, 1) +
  draw_plot(pB, 0.4, 0, 0.6, 1) +
  draw_plot_label(c("A", "B"), c(0, 0.4), c(1, 1))
ggsave(plot, file = "response_letter/letter_figure1.pdf", width = 14, height = 6)
