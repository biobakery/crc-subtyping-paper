rm(list = ls())
source("src/misc.R")
library(tidyverse)
library(CMSclassifier)
library(org.Hs.eg.db)

load('data/eSets/setNames.RData')
for (set in setNames) {
  load(set %>% paste0('data/eSets/new/', ., '.RData'))
}

data("centroids")
df_mapping_SSP <- AnnotationDbi::select(org.Hs.eg.db, 
                                        rownames(centroids), 
                                        columns = "SYMBOL")

# Panel A -----------------------------------------------------------------
# data frame for silhouette width across all samples
df_silWidth <- lapply(setNames, function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  exprs_tmp <- exprs_tmp %>% 
    rmNaInf %>% 
    apply(1, function(x) x - mean(x)) %>% 
    t
  
  exprs_tmp_SSP <- exprs_tmp[intersect(rownames(exprs_tmp), df_mapping_SSP$SYMBOL),]
  
  distance <- dist(exprs_tmp_SSP %>% t, 'eucl') %>% as.matrix
  ind_lbl <- pdata_tmp$cms_label_SSP != "unlabeled"
  
  rbind(data.frame(sw=silhouette( pdata_tmp$cms_label_SSP[ind_lbl] %>% 
                                    factor %>% 
                                    as.numeric, 
                                  dist=distance[ind_lbl, ind_lbl] %>% 
                                    as.dist )[, 'sil_width'], 
                   study=set,
                   samples='labeled'),
        data.frame(sw=silhouette( pdata_tmp$cms_label_SSP %>% 
                                    factor %>% 
                                    as.numeric, 
                                  dist=distance %>% 
                                    as.dist )[, 'sil_width'], 
                   study=set,
                   samples='all')
  ) %>% return
}) %>% Reduce("rbind", .)

load("data/eSets/trainingSetNames.RData")
setNames.crcsc <- c('GSE42284', 'GSE33113', 'GSE39582', 'GSE35896', 'GSE13067', 
                    'GSE13294', 'GSE14333', 'GSE17536', 'GSE20916', 'GSE2109', 
                    'GSE2109', 'TCGA.RNASeqV2', 'TCGA.COAD')
df_silWidth <- df_silWidth %>% 
  mutate(training = ifelse(study %in% trainingSetNames, 
                           "training",
                           "validation"),
         study = study %>% 
           gsub("_eset", "", ., fixed = T),
         inCRCSC = ifelse(study %in% setNames.crcsc,
                          "Used in CMS classification paper",
                          "Others"))
df_silWidth <- df_silWidth %>% 
  mutate(training = training %>% 
           factor(levels = c("training", "validation")),
         samples = samples %>% 
           factor(levels = c("labeled", "all")),
         inCRCSC = inCRCSC %>% 
           factor(levels = c("Used in CMS classification paper", "Others")))
df_avgWidth_ranked <- df_silWidth %>% 
  filter(samples == "all") %>% 
  group_by(study) %>% 
  summarise(avg_width = mean(sw, na.rm = T),
            training = unique(training)) %>% 
  arrange(training, desc(avg_width))
df_silWidth <- df_silWidth %>% 
  mutate(study = study %>% 
           factor(levels = df_avgWidth_ranked$study))

# boxplots
cols.tick <- rep( 'black', 18 )
cols.tick[levels(df_silWidth$study) %in% setNames.crcsc] <- 'red'
pA <- ggplot(data=df_silWidth, aes(study, sw)) + 
  geom_boxplot(aes(color=inCRCSC, fill=samples)) + 
  stat_summary(aes(fill=samples), fun.y=mean, geom="point", shape=18, 
               size=5, show.legend = FALSE,
                position = position_dodge(width = 0.75))  + 
  theme_bw() + 
  labs(x = "Study", y = "Silhouette Width") + 
  theme(axis.text.x = element_text(angle = 30, hjust=0.9, color=cols.tick)) + 
  scale_colour_manual(values=c('Used in CMS classification paper'='red', 'Others'='black'), 
                      name='Datasets' ) +
  scale_fill_manual(values=c('labeled' = 'white', 'all' = 'grey'), name = 'Samples included')
# add in annotation horizontal and vertical ines
pA <- pA + geom_hline(aes(yintercept=0.25), color='grey', linetype='dashed') + 
  geom_hline(aes(yintercept=0.5), color='grey', linetype='dashed') + 
  geom_hline(aes(yintercept=0.75), color='grey', linetype='dashed') + 
  geom_vline(aes(xintercept=8.5), color='black', linetype='solid')
# add in annotationt texts
pA <- pA + annotate('text', x=17.55, y=0.8, color='grey', label ='Strong Support') + 
  annotate('text', x=17.4, y=0.55, color='grey', label ='Moderate Support') + 
  annotate('text', x=16.85, y=0.3, color='grey', label ='Weak and could be artificial') + 
  annotate('text', x=17.05, y=-0.3, color='grey', label ='No substantial structure') + 
  annotate('text', x=7.3, y=0.5, color='black', label ='Training') + 
  annotate('text', x=10, y=0.5, color='black', label ='Validation') + 
  theme(legend.box = "horizontal", legend.justification=c(0,1), legend.position=c(0,1), 
        legend.background = element_rect( colour='black' ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  guides(color = guide_legend(order = 1))


# panel B -----------------------------------------------------------------
set.toplot <- 'GSE17536_eset'
eSet <- get(set.toplot)
pdata <- pData(eSet)
eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
exprs_tmp <- exprs(eSet_tmp)
pdata_tmp <- pData(eSet_tmp)

exprs_tmp <- exprs_tmp %>% 
  rmNaInf %>% 
  apply(1, function(x) x - mean(x)) %>% 
  t
exprs_tmp_SSP <- exprs_tmp[intersect(rownames(exprs_tmp), df_mapping_SSP$SYMBOL),]
fit.prcomp <- prcomp(t(exprs_tmp_SSP))
PC.toplot <- fit.prcomp$x[, 1:4]

pdata_tmp <- pdata_tmp %>% 
  cbind(PC.toplot) %>% 
  mutate(cms_label_SSP = cms_label_SSP %>% 
           recode("unlabeled" = "not labeled"))
colors <- gg_color_hue(4)
colors.toplot <- c(colors, 'grey')
names(colors.toplot) <- c(paste0('CMS', 1:4), 'not labeled')
pB1 <- pdata_tmp %>% 
  ggplot(aes(x = PC1, y = PC2, color = cms_label_SSP)) +
  geom_point(aes(alpha = (cms_label_SSP == 'not labeled'))) +
  scale_color_manual(values = colors.toplot, name = 'CMS subtype') +
  scale_alpha_manual(values = c('TRUE' = 0.5, 'FALSE' = 1), guide = F) +
  theme_bw() +
  theme(legend.direction = 'horizontal',
        legend.position = c(0, 1),
        legend.justification=c(0,1),
        legend.background = element_rect(colour='black')) +
  ggtitle('GSE17536, PC1 vs. PC2') +
  xlab(paste0('PC1 (', round(fit.prcomp$sdev[1]^2 / sum(fit.prcomp$sdev^2) * 100, 2), '%)')) +
  ylab(paste0('PC2 (', round(fit.prcomp$sdev[2]^2 / sum(fit.prcomp$sdev^2) * 100, 2), '%)'))
pB2 <- pdata_tmp %>% 
  ggplot(aes(x = PC3, y = PC4, color = cms_label_SSP)) +
  geom_point(aes(alpha = (cms_label_SSP == 'not labeled'))) +
  scale_color_manual(values = colors.toplot, name = 'CMS subtype') +
  scale_alpha_manual(values = c('TRUE' = 0.5, 'FALSE' = 1), guide = F) +
  theme_bw() +
  theme(legend.direction = 'horizontal',
        legend.position = c(0, 1),
        legend.justification=c(0,1),
        legend.background = element_rect(colour='black')) +
  ggtitle('GSE17536, PC3 vs. PC4') +
  xlab(paste0('PC3 (', round(fit.prcomp$sdev[3]^2 / sum(fit.prcomp$sdev^2) * 100, 2), '%)')) +
  ylab(paste0('PC4 (', round(fit.prcomp$sdev[4]^2 / sum(fit.prcomp$sdev^2) * 100, 2), '%)'))


# arranging plots ---------------------------------------------------------
library(cowplot)
# plot <- ggdraw() +
#     draw_plot(pA, 0, 0.5, 1, .5) +
#     draw_plot(pB1, 0, 0, .5, .5) +
#     draw_plot(pB2, .5, 0, .5, .5) +
#     draw_plot_label(c("A", "B", ""), c(0, 0, 0.5), c(1, 0.5, 0.5))
# ggsave(plot, file = "manuscript/Figure2.pdf", width = 12, height = 10)
ggsave(pA, file = "manuscript/Figure2.pdf", width = 12, height = 5)
