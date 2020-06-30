rm(list=ls())
source("src/misc.R")
library(Biobase)
library(tidyverse)
library(survival)
library(metafor)
library(logistf)

load('data/eSets/setNames.RData')
for (set in setNames) {
  load(set %>% paste0('data/eSets/new/', ., '.RData'))
}


# dist’n of PCSS scores ---------------------------------------------------
df.results <- setNames %>% lapply(function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  pdata_tmp$study <- set
  return(pdata_tmp)
}) %>% Reduce('rbind', .)
df.results <- df.results %>% 
  mutate(cms_label_SSP = cms_label_SSP %>% 
           recode("unlabeled" = "not labeled"))
df.results <- df.results %>% 
  group_by(study, cms_label_SSP) %>% 
  dplyr::summarise(mean_PCSS1 = mean(PCSS1),
                   mean_PCSS2 = mean(PCSS2),
                   sd_PCSS1 = sd(PCSS1),
                   sd_PCSS2 = sd(PCSS2))

colors <- gg_color_hue(4)
colors.toplot <- c(colors, 'grey')
names(colors.toplot) <- c(paste0('CMS', 1:4), 'not labeled')

pA <- ggplot(df.results, 
             aes(x = mean_PCSS1, y = mean_PCSS2, color = cms_label_SSP)) +
  geom_point() +
  geom_errorbar(aes(x = mean_PCSS1, 
                    ymin = mean_PCSS2 - sd_PCSS2, 
                    ymax = mean_PCSS2 + sd_PCSS2)) +
  geom_errorbarh(aes( y = mean_PCSS2, 
                      xmin = mean_PCSS1 - sd_PCSS1, 
                      xmax = mean_PCSS1 + sd_PCSS1)) +
  scale_color_manual(values = colors.toplot, name = "CMS Subtype") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('PCSS1') + ylab('PCSS2') +
  coord_cartesian(xlim = c(-2, 2.5)) +
  theme(legend.direction = "horizontal", legend.justification=c(0,1), legend.position=c(0,1), 
        legend.background = element_rect(colour='black'))

# panel for discrete subtypes AIC likelihood ------------------------------
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
             "summarylocation" = 'tumor location', 
             "dfs" = 'DFS')
  )
# df.toplot %>% 
#   ggplot(aes(x = log10(LR.disc), y = log10(LR.cont))) +
#   geom_point() +
#   facet_wrap(~Variable, scales = "free") +
#   geom_abline(intercept = 0, slope = 1) +
#   coord_fixed() +
#   theme_bw()
df.toplot <- df.toplot %>% 
  gather(key = Model,
         value = log10_pvalue,
         LR.cont, LR.disc) 
df.toplot$Model[df.toplot$Model == 'LR.cont'] <- 'score'
df.toplot$Model[df.toplot$Model == 'LR.disc'] <- 'subtype'
df.toplot$Model <- factor(df.toplot$Model, levels = c('subtype', 'score'))
pB <- df.toplot %>%
  ggplot(aes( x= Model, y = -log10(log10_pvalue) )) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2)) +
  facet_wrap(~Variable, nrow = 1, scales = 'free') +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed') +
  ylab('-log-10 p value')


# panel for survival ------------------------------------------------------
library(survival)
library(survminer)

df.tmp <- lapply(c("GSE12945_eset",
                   "GSE14333_eset",
                   "GSE39582_eset",
                   "GSE33113_eset",
                   "GSE17536_eset",
                   "TCGA.COAD_eset",
                   "TCGA.RNASeqV2_eset"
),
function(set) {
  pData(get(set)) %>% 
    subset(sample_type == "tumor" &
             !is.na(dfs_status) &
             !is.na(days_to_recurrence_or_death)) %>%
    data.frame(study = set)
}) %>% 
  Reduce("rbind", .) %>% 
  mutate(dfs_status = as.character(dfs_status)) %>% 
  mutate(dfs_status = ifelse(days_to_recurrence_or_death > 365*5,
                             "living_norecurrence", 
                             dfs_status),
         days_to_recurrence_or_death = ifelse(days_to_recurrence_or_death > 365*5,
                                              365*5, 
                                              days_to_recurrence_or_death)
  )
colors_cms <- c(gg_color_hue(4))
names(colors_cms) <- c(paste0("CMS", 1:4))
pC1 <- ggsurvplot(survfit(Surv(days_to_recurrence_or_death, 
                               dfs_status=='deceased_or_recurrence') ~ cms_label_SSP,
                          data = df.tmp %>% subset(cms_label_SSP != "unlabeled")), 
                  legend.labs = names(colors_cms),
                  legend.title = "CMS Subtype")$plot +
  scale_color_manual(values = colors_cms) +
  theme_bw() +
  theme(legend.direction = "horizontal", legend.justification=c(1,0), legend.position=c(1,0), 
        legend.background = element_rect( colour='black' ),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.2, "cm")) +
  coord_cartesian(xlim = c(0, 5*365))  +
  xlab('Days to recurrence or death')

df.tmp.CMS4 <- subset(df.tmp, cms_label_SSP == "CMS4")
df.tmp.CMS4 <- df.tmp.CMS4 %>% 
  group_by(study) %>% 
  mutate(Strata = 
           ifelse(PCSS1 > quantile(PCSS1, 0.75) |
                    PCSS2 > quantile(PCSS1, 0.75),
                  "High PCSS1/PCSS2",
                  "Others") %>% 
           factor(levels = c("Others", "High PCSS1/PCSS2"))) %>% 
  ungroup
colors_cms4 <- c("orchid1", "purple4")
names(colors_cms4) <- c("Others", "High PCSS1/PCSS2")
pvalue <- (coxph(Surv(days_to_recurrence_or_death, 
                      dfs_status=='deceased_or_recurrence') ~ Strata,
                 data = df.tmp.CMS4) %>% summary)$coef[1, 5]
pC2 <- ggsurvplot(survfit(Surv(days_to_recurrence_or_death, 
                               dfs_status=='deceased_or_recurrence') ~ Strata,
                          data = df.tmp.CMS4), 
                  legend.labs = c("Others", "High PCSS1/PCSS2"),
                  legend.title = "CMS4 samples")$plot +
  scale_color_manual(values = colors_cms4) +
  theme_bw() +
  theme(legend.direction = "horizontal", legend.justification=c(1,0), legend.position=c(1,0), 
        legend.background = element_rect( colour='black' ),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.2, "cm")) +
  annotate("text", x = 1800, y = 0.45, 
           label = paste0("Cox PH p = ",
                          round(pvalue, digits = 4)),
           hjust = 1) +
  coord_cartesian(xlim = c(0, 5*365))  +
  xlab('Days to recurrence or death')

# plotting the figure -----------------------------------------------------
library(cowplot)
plot <- ggdraw() +
  draw_plot(pA, 0, 0, 0.5, 1) +
  draw_plot(pB, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(pC1, 0.5, 0, 0.25, 0.5) +
  draw_plot(pC2, 0.75, 0, 0.25, 0.5) +
  draw_plot_label(c("A", "B", "C", ""), c(0, 0.5, 0.5, 0.75), c(1, 1, 0.5, 0.5))
ggsave(plot, file = "manuscript/Figure4.pdf", width = 15, height = 7.5)
