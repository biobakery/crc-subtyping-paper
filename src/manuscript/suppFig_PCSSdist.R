rm(list = ls())
source("src/misc.R")
library(tidyverse)
library(CMSclassifier)
library(org.Hs.eg.db)

load("data/eSets/setNames.RData")
for (set in setNames) {
  load(set %>% paste0("data/eSets/new/", ., ".RData"))
}

# Panel A -----------------------------------------------------------------
# PCSS1 vs PCSS2 for GSE17536
set.toplot <- "GSE17536_eset"
df_plot <- setNames %>% lapply(function(set) {
  set %>% get %>% pData
}) %>%
  Reduce("rbind", .) %>% 
  filter(sample_type %in% "tumor")

colors <- gg_color_hue(4)
colors.toplot <- c(colors, "grey")
names(colors.toplot) <- c(paste0("CMS", 1:4), "unlabeled")

pA <- df_plot %>% 
  ggplot(aes(x = PCSS1, y = PCSS2)) +
  geom_point(alpha = 0.3) +
  geom_density_2d() +
  theme_bw()

pB <- df_plot %>% 
  ggplot(aes(x = PCSS1, y = PCSS2, color = cms_label_SSP)) +
  geom_point(aes(alpha = (cms_label_SSP == "not labeled"))) +
  scale_color_manual(values = colors.toplot, name = "CMS subtype") +
  scale_alpha_manual(values = c("TRUE" = 0.5, "FALSE" = 1), guide = F) +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.position = c(0, 1),
        legend.justification=c(0,1),
        legend.background = element_rect(colour="black"))


# panel B -----------------------------------------------------------------
load("data/eSets/trainingSetNames.RData")
df_silWidth <- lapply(setNames, function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp)
  pdata_tmp <- pData(eSet_tmp)
  
  distance <- dist(pdata_tmp[, c("PCSS1", "PCSS2")], "eucl") %>% as.matrix
  ind_lbl <- pdata_tmp$cms_label_SSP != "unlabeled"
  
  rbind(data.frame(sw=silhouette((pdata_tmp$cms_label_SSP[ind_lbl]) %>% 
                                   factor %>% 
                                   as.numeric, 
                                 dist=distance[ind_lbl, ind_lbl] %>% 
                                   as.dist )[, "sil_width"], 
                   study=set,
                   samples="labeled"),
        data.frame(sw=silhouette((pdata_tmp$cms_label_SSP) %>% 
                                   factor %>% 
                                   as.numeric, 
                                 dist=distance %>% 
                                   as.dist )[, "sil_width"], 
                   study=set,
                   samples="all")
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
pC <- df_avgWidth_ranked %>% 
  ggplot(aes(study, avg_width)) + 
  geom_point(shape=18, 
             size=5) +
  geom_errorbar(aes(ymin = avg_width - sd_width,
                    ymax = avg_width + sd_width),
                alpha = 0.5) +
  theme_bw() + 
  labs(x = "Study", y = "Average Silhouette Width") + 
  theme(axis.text.x = element_text(angle = 30, hjust=0.9)) +
  geom_hline(aes(yintercept=0.25), color="black", linetype="dashed") +
  annotate("text", x=16, y=0.2, color="black", label ="No substantial structure") 

# arranging plots ---------------------------------------------------------
library(cowplot)
plot <- ggdraw() +
  draw_plot(pA, 0, 0.5, 0.5, 0.5) +
  draw_plot(pB, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(pC, 0, 0, 1, 0.5) +
  draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0), c(1, 1, 0.5))
ggsave(plot, file = "manuscript/SuppFigure_PCSSdist.pdf", width = 10.5, height = 10.5)
