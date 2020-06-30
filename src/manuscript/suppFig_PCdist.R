rm(list=ls())
source("src/misc.R")
library(tidyverse)
library(CMSclassifier)
library(org.Hs.eg.db)
library(cowplot)

load("data/eSets/setNames.RData")
for (set in setNames) {
  load(set %>% paste0("data/eSets/new/", ., ".RData"))
}

colors.toplot <- gg_color_hue(4)
colors.toplot <- c(colors.toplot, "grey")
names(colors.toplot) <- c(paste0("CMS", 1:4), "unlabeled")

data("centroids")
df_mapping_SSP <- select(org.Hs.eg.db, 
                         rownames(centroids), 
                         columns = "SYMBOL")

l_p <- list()
i <- 1
for (set.toplot in setNames) {
  exprs.toplot <- set.toplot %>% get %>% exprs
  genes.common <- intersect(rownames(exprs.toplot), df_mapping_SSP$SYMBOL)
  pdata.toplot <- pData(get(set.toplot))
  exprs.toplot <- exprs.toplot[genes.common, ]
  exprs.toplot <- exprs.toplot[!apply(is.na(exprs.toplot), 1, any), ]
  exprs.toplot <- exprs.toplot[!apply(abs(exprs.toplot) == Inf, 1, any), ]
  fit.prcomp <- prcomp(t(exprs.toplot))
  PC.toplot <- fit.prcomp$x[, 1:4]
  
  p_tmp <- cbind(pdata.toplot, PC.toplot) %>% 
    ggplot(aes(x = PC1, y = PC2, color = cms_label_SSP)) +
    geom_point(aes(alpha = cms_label_SSP == "unlabeled")) +
    scale_color_manual(values = colors.toplot, name = "CMS subtype") +
    scale_alpha_manual(values = c("TRUE" = 0.5, "FALSE" = 1), guide = F) +
    theme_bw() + ggtitle(paste0(gsub("_eset", "", set.toplot, fixed = T), ", PC1 vs. PC2")) +
    theme(
      legend.direction = "horizontal",
      legend.position = c(0, 0),
      legend.justification=c(0,0),
      legend.background = element_rect(colour="black")) +
    xlab(paste0("PC1 (", round(fit.prcomp$sdev[1]^2 / sum(fit.prcomp$sdev^2) * 100, 2), "%)")) +
    ylab(paste0("PC2 (", round(fit.prcomp$sdev[2]^2 / sum(fit.prcomp$sdev^2) * 100, 2), "%)"))
  df_mock <- data.frame(PC1 = rep(0, 5),
                        PC2 = rep(0, 5),
                        cms_label_SSP = factor(c("CMS1", "CMS2", "CMS3", "CMS4", "unlabeled"),
                                               levels = c("CMS1", "CMS2", "CMS3", "CMS4", "unlabeled")))
  p_tmp <- p_tmp + geom_point(aes(x = PC1, y = PC2, color = cms_label_SSP), alpha = 0, data = df_mock)
  if(i > 1) p_tmp <- p_tmp + guides(color = FALSE)
  if(i == 1) p_tmp <- p_tmp + guides(color = guide_legend(nrow = 2, byrow = TRUE))
  l_p[[i]] <- p_tmp
  i <- i + 1
  p_tmp <- cbind(pdata.toplot, PC.toplot) %>% 
    ggplot(aes(x = PC3, y = PC4, color = cms_label_SSP)) +
    geom_point(aes(alpha = cms_label_SSP == "unlabeled")) +
    scale_color_manual(values = colors.toplot, name = "CMS subtype") +
    scale_alpha_manual(values = c("TRUE" = 0.5, "FALSE" = 1), guide = F) +
    theme_bw() + ggtitle(paste0(gsub("_eset", "", set.toplot, fixed = T), ", PC3 vs. PC4")) +
    theme(
      legend.direction = "horizontal",
      legend.position = c(0, 1),
      legend.justification=c(0,1),
      legend.background = element_rect(colour="black")) +
    xlab(paste0("PC3 (", round(fit.prcomp$sdev[3]^2 / sum(fit.prcomp$sdev^2) * 100, 2), "%)")) +
    ylab(paste0("PC4 (", round(fit.prcomp$sdev[4]^2 / sum(fit.prcomp$sdev^2) * 100, 2), "%)"))
  if(i > 1) p_tmp <- p_tmp + guides(color = FALSE)
  l_p[[i]] <- p_tmp
  i <- i + 1
}
pdf("manuscript/SuppFigure_PCdist.pdf", 
    width = 24, height = 30)
plot_grid(plotlist = l_p, ncol = 6)
dev.off()