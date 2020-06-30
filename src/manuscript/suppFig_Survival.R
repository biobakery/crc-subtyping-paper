rm(list=ls())
source("src/misc.R")
library(tidyverse)

load("data/eSets/setNames.RData")

df_CMS <- read_csv("results/model_comparison/CMS/dfs.csv")
p1 <- df_CMS %>% 
  mutate(study = study %>% 
           gsub("_eset", "", ., fixed = T)) %>% 
  mutate(study = study %>% 
           factor(levels = c("RE Model", rev(study[-1])))) %>% 
  ggplot(aes(x = study, y = beta)) +
  geom_point(aes(shape = (study == "RE Model"),
                 size = 1 / sd)) +
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 16)) +
  geom_errorbar(aes(ymin = beta - 1.96*sd,
                    ymax = beta + 1.96*sd),
                width = 0.25) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("log Hazard Ratio") +
  xlab("Study") +
  theme_bw() +
  ggtitle("CMS4 differentiation") +
  coord_flip() +
  theme(legend.position="none")
df_CRIS <- read_csv("results/model_comparison/CRIS/dfs.csv")
p2 <- df_CRIS %>% 
  mutate(study = study %>% 
           gsub("_eset", "", ., fixed = T)) %>% 
  mutate(study = study %>% 
           factor(levels = c("RE Model", rev(study[-1])))) %>% 
  ggplot(aes(x = study, y = beta)) +
  geom_point(aes(shape = (study == "RE Model"),
                 size = 1 / sd)) +
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 16)) +
  geom_errorbar(aes(ymin = beta - 1.96*sd,
                    ymax = beta + 1.96*sd),
                width = 0.25) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("log Hazard Ratio") +
  xlab("Study") +
  theme_bw() +
  ggtitle("CRIS-B differentiation") +
  coord_flip() +
  theme(legend.position="none")

library(cowplot)
plot_grid(p1, p2, labels = c("A", "B"), nrow = 1) %>% 
  ggsave(
    file = "manuscript/SuppFigure_Survival.pdf",
    width = 8, height = 3
  )
ggforest = function(x, bug){
    require("ggplot2")
    # Function to convert REM results in `rma`-format into a data.frame
    rma2df = function(x){
      rbind(
        data.frame(Study = c("MSH", "PRISM", "RISK biopsy", "RISK stool"), LogFC = x$yi,
                   CILB=x$yi - 2*sqrt(x$vi),
                   CIUB=x$yi + 2*sqrt(x$vi),
                   p = x$pval,
                   stringsAsFactors = FALSE),
        data.frame(Study = "RE Model", LogFC = x$b, CILB=x$ci.lb, CIUB=x$ci.ub,
                   p = x$pval,
                   stringsAsFactors = FALSE)
      ) %>% mutate(Study = factor(Study, levels = c("RE Model",
                                                    "RISK stool",
                                                    "RISK biopsy",
                                                    "PRISM",
                                                    "MSH")),
                   aggregate = ifelse(Study == "RE Model", "all", "individual"))
    }
    remresdf = rma2df(x)
    remresdf <- transform(remresdf, interval = CIUB - CILB)
    remresdf <- transform(remresdf, RelConf = 1/interval)
    p = ggplot(remresdf,
               aes(LogFC, Study, xmax=CIUB, xmin=CILB)) +
      # coord_cartesian(xlim=c(-2, 2)) +
      scale_alpha_discrete(range = c(0.2, 1)) +
      geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) +
      geom_errorbarh(alpha=0.5, color="black", height = 0.3) +
      geom_point(aes(size = RelConf, shape = aggregate)) +
      geom_point(data = subset(remresdf, Study=="RE Model"), aes(shape = aggregate), size=7) +
      scale_size(range = c(2, 5), guide=FALSE) +
      theme_bw() +
      theme(text = element_text(size=14)) +
      xlab("Differnetial Abundance (IBD vs. control)") +
      ggtitle(bug) +
      scale_shape_manual(values = c(15, 16), guide = F) +
      theme(plot.title=element_text(face="italic"))
    return(p)
  }