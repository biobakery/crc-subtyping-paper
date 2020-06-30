rm(list = ls())
source("src/misc.R")
library(GEOquery)
library(CRISclassifier)
library(org.Hs.eg.db)


# GSE59857 ----------------------------------------------------------------
dir.create("data/CRIS/", recursive = T, showWarnings = F)
gse59857 <- getGEO(GEO = "GSE59857", destdir = "data/CRIS/")
gpl10558 <- getGEO(GEO = "GPL10558", destdir = "data/CRIS/")

# expression matrix
exprs <- exprs(gse59857[[1]])
mapping <- Table(gpl10558)
mapping <- mapping %>% 
  dplyr::filter(ID %in% rownames(exprs))
exprs_collapsed <- sapply(setdiff(unique(mapping$Entrez_Gene_ID), NA), 
                          function(gene) {
                            exprs_tmp <- exprs[mapping$ID[mapping$Entrez_Gene_ID %in% gene], , 
                                               drop = F] %>% log
                            mean_tmp <- apply(exprs_tmp, 1, mean, na.rm = T)
                            return(exprs_tmp[order(mean_tmp, decreasing = T)[1], ])
                          }) %>% t
rownames(exprs_collapsed) <- setdiff(unique(mapping$Entrez_Gene_ID), NA)

data("features")
df_mapping_CRIS <- select(org.Hs.eg.db, 
                          features$`Gene ID` %>% as.character, 
                          keytype = "SYMBOL",
                          columns = "ENTREZID")
exprs_features <- exprs_collapsed %>%
  as.data.frame %>% 
  rownames_to_column("ENTREZID") %>% 
  right_join(df_mapping_CRIS, by = "ENTREZID") %>% 
  filter(!is.na(ENTREZID)) %>% 
  column_to_rownames("SYMBOL") %>% 
  dplyr::select(-ENTREZID) %>% 
  as.matrix
exprs_features %>%
  data.frame(symbol = rownames(exprs_features),
             .) %>%
  write.table(file = "results/discrete_evaluation/additional/GSE59857.txt",
              quote = F,
              sep = "\t",
              row.names = F)

# pdata
set.seed(1)
cris_classifier("results/discrete_evaluation/additional/GSE59857.txt",
                "results/discrete_evaluation/additional/GSE59857_classify.txt",
                nresmpl = 1)
pdata <- read.table("results/discrete_evaluation/additional/GSE59857_classify.txt_sample_info.txt",
                    sep = "\t",
                    header = T) %>% 
  mutate(CRIS_label = predict.label %>% 
           recode(`1` = "CRIS-A",
                  `2` = "CRIS-B",
                  `3` = "CRIS-C",
                  `4` = "CRIS-D",
                  `5` = "CRIS-E") %>% 
           factor(levels = c(
             "CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E", "unlabeled"
           ))) %>% 
  column_to_rownames("sample.names")

# assign continuous score
avg_loadings <- read.csv("data/loadings/avg_loadings.csv",
                         row.names = 1) %>% 
  as.matrix
df_mapping_score <- select(org.Hs.eg.db, 
                           rownames(avg_loadings), 
                           keytype = "SYMBOL",
                           columns = "ENTREZID") %>% 
  filter(!duplicated(SYMBOL))
exprs_score <- exprs_collapsed %>% 
  as.data.frame %>% 
  rownames_to_column(("ENTREZID")) %>% 
  right_join(df_mapping_score, by = "ENTREZID") %>% 
  filter(!is.na(ENTREZID)) %>% 
  column_to_rownames("SYMBOL") %>% 
  dplyr::select(-ENTREZID) %>% 
  as.matrix
exprs_score <- exprs_score %>% rmNaInf
exprs_score <- apply(exprs_score, 1, function(x) x - mean(x)) %>% t
genesCommon <- intersect(rownames(exprs_score), rownames(avg_loadings))
pcss <- t(exprs_score[genesCommon, ]) %*% avg_loadings[genesCommon, ]
pcss <- (t(pcss) / apply(pcss, 2, sd)) %>% t
colnames(pcss) <- paste0( 'PCSS', 1:4 ) 
pdata[, paste0("PCSS", 1:4)] <- pcss
pdata$sample_type <- "tumor"

GSE59857_eset <- ExpressionSet(assayData = exprs_features %>% as.matrix,
                               phenoData = AnnotatedDataFrame(pdata))
save(GSE59857_eset, file = "data/CRIS/GSE59857_eset.RData")


# GSE76402 ----------------------------------------------------------------
gse76402 <- getGEO(GEO = "GSE76402", destdir = "data/CRIS/")
gpl10558 <- getGEO(GEO = "GPL10558", destdir = "data/CRIS/")

# expression matrix
exprs <- exprs(gse76402[[1]])
mapping <- Table(gpl10558)
mapping <- mapping %>% 
  dplyr::filter(ID %in% rownames(exprs))
exprs_collapsed <- sapply(setdiff(unique(mapping$Entrez_Gene_ID), NA), 
                          function(gene) {
                            exprs_tmp <- exprs[mapping$ID[mapping$Entrez_Gene_ID %in% gene], , 
                                               drop = F] %>% log
                            mean_tmp <- apply(exprs_tmp, 1, mean, na.rm = T)
                            return(exprs_tmp[order(mean_tmp, decreasing = T)[1], ])
                          }) %>% t
rownames(exprs_collapsed) <- setdiff(unique(mapping$Entrez_Gene_ID), NA)
exprs_features <- exprs_collapsed %>%
  as.data.frame %>% 
  rownames_to_column("ENTREZID") %>% 
  right_join(df_mapping_CRIS, by = "ENTREZID") %>% 
  filter(!is.na(ENTREZID)) %>% 
  column_to_rownames("SYMBOL") %>% 
  dplyr::select(-ENTREZID) %>% 
  as.matrix

# pdata
pdata <- pData(gse76402[[1]]) %>% 
  mutate(CRIS_label = characteristics_ch1.1 %>% 
           recode("cris class: CRISA" = "CRIS-A",
                  "cris class: CRISB" = "CRIS-B",
                  "cris class: CRISC" = "CRIS-C",
                  "cris class: CRISD" = "CRIS-D",
                  "cris class: CRISE" = "CRIS-E",
                  "cris class: Not assigned" = "unlabeled") %>% 
           factor(levels = c(
             "CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E", "unlabeled"
           ))) %>% 
  column_to_rownames("geo_accession")

# assign continuous score
avg_loadings <- read.csv("data/loadings/avg_loadings.csv",
                         row.names = 1) %>% 
  as.matrix
df_mapping_score <- select(org.Hs.eg.db, 
                           rownames(avg_loadings), 
                           keytype = "SYMBOL",
                           columns = "ENTREZID") %>% 
  filter(!duplicated(SYMBOL))
exprs_score <- exprs_features
exprs_score <- exprs_score %>% rmNaInf
exprs_score <- apply(exprs_score, 1, function(x) x - mean(x)) %>% t
genesCommon <- intersect(rownames(exprs_score), rownames(avg_loadings))
pcss <- t(exprs_score[genesCommon, ]) %*% avg_loadings[genesCommon, ]
pcss <- (t(pcss) / apply(pcss, 2, sd)) %>% t
colnames(pcss) <- paste0( 'PCSS', 1:4 ) 
pdata[, paste0("PCSS", 1:4)] <- pcss
pdata$sample_type <- "tumor"

GSE76402_eset <- ExpressionSet(assayData = exprs_features %>% as.matrix,
                               phenoData = AnnotatedDataFrame(pdata))
save(GSE76402_eset, file = "data/CRIS/GSE76402_eset.RData")
