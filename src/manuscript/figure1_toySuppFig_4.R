rm( list=ls() )
library( Biobase )
library( ggplot2 )
library( cluster )
library( fpc )
source( 'src/source/arrange_ggplot2.R' )
data.dir <- "data/clean_new/"
esets.dir <- "data/clean_new/"

load( paste0(data.dir, "setNames.RData") )
load( paste0(data.dir, "trainingSetNames.RData") )
load('data/CRCSC/CRCSC_eset.RData')

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

eSets <- list()
for( setName in c('GSE13067_eset', 'GSE14095_eset', 'GSE14333_eset') ){
  load(setName %>% paste0('data/clean_new/CRCSC_filled/', ., '.RData'))
  eSets[[setName]] <- get( setName )
}

K_max <- 8
B <- 10
M <- 10
l_df_toplot <- lapply( names( eSets ), function( name_eset ) {
  eset_tmp <- eSets[[name_eset]]
  mat_exprs <- exprs( eset_tmp )
  genes.common <- intersect(rownames(mat_exprs), featureNames(CRCSC_eset))
  mat_exprs <- mat_exprs[genes.common, ]
  
  mat_exprs <- mat_exprs[order( apply( mat_exprs, 1, sd, na.rm=T ), 
                                decreasing = T )[1:1000], ]
  mat_exprs <- mat_exprs[!apply( mat_exprs, 1, function( x ) any( is.na( x ) ) ), ]
  mat_exprs <- mat_exprs - apply( mat_exprs, 1, median, na.rm=T )
  # dist_tmp <- dist( t( mat_exprs ), method = 'eucl' )
  predstr_tmp <- prediction.strength( t( mat_exprs ), Gmin=2, Gmax=K_max, M=M, clustermethod=claraCBI )
  df_toreturnPredStr <- data.frame( statistics=predstr_tmp$mean.pred[-1], 
                                    se=sapply( predstr_tmp$predcorr[2:K_max], sd ),
                                    K=2:K_max, 
                                    measure='Prediction Strength', 
                                    ref=NA )
  df_toreturn <- df_toreturnPredStr
  df_toreturn$set <- name_eset
  return( df_toreturn )
} ) 
colnames(perc.var) <- gsub( '_eset', '',  colnames(perc.var), fixed=T )
rownames(perc.var) <- paste0('PC', 1:20)
write.csv(perc.var, file = 'CRC Manuscript/GenomeBiology/Resubmission/SuppMaterials/SuppTable 1/SuppTable1.csv',
          quote = F)
