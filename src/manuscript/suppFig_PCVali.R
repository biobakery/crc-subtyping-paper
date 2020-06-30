rm( list=ls() )
library( Biobase )
library( ggplot2 )
library( cluster )
library( fpc )
library( tidyr )
library( pamr )
source( 'src/source/arrange_ggplot2.R' )
data.dir <- "data/clean_new/"
esets.dir <- "data/clean_new/"
dir_output <- 'Resubmission/SuppMaterials/SuppFigure New 2/'

load( paste0(data.dir, "setNames.RData") )
load( paste0(data.dir, "trainingSetNames.RData") )
load( paste0(data.dir, "topGenes.RData") )
load( 'data/clean/trainingSetNames.RData' )

eSets <- list()
for( setName in setNames ){
  load( paste0( esets.dir, setName, ".RData") )
  eSets[[setName]] <- get( setName )
}

eSets_sub <- eSets[unlist( lapply( eSets, function( eset ) {
  df_meta_tmp <- pData( eset )
  sum( df_meta_tmp$sample_type == 'adjacentnormal' ) > 8
} ) )]

# unlist( lapply( eSets, function( eset ) {
#   df_meta_tmp <- pData( eset )
#   sum( df_meta_tmp$sample_type == 'adjacentnormal' )
# } ) )

eSets_sub <- lapply( eSets_sub, function( eset ) {
  return( eset[, pData( eset )$sample_type %in% 'adjacentnormal'] )
} )

mat_avgLoading <- read.csv( 'results/PCA/loadings/avg_loadings.csv', header=T,
                            row.names = 1 )
# continuous score validation
{
mat_corMaxNormal <- t( sapply( eSets_sub, function( eset )  {
  mat_exprs <- exprs( eset )
  mat_exprs <- mat_exprs[rownames( mat_exprs ) != 'MEOX2', ]
  mat_exprs <- mat_exprs[apply(mat_exprs, 1, function(x){!any(is.na(x)|(x==Inf)|(x==-Inf))}) ,]
  
  gene_common <- intersect( rownames( mat_avgLoading ), rownames( mat_exprs ) )
  # print( length( gene_common ) )
  prcomRes <- prcomp( t( mat_exprs[gene_common, ] ) )
  loadings <- prcomRes$rotation[,1:8]
  matCor <- abs( cor( mat_avgLoading[gene_common, ], loadings[gene_common, ], 
                      use = 'pairwise.complete.obs' ) )
  return( apply( matCor, 1, max ) )
} ) )

mat_corMaxTraining <- t( sapply( eSets[trainingSetNames], function( eset )  {
  mat_exprs <- exprs( eset )
  mat_exprs <- mat_exprs[rownames( mat_exprs ) != 'MEOX2', ]
  mat_exprs <- mat_exprs[apply(mat_exprs, 1, function(x){!any(is.na(x)|(x==Inf)|(x==-Inf))}) ,]
  
  gene_common <- intersect( rownames( mat_avgLoading ), rownames( mat_exprs ) )
  # print( length( gene_common ) )
  prcomRes <- prcomp( t( mat_exprs[gene_common, ] ) )
  loadings <- prcomRes$rotation[,1:8]
  matCor <- abs( cor( mat_avgLoading[gene_common, ], loadings[gene_common, ], 
                      use = 'pairwise.complete.obs' ) )
  return( apply( matCor, 1, max ) )
} ) )

mat_corMaxVali <- t( sapply( eSets[setdiff( setNames, trainingSetNames )], function( eset )  {
  mat_exprs <- exprs( eset )
  mat_exprs <- mat_exprs[rownames( mat_exprs ) != 'MEOX2', ]
  mat_exprs <- mat_exprs[apply(mat_exprs, 1, function(x){!any(is.na(x)|(x==Inf)|(x==-Inf))}) ,]
  
  gene_common <- intersect( rownames( mat_avgLoading ), rownames( mat_exprs ) )
  # print( length( gene_common ) )
  prcomRes <- prcomp( t( mat_exprs[gene_common, ] ) )
  loadings <- prcomRes$rotation[,1:8]
  matCor <- abs( cor( mat_avgLoading[gene_common, ], loadings[gene_common, ], 
                      use = 'pairwise.complete.obs' ) )
  return( apply( matCor, 1, max ) )
} ) )


R <- 30
mat_exprsTest <- exprs( eSets[['GSE13294_eset']] )
gene_common <- intersect( rownames( mat_exprsTest ), rownames( mat_avgLoading ) )
mat_exprsTest <- mat_exprsTest[gene_common, ]
n <- ncol( mat_exprsTest )
mat_corMaxTest <- t( sapply( 1:R, function( r ) {
  mat_order <- t( sapply( 1:nrow( mat_exprsTest ), function( i ) return( sample.int( n, n ) ) ) )
  mat_exprsTmp <- t( sapply( 1:nrow( mat_exprsTest ), function( i ) return( mat_exprsTest[i, ][mat_order[i, ]] ) ) )
  rownames( mat_exprsTmp ) <- gene_common
  
  prcomRes <- prcomp( t( mat_exprsTmp ) )
  loadings <- prcomRes$rotation[,1:8]
  matCor <- abs( cor( mat_avgLoading[gene_common, ], loadings[gene_common, ], 
                      use = 'pairwise.complete.obs' ) )
  return( apply( matCor, 1, max ) )
} ) )

df_toplot <- rbind( data.frame( mat_corMaxTraining[, 1:2], type='Training datasets' ), 
                    data.frame( mat_corMaxVali[, 1:2], type='Validation datasets' ), 
                    data.frame( mat_corMaxNormal[, 1:2], type='Normal tissues' ),
                    data.frame( mat_corMaxTest[, 1:2], type='Randomized datasets' ) )
colnames( df_toplot )[1:2] <- c( 'PCSS1', 'PCSS2' )
df_toplot <- gather( df_toplot, key = 'score_name', value = 'max_correlation', PCSS1, PCSS2 )
df_toplot$type <- factor( df_toplot$type, levels=c( 'Training datasets', 'Validation datasets', 'Normal tissues', 'Randomized datasets' ) )
pdf( paste0( dir_output, 'figure.pdf' ), width=12, height=4 )
print( ggplot( df_toplot, aes( x=type, y=max_correlation ) ) + 
         geom_boxplot( ) +
         geom_hline( yintercept=0.5, color='red', linetype=2 ) +
         facet_grid( .~score_name ) +
         theme_bw( ) )
dev.off()
}

# version 2: randomly select PCs from training datasets
l_loadings <- lapply( trainingSetNames, function( set ) {
  df_loading <- read.csv( paste0( 'results/PCA/loadings/', set, '.csv' ), header=T, 
                          row.names=1 )
  return( df_loading )
} )
df_loadings <- Reduce( 'cbind', l_loadings )
gene_common <- intersect( rownames( df_loadings ), rownames( mat_avgLoading ) )
mat_corMaxTest2 <- abs( cor( mat_avgLoading[gene_common, ], 
                             df_loadings[gene_common, sample.int( ncol( df_loadings ), size=R )] ) )

df_toplot <- rbind( data.frame( mat_corMaxTraining[, 1:2], type='Training datasets' ), 
                    data.frame( mat_corMaxVali[, 1:2], type='Validation datasets' ), 
                    data.frame( mat_corMaxNormal[, 1:2], type='Normal tissues' ),
                    data.frame( t( mat_corMaxTest2[1:2, ] ), type='Randomized datasets' ) )
colnames( df_toplot )[1:2] <- c( 'PCSS1', 'PCSS2' )
df_toplot <- gather( df_toplot, key = 'score_name', value = 'max_correlation', PCSS1, PCSS2 )
df_toplot$type <- factor( df_toplot$type, levels=c( 'Training datasets', 'Validation datasets', 'Normal tissues', 'Randomized datasets' ) )
pdf( paste0( dir_output, 'figure_v2.pdf' ), width=12, height=4 )
print( ggplot( df_toplot, aes( x=type, y=max_correlation ) ) + 
         geom_boxplot( ) +
         geom_hline( yintercept=0.5, color='red', linetype=2 ) +
         facet_grid( .~score_name ) +
         theme_bw( ) )
dev.off()

# version 3: both random cases
df_toplot <- rbind( data.frame( mat_corMaxTraining[, 1:2], type='Training datasets' ), 
                    data.frame( mat_corMaxVali[, 1:2], type='Validation datasets' ), 
                    data.frame( mat_corMaxNormal[, 1:2], type='Normal tissues' ),
                    data.frame( t( mat_corMaxTest2[1:2, ] ), type='Random PCs' ),
                    data.frame( mat_corMaxTest[, 1:2], type='Permuted datasets' ) )
colnames( df_toplot )[1:2] <- c( 'PCSS1', 'PCSS2' )
df_toplot <- gather( df_toplot, key = 'score_name', value = 'max_correlation', PCSS1, PCSS2 )
df_toplot$type <- factor( df_toplot$type, levels=c( 'Training datasets', 'Validation datasets', 'Normal tissues', 
                                                    'Random PCs', 'Permuted datasets' ) )
pdf( paste0( dir_output, 'figure_v3.pdf' ), width=13, height=4 )
print( ggplot( df_toplot, aes( x=type, y=max_correlation ) ) + 
         geom_boxplot( ) +
         geom_hline( yintercept=0.5, color='red', linetype=2 ) +
         facet_grid( .~score_name ) +
         theme_bw( ) )
dev.off()

# discrete score validation
{
# data(classifier, package="DeSousa2013")
# data( diffGenes.f, package='DeSousa2013' )
# R <- 30
# mat_exprsTest <- exprs( eSets[['GSE13294_eset']] )
# gene_common <- intersect( rownames( mat_exprsTest ), rownames( mat_avgLoading ) )
# mat_exprsTest <- mat_exprsTest[gene_common, ]
# n <- ncol( mat_exprsTest )
# swTest <- sapply( 1:R, function( r ) {
#   mat_order <- t( sapply( 1:nrow( mat_exprsTest ), function( i ) return( sample.int( n, n ) ) ) )
#   mat_exprsTmp <- t( sapply( 1:nrow( mat_exprsTest ), function( i ) return( mat_exprsTest[i, ][mat_order[i, ]] ) ) )
#   rownames( mat_exprsTmp ) <- gene_common
#   
#   mat_exprsTmp <- mat_exprsTmp - apply( mat_exprsTmp, 1, median, na.rm=T )
#   
#   newExprs <- matrix( 0, length( setdiff( diffGenes.f, rownames( mat_exprsTmp ) ) ), ncol( mat_exprsTmp ) )
#   dimnames( newExprs ) <- list( setdiff( diffGenes.f, rownames( mat_exprsTmp ) ), colnames( mat_exprsTmp ) )
#   mat_exprsTmp <- rbind( mat_exprsTmp, newExprs )
#   mat_exprsTmp <- mat_exprsTmp[diffGenes.f, ]
#   
#   # Posterior classifying probabilities
#   clst_tmp <- as.numeric( pamr.predict( pam.rslt, mat_exprsTmp, thresh, type='class' ) )
#   dist_tmp <- dist( t( mat_exprsTmp[intersect( signature, rownames( mat_exprsTmp ) ), ] ), method='eucl' )
#   sw_tmp <- silhouette( clst_tmp, dist=dist_tmp )[, 'sil_width']
#   
#   return( mean( sw_tmp ) )
# } )
# 
# swAll <- sapply( eSets, function( eset ) {
#   mat_exprsTmp <- exprs( eset ) 
#   exprs.toplot <- mat_exprsTmp[intersect( signature, rownames( mat_exprsTmp ) ), ]
#   exprs.toplot <- exprs.toplot - apply( exprs.toplot, 1, median, na.rm=T )
#   dist_tmp <- dist( t( exprs.toplot ), method="euclidean" )
#   sw_tmp <- silhouette( as.numeric( pData( eset )$class ), dist=dist_tmp )[, 'sil_width']
#   return( mean( sw_tmp ) )
# } )
# 
# data( silh, package='DeSousa2013' )
# exprs.toplot <- sdat.f[signature, ]
# exprs.toplot <- exprs.toplot - apply( exprs.toplot, 1, median, na.rm=T )
# dist_tmp <- dist( t( exprs.toplot ), method="euclidean" )
# clus.f <- clus.f[ colnames( sdat.f ) ]
# swAMC <- mean( silhouette( clus.f, dist=dist_tmp )[, 'sil_width'] )
# 
# df_toplot <- rbind( data.frame( silhouette_width=swTest, type='Randomized datasets' ),
#                     data.frame( silhouette_width=swAll, type='All datasets' ),
#                     data.frame( silhouette_width=swAMC, type='AMC' ) )
# df_toplot$type <- factor( df_toplot$type, levels=c( 'AMC', 'All datasets', 'Randomized datasets' ) )
# pdf( paste0( dir_output, 'validation_discrete.pdf' ), width=7, height=4 )
# print( ggplot( df_toplot, aes( x=type, y=silhouette_width ) ) + 
#          geom_boxplot( ) +
#          theme_bw( ) )
# dev.off()
}
