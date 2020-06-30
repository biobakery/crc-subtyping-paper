library( Biobase )
library( ggplot2 )
library( grid )
data.dir <- 'data/clean_new/'
esets.dir <- 'data/clean_new/'

load( paste0( data.dir, 'setNames.RData' ) )
load( paste0( data.dir, 'trainingSetNames.RData' ) )
eSets <- list()
for ( set in setNames ) {
  load( paste( esets.dir, set,".RData",sep="" ) )
  eSets[[set]] <- get( set )
}

loadings <- read.csv( 'results/PCA/loadings/avg_loadings.csv', header=T, row.names=1 )
ns <- seq( 100, 3000, by=100 )
ps <- list(  )
for ( i_score in 1:2 ) {
  score <- colnames( loadings )[i_score]
  results <- sapply(ns, function(n) {
    genes.sig <- rownames( loadings )[order( abs( loadings[, score] ), decreasing=T )[1:n]]
    apply( sapply( setNames, function( set ){ 
      genes.sig.tmp <- intersect( genes.sig, featureNames( eSets[[set]] ) )
      scores.sig.tmp <- t( exprs( eSets[[set]] )[ genes.sig.tmp, ] ) %*% loadings[genes.sig.tmp, score]
      c(cor( scores.sig.tmp, 
             pData( eSets[[set]] )[, score], 
             method='pearson', use='pairwise.complete.obs' ),
        cor.test( scores.sig.tmp, 
                  pData( eSets[[set]] )[, score], 
                  method='pearson', use='pairwise.complete.obs' )$conf.int)
    } ), 1, mean, na.rm=T )
  }) %>% t %>% data.frame(n = ns)
  ps[[score]] <- ggplot( data=results, aes( n, X1 ) ) + geom_point( size=1.5) + 
    geom_line() + theme_bw() + geom_errorbar(aes(ymin = X2, ymax = X3)) +
    labs(title=paste0( 'PCSS', i_score ), x='Number of Genes', y='Pearson Correlation') + 
    scale_y_continuous(breaks=seq(c(0.875, 0.825)[i_score], 1, by=0.025), limits=c(c(0.875, 0.825)[i_score], 1)) + 
    scale_x_continuous(breaks=seq(500, 3000, by=500))
}

pdf( 'CRC Manuscript/GenomeBiology/Resubmission/SuppMaterials/SuppFigure 6/figure.pdf', width=10, height=5 )
print ( ps[[1]], vp=viewport( x=0.25, y=0.5, width=0.5, height=1 ) )
print ( ps[[2]], vp=viewport( x=0.75, y=0.5, width=0.5, height=1 ) )
dev.off( )