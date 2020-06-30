rm( list=ls() )

library( Biobase )
library( GEOquery )
library( gplots )
library( labdsv )
library( ggplot2 )
library( tidyr )
library( grid )

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
colors <- gg_color_hue( 5 )

exprs <- read.table( 'data/breast/therese_2003.txt', sep='\t', header=F, row.names=1, skip=3 )
exprs <- as.matrix( exprs[ , -(1:3)] )
samples <- read.table( 'data/breast/therese_2003.txt', sep='\t', header=F, row.names=1, nrow=1 )
samples <- as.vector( unlist( samples[1, -(1:3)] ) )
colnames( exprs ) <- samples

genes <- read.table( 'data/breast/therese_2003.txt', sep='\t', header=F, skip=3 )[, 1:2]
# write.table( genes[, -1], file=paste0(output.dir, 'gene_clone_id.txt' ), row.names=F, col.names=F, quote=F )
genes.intrin <- read.table( 'data/breast/gene_symbol.txt', header=T, sep='\t', stringsAsFactors=F )[, 2]

class <- read.csv( 'data/breast/therese_2003_class.csv', header=T, stringsAsFactors=F )
class <- gather( class, key=class, value=sample )
class <- class[class$sample!='', ]
samples.common <- intersect( class$sample, samples )
rownames(class) <- class$sample
clus.f <- as.vector( class[samples.common, -2] )
clus.f[clus.f == "Luminal.A...0.32."] <- 'Luminal A'
clus.f[clus.f == "Luminal.B...0.28."] <- 'Luminal B'
clus.f[clus.f == "ERRB2....0.34."] <- 'ERRB2'
clus.f[clus.f == "Basal...0.41."] <- 'Basal'
clus.f[clus.f == "Normal.like...0.31."] <- 'Normal like'
clus.f <- factor(clus.f)
(silhouette(as.numeric(clus.f), dist(t(exprs.toplot))) %>% summary)$avg.width
exprs.toplot <- exprs[, samples.common]
cood.pcoa <- prcomp( dist( t( exprs.toplot ), method="euclidean" ) )$x
cood.toplot <- data.frame( PC1=cood.pcoa[, 1], PC2=cood.pcoa[, 2], 
                           PC3=cood.pcoa[, 3], PC4=cood.pcoa[, 4],
                           class=clus.f )


p1 <- ggplot( cood.toplot, aes(PC1, PC2) ) + 
  labs(title='Sorlie et al., 2003') + 
  aes(shape = factor(class)) + 
  geom_point(aes(colour = factor(class)) ) + 
  theme_bw() + 
  scale_colour_manual(name  ="Subtype", values=colors ) + 
  scale_shape_discrete( name  ="Subtype" ) + 
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        legend.direction = 'horizontal',
        legend.position = c(0, 0),
        legend.justification=c(0,0),
        legend.background = element_rect( colour='black' )) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
p2 <- ggplot( cood.toplot, 
              aes(PC3, PC4) ) + 
  labs(title='') + 
  aes(shape = factor(class)) + 
  geom_point(aes(colour = factor(class), alpha = class %in% c('Normal like', 'Luminal B', 'ERRB2') ) ) + 
  theme_bw() + 
  scale_colour_manual(name  ="Subtype", values=colors ) + 
  scale_shape_discrete( name  ="Subtype" ) + 
  scale_alpha_manual(values = c('TRUE' = 1, 'FALSE' = 0)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  guides(color = FALSE, shape = FALSE, alpha = FALSE)

library(cowplot)
pdf( 'CRC Manuscript/GenomeBiology/Resubmission/SuppMaterials/SuppFigure 3/SuppFigure3.pdf', width=10, height=5 )
plot_grid(p1, p2, ncol = 2)
dev.off()