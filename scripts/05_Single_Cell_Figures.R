# Single cell analysis - Part II: figures and analysis

#----------------------------------------#
# ------------ configuration ------------#
#----------------------------------------#

library(data.table)
library(dplyr)
library(Matrix)
library(DropletUtils)
library(scran)
library(scater)
library(vegan)
library(PCAtools)
library(cowplot)
library(ggplot2)
library(lme4)
library(ggtext)
library(grid)
library(gridExtra)

# set wd
wd <- '~/Documents/R/onecarbon/'
setwd(wd)

theme_set(theme_cowplot())

# read matched data
sce_match <- readRDS(file='data/Ascites_Tumor_Cells_Matched_for_Metabo_with_Cluster_Feb2020_clean.rds')


#----------------------------------------#
# --------- cluster and NNMT plot -------#
#----------------------------------------#
sce_match <- readRDS(file='data/Ascites_Tumor_Cells_Matched_for_Metabo_with_Cluster_Feb2020.rds')
dset <- logcounts(sce_match)
dset <- as.matrix(t(dset[c('AOX1', 'COL1A1', 'PTPRC', 'MYCN', 'GNLY', 'KLRB1','NNMT', 'EPCAM', 'CD8A', 'CD4', 'NKG7', 'CD3D', 'GZMA', 'RPS12'), ]))
dset <- data.frame(dset, 
                   clust=sce_match$cl_map, 
                   sce_match@int_colData@listData$reducedDims@listData$TSNE, 
                   comp=ifelse(grepl('_a', sce_match$tsamp), 'ascites', 'tumor'))

# tabulate NNMT expression
tab <- table(sign(dset$NNMT), sce_match$cl_map)
prop.tab <- data.frame(t(apply(tab, 2, function(x) x/sum(x))*100))
prop.tab$clust <- as.factor(rownames(prop.tab))

prop.tab <- prop.tab[(order(prop.tab$'X1')), ]
prop.tab$clust <- factor(prop.tab$clust, levels = prop.tab$clust)

meds <- setDT(dset[, c('clust', 'X1','X2')]) %>% .[ ,.(medX=median(X1), medY=median(X2)), by='clust']

# plot 
ppal <- rainbow(9, s=0.7, v=0.7, alpha = 0.4)
names(ppal) <- levels(prop.tab$clust)

g0 <- ggplot(dset, aes(x=X1, y=X2, color=CD8A)) + geom_point(size=2) +
  xlab('t-SNE 1') + ylab('t-SNE 2') +
  scale_color_manual('', values=c(rgb(0.7,0.7,0.7,0.3), rgb(1,0,0,0.3))) +
  guides(color=FALSE) 
g1 <- ggplot(dset, aes(x=X1, y=X2, color=clust)) + geom_point(size=2) +
  xlab('t-SNE 1') + ylab('t-SNE 2')  +
  scale_color_manual(values=ppal[prop.tab$clust]) + guides(color=FALSE) +
  annotate('text', x=meds$medX, y = meds$medY, label=meds$clust)
g2 <- ggplot(prop.tab, aes(x=clust, y=X1, fill=clust)) + 
  geom_bar(stat='identity') + 
  coord_flip() + xlab('') +
  ylab(expression(paste('% expressing ', italic(NNMT)))) +
  scale_fill_manual('', values=rainbow(9, s=0.7, v=0.7, alpha = 0.8)) +
  guides(fill=FALSE)
plot_grid(g0, g1, g2, rel_widths = c(1.6, 1.6, 1.3), ncol = 3)
dev.copy2pdf(file='figures/TSNE_clusters_w_NNMT_revised.pdf', width=12, height=4)

#----------------------------------------#
# ------------- marker plots ------------#
#----------------------------------------#

g1 <- ggplot(dset, aes(x=X1, y=X2, color=CD3D)) + geom_point() + 
  scale_color_gradient("*CD3D*", low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) +
  xlab('') + ylab('') +
  theme(legend.title = element_markdown())

g2 <- ggplot(dset, aes(x=X1, y=X2, color=EPCAM)) + geom_point() + 
  scale_color_gradient("*EPCAM*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) +
  xlab('') + ylab('') +
  theme(legend.title = element_markdown())

g3 <- ggplot(dset, aes(x=X1, y=X2, color=COL1A1)) + geom_point() + 
  scale_color_gradient("*COL1A1*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) +
  xlab('') + ylab('') +
  theme(legend.title = element_markdown())
 

g4 <- ggplot(dset, aes(x=X1, y=X2, color=PTPRC)) + geom_point() + 
  scale_color_gradient("*PTPRC*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

g5 <- ggplot(dset, aes(x=X1, y=X2, color=NNMT)) + geom_point() + 
  scale_color_gradient("*NNMT*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

g6 <- ggplot(dset, aes(x=X1, y=X2, color=AOX1)) + geom_point() + 
  scale_color_gradient("*AOX1*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

g7 <- ggplot(dset, aes(x=X1, y=X2, color=CD8A)) + geom_point() + 
  scale_color_gradient("*CD8A*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

g8 <- ggplot(dset, aes(x=X1, y=X2, color=RPS12)) + geom_point() + 
  scale_color_gradient("*RPS12*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

g9 <- ggplot(dset, aes(x=X1, y=X2, color=GZMA)) + geom_point() + 
  scale_color_gradient("*GZMA*",low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

p_all <- plot_grid(g5, g1, g2, g3, g4, g6, g7, g9, g8, ncol = 3)

y_grob <- textGrob('t-SNE2', 
                   gp=gpar(col="black", fontsize=15), rot=90)
x_grob <- textGrob('t-SNE1', 
                   gp=gpar(col="black", fontsize=15))
grid.arrange(arrangeGrob(p_all, left = y_grob, bottom=x_grob)) 
#dev.copy2pdf(file='figures/single_cell_marker_tSNEs_w_CD8.pdf', height=8, width=10)

#-----------------------------------------#
# ----------- TNF expression ------------ #
#-----------------------------------------#

tcell <- sce_match[ , sce_match$cl_map=='T cell' ]
tcell <- sce_match[ , sce_match$cl_map=='T cell' & sce_match$cl_cluster=='1']
mtcell <- findMarkers(tcell, tcell$Sample)


pdf <- data.frame(TNF=as.numeric(counts(tcell)['TNF', ]),
                  logTNF=as.numeric(logcounts(tcell)['TNF', ]),
                  pat=tcell$pat,
                  samp=ifelse(grepl('tumor', tcell$Sample), 'tumor', 'ascites')
)

# random sample for plot - to equalize sample size
# set.seed
set.seed(111)
spdf <- setDT(pdf)[, .SD[sample(.N, 300, replace = FALSE)], by = samp]

# get pvals with both find amrkers and neg bin glmer
pvals <- findMarkers(tcell, tcell$Sample)
tmod <- glmer.nb(TNF ~ samp + (1|pat), data=pdf)
pdf <- data.frame(group1='ascites', group2='tumor', p='**', y.position=5.5)


ggplot(spdf, aes(x=samp, y=logTNF)) + 
  geom_boxplot(outlier.shape = NA, fatten=0, color='black') +
  geom_jitter(width=0.1, height = 0.1, size=2, aes(color=pat)) +
  scale_color_manual('', values = rainbow(3, s=0.4, v=0.7, alpha = 0.5)) + 
  guides(color=FALSE) +
  xlab('') + 
  ylab(expression(log[2]~italic(TNF)~UMI~counts)) +
  stat_pvalue_manual(pdf, tip.length = 0)
#dev.copy2pdf(file='figures/TNF_in_T_cells_scRNA_w_Signif.pdf', height=4, width=3)
