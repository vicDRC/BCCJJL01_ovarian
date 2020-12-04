# Addressing reviewer comments 
# Sept 2020

library(data.table)
library(dplyr)
library(Matrix)
library(DropletUtils)
library(scran)
library(biomaRt)
library(scater)
library(vegan)
library(PCAtools)
library(cowplot)
library(ggplot2)
library(ggtext)
library(grid)
library(gridExtra)


# get data from Izar et al. 2020 - download as GSE
# make metadata and expression objects

ov_sce_meta <- t(head(fread('~/Documents/R/onecarbon/data/GSE146026_Izar_HGSOC_ascites_10x_log.tsv'), 7))
colnames(ov_sce_meta) <- ov_sce_meta[1, ]
ov_sce_meta <- data.frame(ov_sce_meta[-1, ])

ov_sce <- fread('~/Documents/R/onecarbon/data/GSE146026_Izar_HGSOC_ascites_10x_log.tsv', skip=8)
ov_sce_expr <- Matrix(data.matrix(ov_sce)[, -1])
rownames(ov_sce_expr) <- ov_sce$V1
colnames(ov_sce_expr) <- ov_sce_meta$X10x_barcode

# extract relevant markers, as well as tSNE embeddings
dset <- t(as.matrix(ov_sce_expr[c('GZMB', 'AOX1','COL1A1', 'PTPRC', 'GNLY', 'KLRB1', 'NNMT', 'EPCAM', 'CD8A', 'CD4', 'NKG7', 'CD3D'), ]))
dset <- data.table(X1=ov_sce_meta$TSNE_x, X2=ov_sce_meta$TSNE_y,  dset)
dset <- dset[ , lapply(.SD, as.numeric)]


# make supp figure
theme_set(theme_cowplot())

g1 <- ggplot(dset, aes(x=X1, y=X2, color=CD3D)) + geom_point() + 
  scale_color_gradient("*CD3D*", 
                       low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) +
  xlab('') + ylab('') +
  theme(legend.title = element_markdown())

g2 <- ggplot(dset, aes(x=X1, y=X2, color=EPCAM)) + geom_point() + 
  scale_color_gradient("*EPCAM*",
                       low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) +
  xlab('') + ylab('') +
  theme(legend.title = element_markdown())

g3 <- ggplot(dset, aes(x=X1, y=X2, color=COL1A1)) + geom_point() + 
  scale_color_gradient("*COL1A1*",
                       low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) +
  xlab('') + ylab('') +
  theme(legend.title = element_markdown())

g4 <- ggplot(dset, aes(x=X1, y=X2, color=PTPRC)) + geom_point() + 
  scale_color_gradient("*PTPRC*",
                       low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

g5 <- ggplot(dset, aes(x=X1, y=X2, color=NNMT)) + geom_point() + 
  scale_color_gradient("*NNMT*",
                       low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

g6 <- ggplot(dset, aes(x=X1, y=X2, color=AOX1)) + geom_point() + 
  scale_color_gradient("*AOX1*",
                       low=rgb(0.7,0.7,0.7,0.1),
                       high=rgb(0,0,0.7,0.8)) + 
  theme(legend.title = element_markdown()) +
  xlab('') + ylab('')

p_all <- plot_grid(g5, g1, g2, g3, g4, g6)

y_grob <- textGrob('t-SNE2', 
                   gp=gpar(col="black", fontsize=15), rot=90)
x_grob <- textGrob('t-SNE1', 
                   gp=gpar(col="black", fontsize=15))

grid.arrange(arrangeGrob(p_all, 
                         left = y_grob, 
                         bottom = x_grob))

dev.copy2pdf(file='~/Documents/R/onecarbon/figures/Supp_Izar_et_al_NNMT_expression.pdf', width=10, height=5)
