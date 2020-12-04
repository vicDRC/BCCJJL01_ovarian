# basic functions for metabolomics analysis
# Sept 2019
# Revised March 2020

# Note not all functions used in final analysis

# run MetaboAnalyst ORA
pathMetab <- function(metabs, 
                      set.backgrnd=TRUE, 
                      backgrnd="all_kegg.csv") {
  oSet <- InitDataObjects("conc", "pathora", FALSE)
  oSet<-Setup.MapData(oSet, metabs)
  oSet<-CrossReferencing(oSet, "name")
  oSet<-CreateMappingResultTable(oSet)
  if(set.backgrnd) oSet<-Setup.KEGGReferenceMetabolome(oSet, backgrnd)
  oSet<-SetKEGG.PathLib(oSet, "hsa")
  oSet<-SetMetabolomeFilter(oSet, T)
  oSet<-CalculateOraScore(oSet, "rbc", "hyperg")
  oSet
}

# run MetaboAnalyst ORA
enrichMetab <- function(metabs, 
                        set.backgrnd=TRUE,
                        backgrnd="all_cmpnd.csv") {
  require(MetaboAnalystR)
  
  oSet<-InitDataObjects("conc", "msetora", FALSE)
  oSet<-Setup.MapData(oSet, metabs);
  oSet<-CrossReferencing(oSet, "name")
  oSet<-CreateMappingResultTable(oSet)
  if(set.backgrnd) oSet<-Setup.HMDBReferenceMetabolome(oSet, backgrnd)
  oSet<-SetMetabolomeFilter(oSet, T)
  oSet<-SetCurrentMsetLib(oSet, "smpdb_pathway", 2)
  oSet<-CalculateHyperScore(oSet)
  oSet<-CalculateHyperScore(oSet)
  oSet

}

# make expression volcano plot - requires decide test from limma
plotVolcano <- function(coefnum, 
                        title,
                        res,
                        pval.thresh=0.005, 
                        fc.thresh=1, 
                        xlab=expression(log[10]~fold~change), 
                        ylab=expression(-log[10]~P[adj])) {
  
  tf <- topTable(eb, coef = coefnum, number = nrow(eb$coefficients))
  
  maplabs <- ifelse(abs(tf$logFC) > fc.thresh & tf$adj.P.Val < pval.thresh, rownames(tf), '' )
  
  ggplot(tf, aes(x=logFC, y=-log10(adj.P.Val), 
                 color=ifelse(res[rownames(tf), coefnum]==1 | res[rownames(tf), coefnum]==-1, 'blue', 'black'))) + 
    geom_point() + 
    scale_color_manual('', values=c('black', 'blue')) + guides(color=FALSE) +
    xlab(xlab) +
    ylab(ylab) +
    xlim(c(-3,3))  + ggtitle(title) +
    geom_text_repel(label=maplabs, color='dark gray', force=1, size=4)
    
}

# add global labels to cowplot plot_grid object
# borrowed grob code from stack overflow
globLabs <- function(plot, 
                     ylab=expression(-log[10]~P[adj]), 
                     xlab=expression(log[10]~fold~change), 
                     main='') {
  
  require(gridExtra)
  
  y.grob <- textGrob(ylab, 
                     gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
  x.grob <- textGrob(xlab, 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
  top.grob <- textGrob(main, 
                     gp=gpar(fontface="bold", col="black", fontsize=15))

  grid.arrange(arrangeGrob(plot, left = y.grob, bottom=x.grob, top=top.grob))

}


boxplotMod <- function(mod, plotdf) {
  mod <- sym(mod)
  ggplot(plotdf, aes(x=treat, y=!!mod)) + 
    scale_fill_manual('', values = c(rgb(0,0,0,0.5), rgb(1,0,0,0.5)), labels=c('ascites', 'tumor')) +
    geom_boxplot(aes(fill=compart)) +
    xlab('') + ylab('module eigenvector') +
    scale_x_discrete(labels=c('CD4+', 'CD8+', 'CD45-','CD4+', 'CD8+', 'CD45-')) +
    theme(axis.text.x = element_text(face="bold", angle=45, hjust=1))
}

plotMOI <- function(mod, modLabs, modEigs, dset, toptab, title=mod, p.thresh=3, cor.thresh=0.5) {
  require(ggplot2)
  require(ggrepel)
  
  mset <- getfromMod(module=mod, mLabels = modLabs, dset)
  mcors <- sapply(mset, function(x) cor(dset[, x], modEigs[ , paste0('ME', mod)]) )
  pvals <- -log10(toptab[mset, 'adj.P.Val'])
  
  df <- data.frame(mset, mcors, pvals, stringsAsFactors = FALSE)

  df$labs <- df$mset
  df$labs[df$pvals < p.thresh | abs(df$mcors) < cor.thresh] <- ''
  #df$labs <- gsub(' ', '\n', df$labs)
  
  
  ggplot(df, aes(x=mcors, y=pvals)) +
    geom_point() + 
    geom_text_repel(label=df$labs, color='dark gray', size=3) +
    ylab(expression(-log[10]~p[adj]))+
    xlab('correlation with eigenvector') +
    geom_vline(xintercept=0, size=0.1, color='gray') +
    geom_hline(yintercept=1.3, size=0.1, col='red') +
    ggtitle(title) + xlim(c(-1,1))
}

getfromMod <- function(module, mLabels, dset=dataset) {
  mlab <- data.frame(module=as.factor(mLabels), 
                     row.names=colnames(dset))
  mods <- rownames(mlab[mlab$module==module, ,drop=FALSE])
  mods
}

getAllcontrasts <- function(emp.bayes, metab, trim=TRUE, ymax, incr=0.15) {
  
  df <- data.frame(group1=gsub('-.*$', '', colnames(emp.bayes$p.value)),
                   group2=gsub('^.*-', '', colnames(emp.bayes$p.value)),
                   p=p.adjust(as.numeric(emp.bayes$p.value[metab, ]), method='BH'),
                   y.position=ymax )
  df$symbs <- ifelse(df$p < 0.001, '***', ifelse(df$p < 0.01, '**', ifelse(df$p < 0.05, '*', 'n.s.')))
  
  if(trim) { df <- df[df$p <= 0.05, ] }
  
  df$y.position <- df$y.position - seq(0, incr*nrow(df)-incr, by=incr)
  
  df
                   
}

