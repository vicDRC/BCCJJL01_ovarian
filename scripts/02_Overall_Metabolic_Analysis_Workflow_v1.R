# This script reproduces the main analyses and figures for HGSC metabolomic analyses.
# Note some functions require sourcing
# Update: Sept 2020 metaboanalyst update breaks code at L48: mset initializiation?

#----------------------------------------#
# ------------ configuration ------------#
#----------------------------------------#

library(MetaboAnalystR)
library(dplyr)
library(genefilter)
library(ggplot2)
library(cowplot)
library(vegan)
library(limma)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)
library(ggpubr)
library(genefilter)

# cowplot theme
theme_set(theme_cowplot())

# set working dir
wd <- '~/Documents/R/onecarbon/'
setwd(wd)

# functions:
source('scripts/00_Metabolic_Processing_Functions.R')

#----------------------------------------#
# ------------ preprocessing ------------#
#----------------------------------------#

# note NA values have been previously imputed as 1000
# this workflow is based on MetaboAnalyst's default online workflow
# copied and reproduced here with MetaboAnalystR

# initialize object and normalize
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, 
                      'data/20190707_Cell_Screening_1213_metabo_2_clean.csv', 
                      format='colu', 'disc')
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet) # set at 1000 - has no effect
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "SumNorm", "LogNorm", 
                      "AutoNorm", 
                      ratio=FALSE, 
                      ratioNum=20)

# get normalized data
ndat <- mSet$dataSet$norm

# extract design from rownames
sname <- rownames(ndat)
sname <- gsub('-', 'n', sname)
sname <- gsub('\\+', '', sname)
dcode <- strsplit(sname, '_')
design <- data.frame(patient=as.factor(sapply(dcode, '[', 1)),
                     compart=as.factor(sapply(dcode, '[', 2)),
                     cell=as.factor(sapply(dcode, '[', 3)),
                     row.names=rownames(ndat))

design$treat <- with(design, as.factor(paste(compart, cell, sep='_')))
design$merge <- paste(design$patient, design$compart, design$cell)


#----------------------------------------#
# --------------- Heatmap ---------------#
#----------------------------------------#

# get top quartile of variable metabolites
varmet <- varFilter(data.matrix(t(ndat)), var.cutoff = 0.75)

# make additional design df
design2 <- design
design2$treat <- dplyr::recode(design2$treat, 
                               'A_45n'='A CD45-', 'A_4'='A CD4+', 'A_8'='A CD8+',
                               'T_45n'='T CD45-', 'T_4'='T CD4+' ,'T_8'='T CD8+'
                               )

design2$treat <- factor(design2$treat, 
                        levels = c('A CD4+','A CD8+', 'A CD45-', 'T CD4+' ,'T CD8+', 'T CD45-') 
                        )

# set colors for patients
patcols <- rainbow(6, v=0.6, s=0.6, alpha = 0.6)
names(patcols) <- unique(design$patient)

# set colors for treatments
treatcols <- c('A CD4+'=rgb(.7,.7,.7,0.5), 'A CD8+'=rgb(.4,.4,.4,0.5), 'A CD45-'=rgb(0,0,0,0.5),
               'T CD4+'=rgb(1,.7,.7,0.5), 'T CD8+'=rgb(1,.4,.4,0.6), 'T CD45-'=rgb(1,0,0,0.7)
               )

# set shapes for treatments
treatshapes <- c(15, 16, 17, 15, 16, 17)
names(treatshapes) <- names(treatcols)
               
# complexheatmap annotations: set colors
column_ha = with(design2, HeatmapAnnotation('Cell Type' = treat, 
                                            'Donor' = patient, 
                                          annotation_legend_param = list(`Cell Type` = list(nrow = 3),
                                                                          `Donor` = list(nrow=3)
                                                                         ),
                                           col = list(`Cell Type` = treatcols, `Donor` = patcols) 
                                          ) 
                 )

# Render Heatmap
hm <- Heatmap(varmet, show_column_names = FALSE, 
        clustering_method_columns = 'ward.D2',
        clustering_method_rows = 'ward.D2', 
        top_annotation = column_ha,
        heatmap_legend_param = list(title = expression(normalized~log~abundance), 
                                    title_position = 'leftcenter-rot',
                                    just=c('left', 'bottom')))
draw(hm, annotation_legend_side='top')
#dev.copy2pdf(file='figures/1_heatmap_draft_May2020.pdf', width=8, height=6)

#---------------------------
# heatmaps for tumor samples
# subset total dataset

column_ha = with(design2[design2$compart=='T', ], 
                 HeatmapAnnotation('Cell Type'=factor(treat, levels=c('T CD4+' ,'T CD8+', 'T CD45-')), 
                                           'Donor'=patient, 
                                           annotation_legend_param = list(`Cell Type`=list(nrow = 1),
                                                                          `Donor`=list(nrow=1)),
                                           col=list(`Cell Type`=treatcols, `Donor`=patcols) ) 
)

hm <- Heatmap(varmet[ , design2$compart=='T'], show_column_names = FALSE, 
              clustering_method_columns = 'ward.D2',
              clustering_method_rows = 'ward.D2', 
              top_annotation = column_ha,
              heatmap_legend_param = list(title = expression(normalized~log~abundance), 
                                          title_position = 'leftcenter-rot',
                                          just=c('left', 'bottom')))
draw(hm, annotation_legend_side='top')
#dev.copy2pdf(file='figures/1_heatmap_Tumor_May2020.pdf', width=6, height=6)

#---------------------------
# heatmap for ascites samples
# subset total dataset

column_ha = with(design2[design2$compart=='A', ], 
                 HeatmapAnnotation('Cell Type'=factor(treat, levels=c('A CD4+' ,'A CD8+', 'A CD45-')), 
                                                                    'Donor'=patient, 
                                                                    annotation_legend_param = list(`Cell Type`=list(nrow = 1),
                                                                                                   `Donor`=list(nrow=1)),
                                                                    col=list(`Cell Type`=treatcols, `Donor`=patcols) ) 
)

hm <- Heatmap(varmet[, design$compart=='A'], show_column_names = FALSE, 
              clustering_method_columns = 'ward.D2',
              clustering_method_rows = 'ward.D2', 
              top_annotation = column_ha,
              heatmap_legend_param = list(title = expression(normalized~log~abundance), 
                                          title_position = 'leftcenter-rot',
                                          just=c('left', 'bottom'))
              )
draw(hm, annotation_legend_side='top')
#dev.copy2pdf(file='figures/1_heatmap_Ascites_May2020pdf', width=4, height=4)



#----------------------------------------#
# ---------------- PCA ------------------#
#----------------------------------------#

# run PCA
rr <- rda(ndat)

# PCA plot - unsupervised
p <- ordiplot(rr, type='none', 
              scaling = 2, xlim=c(-4.25,2), bty = 'n', 
              xlab = 'PC1 (33% variance)', ylab = 'PC2 (11% variance)', 
              las = 1)
ordispider(p, groups = design2$patient, ylim=c(-1,2),
           col = patcols, lwd = 1.5)
points(p, what = 'sites', 
       col = treatcols[as.character(design2$treat)], 
       pch = treatshapes[as.character(design2$treat)], 
       cex = 1.5)
legend(-4,3, levels(design2$treat), pch = treatshapes[levels(design2$treat)], 
       pt.cex = 1.2, col = treatcols[levels(design2$treat)], bty = 'n', 
       cex = 0.7, title.adj = 0.2, title = 'Cell Type')
legend(1.5,3, names(patcols), pt.cex = 1.2, col=patcols, 
       bty = 'n', lty = 1, lwd=2, xpd = TRUE, cex=0.7, title.adj = -.1, title = 'Donor')

#dev.copy2pdf(file='figures/PCA_not_controlling_for_donor_Feb2020.pdf', height=4.5, width=4.5)


# PCA plot - conditioned on patient
rr <- rda(ndat ~ Condition(design2$patient))

p <- ordiplot(rr, type = 'none', scaling = 2, bty = 'n', 
              xlab = 'PC1 (34% variance)', ylab = 'PC2 (10% variance)', 
              las = 1)
ordihull(p, groups = design2$treat, 
           col = treatcols[levels(design2$treat)], lwd = 1.5)
points(p, what = 'sites', 
       col = treatcols[as.character(design2$treat)], 
       pch = treatshapes[as.character(design2$treat)], 
       cex = 1.5)

#dev.copy2pdf(file='figures/PCA_conditioned_for_patient_Feb2020.pdf', height=4.5, width=4.5)


#----------------------------------------#
#- Differential Analysis with Limma -----#
#----------------------------------------#

# make design matrix
# use group means parameterization
# average reps plots
y <- new("EList", list(E = t(ndat), targets = design))
y2 <- avearrays(y, ID = design$merge)
treat <- factor(y2$targets$treat)
pat <- factor(y2$targets$patient)
compart <- factor(y2$targets$compart)
cell <- factor(y2$targets$cell)

# for group means paramterization
dmat_np <- model.matrix(~ 0 + treat)
colnames(dmat_np) <- gsub('treat', '', colnames(dmat_np))
corfit <- duplicateCorrelation(y2, dmat_np, block=pat)
fitb <- eBayes(lmFit(y2, dmat_np, block=pat, correlation=corfit$consensus))


cont <- makeContrasts('T_8-A_8', 
                      'T_4-A_4',
                      'T_45n-A_45n',
                      'T_8-T_4',
                      'T_8-T_45n',
                      'T_4-T_45n',
                      'A_8-A_4',
                      'A_8-A_45n',
                      'A_4-A_45n',
                      levels=dmat_np
)        

fitc <- contrasts.fit(fitb, contrasts = cont)
ebc <- eBayes(fitc, robust = TRUE)

decideTests(ebc, method='separate')
topTable(ebc, n=10, coef=1:9)

# compute overall F stats for P-vals
Fstats <- topTableF(ebc, number=nrow(ebc))

#---------------------------------------#
#--- inter-patient differences ---------#
#---------------------------------------#

# include patient as fixed effect to get metabolites that differentiate patients

dmat_p <- model.matrix(~ 0 + pat + treat)
colnames(dmat_p) <- gsub('treat', '', colnames(dmat_p))
fit <- eBayes(lmFit(y2, dmat_p))
tt <- topTable(fit, coef=1:6)
topM <- rownames(tt)[1:4] # get top 4 metabolites

plot_df <- data.frame(t(y2$E), y2$targets, check.names = FALSE)

# make plots
g1 <- ggplot(plot_df, aes(x=patient, y=Dodecanoylcarnitine)) + 
  geom_boxplot() + theme_cowplot() + xlab('')
g2 <- ggplot(plot_df, aes(x=patient, y=`Hydroxyisocaproic acid`)) +
  geom_boxplot() + theme_cowplot() + xlab('')
g3 <- ggplot(plot_df, aes(x=patient, y=Putrescine)) +
  geom_boxplot() + theme_cowplot() + xlab('')
g4 <- ggplot(plot_df, aes(x=patient, y=`L-3-Phenyllactic acid`)) +
  geom_boxplot() + theme_cowplot() + xlab('')

p_group <- plot_grid(g1, g2, g3, g4)
p_group <- globLabs(p_group, xlab=expression(patient), ylab='')
#dev.copy2pdf(file='figures/Patient_Boxplots_April2020.pdf', height=6, width=8)


# Make heatmap of top tumor diff metabolites
# tumor top metabo HM

tumtop <- rownames(topTable(ebc, coef=4:6, number = 20))

column_ha = with(design2[design2$compart=='T', ], HeatmapAnnotation('Cell Type' = factor(treat), 
                                                                    'Patient' = patient, 
                                                                    annotation_legend_param = list(`Cell Type`=list(nrow = 1),
                                                                                                   `Patient`=list(nrow=1)),
                                                                    col=list(`Cell Type`=treatcols, `Patient`=patcols) ) 
)

hm <- Heatmap(t(ndat[design2$compart=='T', tumtop]), show_column_names = FALSE, 
              clustering_method_columns = 'ward.D2',
              clustering_method_rows = 'ward.D2', 
              top_annotation = column_ha,
              heatmap_legend_param = list(title = expression(normalized~log~abundance), 
                                          title_position = 'leftcenter-rot',
                                          just=c('left', 'bottom')))
draw(hm, annotation_legend_side='top')
#dev.copy2pdf(file='figures/Heatmap_Tumor_sigs_May2020.pdf', width=5, height=6)

#--------------------------------------#

# ascites #
asctop <- rownames(topTable(ebc, coef=7:9, number = 20))

column_ha = with(design2[design2$compart=='A', ], HeatmapAnnotation('Cell Type'=factor(treat), 
                                                                    'Patient'=patient, 
                                                                    annotation_legend_param = list(`Cell Type`=list(nrow = 1),
                                                                                                   `Patient`=list(nrow=1)),
                                                                    col=list(`Cell Type`=treatcols, `Patient`=patcols) ) 
)

hm <- Heatmap(t(ndat[design2$compart=='A', asctop]), show_column_names = FALSE, 
              clustering_method_columns = 'ward.D2',
              clustering_method_rows = 'ward.D2', 
              top_annotation = column_ha,
              heatmap_legend_param = list(title = expression(normalized~log~abundacne), 
                                          title_position = 'leftcenter-rot',
                                          just=c('left', 'bottom')))
draw(hm, annotation_legend_side='top')
#dev.copy2pdf(file='figures/Heatmap_Ascites_sigs_May2020.pdf', width=6, height=6)


#----------------------------------------#
#----- Metabolite Specific Plots --------#
#----------------------------------------#


# quick convenience plotter
plotMetabo <- function(metabo, data=y2, ebayes=eb, ...) {
  mnc <- t(data.frame(data$E))[, metabo]
  design <- data$targets
  design$treat <- factor(design$treat, levels=c("A_4", "A_8", "A_45n", "T_4","T_8", "T_45n"))
  
  sigs <- getAllcontrasts(ebayes, metabo, trim=TRUE, ...)
  
  ggplot(design, aes(y=mnc, x=treat, color=compart))  +
    geom_boxplot() + 
    scale_color_manual('', values=c(rgb(0.1,0.1,0.1,0.8), rgb(1,0,0,0.8)), labels=c('Ascites', 'Tumor')) + 
    scale_x_discrete( labels=c('CD4+', 'CD8+', 'CD45-','CD4+', 'CD8+', 'CD45-')) +
    guides() +
    xlab('') + theme(axis.text.x = element_text(face="bold", angle=45, hjust=1)) +
    ylab(bquote(.(metabo)~(normalized) )) +
    stat_pvalue_manual(data = sigs, label = "symbs", tip.length = 0, size=3, hide.ns = TRUE)
}

# render plots
plotMetabo('1-Methylnicotinamide', ymax=3, incr=0.2, data = y2, ebayes=ebc)
#dev.copy2pdf(file='figures/MNA_boxplot_May2020_new_analysis.pdf', width=4, height=4)

plotMetabo('Adenosine', ymax=3.75, incr=0.2,  ebayes=ebc)
#dev.copy2pdf(file='figures/Adenosine_boxplot_May2020_new_analysis.pdf', width=4, height=4)

plotMetabo('L-Kynurenine', ymax=3.5, incr=0.2, ebayes=ebc)
#dev.copy2pdf(file='figures/L-kynurenine_boxplot_May2020_new_analysis.pdf', width=4, height=4)





