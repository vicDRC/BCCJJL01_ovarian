# analysis of flow data using FAUST
# ran two projects - one with and one without metabolic markers 
# to link events together without clustering on metabolic markers
# Nov 2019

# based on FAUST vignette
# https://github.com/RGLab/FAUST

#----------------------------------------#
# ------------ configuration ------------#
#----------------------------------------#

#library(flowWorkspaceData)
library(flowWorkspace)
library(ggdendro)
library(scamp)
library(ggplot2)
library(cowplot)
library(knitr)
library(dplyr)
library(tidyr)
library(faust)
library(CytoML)
library(openCyto)
library(CytoML)
library(pheatmap)
library(data.table)
library(ComplexHeatmap)

# cowplot theme
theme_set(theme_cowplot())

# set working dir
wd <- '~/Documents/R/onecarbon/'
setwd(wd)

#-----------------------------------------------------------------#
#------------------------ begin analysis -------------------------#
#-----------------------------------------------------------------#

# get flowjo workspace
# manually gated 
ws <- open_flowjo_xml('data/flow_panel/25Nov2019_MK_flow-1.wsp')
gt <- flowjo_to_gatingset(ws) # interactive

gh_get_transformations(gt[[1]])

# gating strategy - starts at live cells
startingNode <- 'Live cells'

# channels to exclude for initial analysis - note viability typo
activeChannelsIn <- markernames(gt)[!(markernames(gt) %in% c('Viaibility', '2-NBDG', 'MT DR'))]
activeChannelsIn2 <- markernames(gt)[!(markernames(gt) %in% c('Viaibility'))]

# define project directories -- make an additional project for all markers
projPath <- file.path('FAUST')
projPath2 <- file.path('FAUSTall') # run second time for metabolic markers

# make dirs
dir.create(projPath, recursive=TRUE)
dir.create(projPath2, recursive=TRUE)

# run initial FAUST: allow automatic bounding
# default params
#faust(
#  gatingSet = gt,
#  experimentalUnit = "name",
#  activeChannels = activeChannelsIn,
#  startingCellPop = startingNode,
#  projectPath = projPath,
#  depthScoreThreshold = 0.05,
#  selectionQuantile = 1.0,
#  debugFlag = FALSE,
#  #set this to the number of threads you want to use on your system
#  threadNum = 2,
#  nameOccuranceNum=1,
#  seedValue = 101,
#  annotationsApproved = FALSE # set to false before we inspect the scores plots.
#) 

# inspect and re-run
faust(
  gatingSet = gt,
  experimentalUnit = "name",
  activeChannels = activeChannelsIn,
  startingCellPop = startingNode,
  projectPath = projPath,
  depthScoreThreshold = 0.05,
  selectionQuantile = 1.0,
  debugFlag = FALSE,
  #set this to the number of threads you want to use on your system
  threadNum = 2,
  nameOccuranceNum=1,
  seedValue = 101,
  annotationsApproved = TRUE # set to false before we inspect the scores plots.
) 

# now run with all markers to provide a coherent dataset for MTDR / 2NBDG analysis
# the annotations from initial runs can be linked to marker expression 
# from the outputs here.
faust(
  gatingSet = gt,
  experimentalUnit = "name",
  activeChannels = activeChannelsIn2,
  startingCellPop = startingNode,
  projectPath = projPath2,
  depthScoreThreshold = 0.05,
  selectionQuantile = 1.0,
  debugFlag = FALSE,
  #set this to the number of threads you want to use on your system
  threadNum = 2,
  nameOccuranceNum=1,
  seedValue = 101,
  annotationsApproved = TRUE # set to false before we inspect the scores plots.
) 

# read output matrix
countMatrix <- readRDS(file.path(projPath,"faustData","faustCountMatrix.rds"))

# write gating plots - takes awhile
for (col in colnames(countMatrix)) {
  for (r in rownames(countMatrix)) {
    faust:::plotFaustGates(col, r, projPath)
  }
}

#---------------------------------------------------------------#
#-------------- integrate with metabolic markers ---------------#
#---------------------------------------------------------------#

# make our own count matrix to allow manual merging of phenotypes if needed

# get sample names/paths
samples <- list.dirs('FAUST/faustData/sampleData/')[-1]

# quick function to read in data
mergeMetabo <- function(samp) {
  texp <- readRDS(paste0(samp, '/exprsMat.rds')) 
  tann <- fread(paste0(samp, '/faustAnnotation.csv'), header=FALSE)
  bexp <- readRDS(paste0(gsub('FAUST', 'FAUSTall', samp), '/exprsMat.rds') )
  data.table(Annotation=tann$V1, bexp, id=gsub('^.*/| .*$', '', samp))
}

# read in samples
all_phenos <- lapply(samples, mergeMetabo)
long_phenos <- rbindlist(all_phenos) %>% setDT


# substitute relevant cell type merges here:
# PD1 expression EPCAM+ cells looks problematic - collapse
long_phenos$Annotation <- gsub('CD3~1~2~AF700~1~2~CD45RO~1~2~CD45~1~2~CD8~1~2~EpCAM~2~2~PD1~2~2~CCR7~1~2~CD25~1~2~', 
                               'CD3~1~2~AF700~1~2~CD45RO~1~2~CD45~1~2~CD8~1~2~EpCAM~2~2~PD1~1~2~CCR7~1~2~CD25~1~2~', 
                                long_phenos$Annotation)
long_phenos$Annotation <- gsub('AF700', 
                               'CD4', 
                               long_phenos$Annotation)
long_phenos$Annotation <- gsub('~2~2~', '+', long_phenos$Annotation)
long_phenos$Annotation <- gsub('~1~2~', '-', long_phenos$Annotation)


# calculate mean metabolic scores per phenotype: lacks power but conservative
med_phenos <- setkey(long_phenos[, lapply(.SD, median), by=c('Annotation', 'id'), .SDcols=c('2-NBDG', 'MT DR')], Annotation)
count_phenos <- setkey(long_phenos[, .(count=.N), by=c('Annotation', 'id')], Annotation)

# cast data table
c_counts <- dcast(count_phenos, Annotation ~ id, value.var='count', fill=0)
cm <- t(data.frame(c_counts, row.names=1))
cm <- apply(cm, 1, function(x) (x/sum(x))*100)

# require a pop >2% on average to retain
cm <- cm[rowMeans(cm) > 2, ]

# remove un-annotated cells from consideration
pops <- rownames(cm)[!rownames(cm)=='0_0_0_0_0']
cm <- cm[pops, ]

pres_pops <- Reduce(intersect, lapply(all_phenos, function(x) x$Annotation))

# wrapper for paired t-test across phenotypes: not much significant, but conservative analysis
wrapT <- function(pop, long_pheno=med_phenos) {
  pop_ac <- long_pheno[long_pheno$Annotation==pop & grepl('AC', long_pheno$id), ]
  pop_tb <- long_pheno[long_pheno$Annotation==pop & grepl('TB', long_pheno$id), ]
  if(nrow(pop_ac)==6 & nrow(pop_tb)==6) {
    tnbdg <- t.test(pop_ac$`2-NBDG`, pop_tb$`2-NBDG`, paired=TRUE)
    tmtdr <- t.test(pop_ac$`MT DR`, pop_tb$`MT DR`, paired=TRUE)
    list(NBDG=tnbdg, MTDR=tmtdr)
  } else {
    message('incomplete observations')
  }
}

# run t tests
tstats <- lapply(pops, wrapT)
names(tstats) <- pops

# make population codes for plotting
popcode <- rep('other', length(pops))
popcode[grepl('CD45-', pops) &  grepl('EpCAM-', pops)] <-  'CD45-/EPCAM-'
popcode[grepl('CD45-', pops) &  grepl('EpCAM\\+', pops)] <-  'CD45-/EPCAM+'
popcode[grepl('CD3\\+', pops) & grepl('CD4\\+', pops) ] <-  'CD4+'
popcode[grepl('CD3\\+', pops) & grepl('CD8\\+', pops) ] <-  'CD8+'

popcols <- rainbow(5, s=0.6, v=0.6, alpha=0.75)

# pull ascites and tumor for MTDR 
ac <- med_phenos[grepl('AC', med_phenos$id), ]
ac <- ac[Annotation %in% pops]
ac_meds <- ac[, lapply(.SD, median), by=c('Annotation'), .SDcols=c('MT DR', '2-NBDG') ] %>% setkey(Annotation)
ac_meds <- ac_meds[pops]
ac_meds <- ac_meds[ , popcode:=popcode] %>% setkey(Annotation)
ac_ord <- ac_meds[order(`MT DR`) , Annotation]
ac$popcode <- ac_meds[as.character(ac$Annotation), 'popcode']
ac$V1 <- factor(ac$Annotation, levels=ac_ord)
names(popcols) <- unique(ac$popcode)

ggplot(ac, aes(group=V1, y=`MT DR`, x=V1, fill=popcode, color=popcode)) + 
  geom_boxplot() + 
  xlab('') + 
  ylab('MT DR (MFI)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) +
  scale_fill_manual('', values=popcols) + 
  scale_color_manual('', values=popcols) + 
  guides(color=FALSE)
dev.copy2pdf(file='figures/MT_DR_ascites_pops.pdf')

ggplot(ac, aes(group=V1, y=`2-NBDG`, x=V1, fill=popcode, color=popcode)) + 
  geom_boxplot(outlier.color=NULL) + 
  xlab('') +
  ylab('2-NBDG (MFI)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) + 
  scale_fill_manual('', values=popcols) + 
  scale_color_manual('', values=popcols) + 
  guides(color=FALSE)
#dev.copy2pdf(file='figures/NBDG_ascites_pops.pdf')
  
#------- repeat for tumor -----------#
tb <- med_phenos[grepl('TB', med_phenos$id)]
tb <- tb[Annotation %in% pops]
tb_meds <- tb[, lapply(.SD, median), by=c('Annotation'), .SDcols=c('MT DR', '2-NBDG') ] %>% setkey(Annotation)
tb_meds <- tb_meds[pops]
tb_meds <- tb_meds[ , popcode:=popcode] %>% setkey(Annotation)
tb_ord <- tb_meds[order(`MT DR`), Annotation]
tb$popcode <- tb_meds[as.character(tb$Annotation), 'popcode']
tb$V1 <- factor(tb$Annotation, levels=ac_ord)

ggplot(tb, aes(group=V1, y=`MT DR`, x=V1, fill=popcode, color=popcode)) + 
  geom_boxplot() + 
  xlab('') +
  ylab('MT DR (MFI)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) + 
  scale_fill_manual('', values=popcols) + 
  scale_color_manual('', values=popcols) + 
  guides(color=FALSE)
#dev.copy2pdf(file='figures/MT_DR_tumor_pops.pdf')

ggplot(tb, aes(group=V1, y=`2-NBDG`, x=V1, fill=popcode, color=popcode)) + 
  geom_boxplot() + 
  xlab('') +
  ylab('2-NBDG (MFI)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) + 
  scale_fill_manual('', values=popcols) + 
  scale_color_manual('', values=popcols) + 
  guides(color=FALSE)
#dev.copy2pdf(file='figures/NBDG_tumor_pops.pdf')

pdf <- data.frame(tb_meds, ac_meds)

g1 <- ggplot(pdf, aes(x=MT.DR, y=X2.NBDG, colour=popcode)) + 
  geom_point(size=4) + 
  ggtitle('tumor populations') + 
  xlab('MT DR (MFI)') + 
  ylab('2 NBDG (MFI)') + 
  scale_color_manual('', values=popcols) + 
  geom_smooth(method='lm', fill=NA, aes(colour=NULL))
g2 <- ggplot(pdf, aes(x=MT.DR.1, y=X2.NBDG.1, colour=popcode))  +
  geom_point(size=4) + 
  ggtitle('ascites populations') + 
  xlab('MT DR (MFI)') + 
  ylab('2 NBDG (MFI)') + 
  scale_color_manual('', values=popcols) +
  geom_smooth(method='lm', fill=NA, aes(colour=NULL))
plot_grid(g1, g2)
#dev.copy2pdf(file='figures/MTDR_vs_2NBDG_pop_comps.pdf', width=8, height=4)


g1 <- ggplot(pdf, aes(x=MT.DR, y=MT.DR.1, colour=popcode)) + 
  geom_point(size=4) + ggtitle('MT DR') + xlab('MT DR tumor (MFI)') + ylab('MT DR ascites (MFI)') +
  scale_color_manual('', values=popcols)+ geom_smooth(method='lm', fill=NA, aes(colour=NULL))
g2 <- ggplot(pdf, aes(x=X2.NBDG, y=X2.NBDG.1, colour=popcode)) + 
  geom_point(size=4) + ggtitle('2 NBDG') + xlab('2 NBDG tumor (MFI)') + ylab('2 NBDG ascites (MFI)') + 
  scale_color_manual('', values=popcols) + geom_smooth(method='lm', fill=NA, aes(colour=NULL))
plot_grid(g1, g2) 
#dev.copy2pdf(file='figures/MTDR_and_2NBDG_pops.pdf', width=8, height=4)

cor.test(pdf$MT.DR, pdf$MT.DR.1, method='sp')
cor.test(pdf$X2.NBDG, pdf$X2.NBDG.1, method='sp')

#---------------------------------------------------#
#----------------- complex heatmap -----------------#
#---------------------------------------------------#

popcols <- popcols[unique(ac$popcode)]

anno_col <- data.frame(sample=ifelse(grepl('AC', colnames(cm)), 'Ascites', 'Tumor'), row.names=colnames(cm))
anno_row <- data.frame(popcode, row.names=pops)

column_ha = with(anno_col, HeatmapAnnotation('Sample Type'=sample, 
                                            annotation_legend_param = list(`Sample Type`=list(nrow = 1)),
                                            col = list(`Sample Type`=c(Ascites=rgb(0.3,.3,.3, 0.75), Tumor=rgb(1, 0,0,0.75))) ) 
  )

row_ha = with(anno_row, rowAnnotation( 'Population'=popcode,
                                             annotation_legend_param = list(`Population`=list(ncol=1)),
                                             col = list(`Population`=popcols) ) 
  )

hm <- Heatmap(cm, show_column_names = FALSE, 
              clustering_method_columns = 'ward.D2',
              clustering_method_rows = 'ward.D2', 
              top_annotation = column_ha,
              right_annotation = row_ha,
              row_names_gp = gpar(fontsize = 8),
              row_names_side = "left", row_dend_side = 'right',
              
              heatmap_legend_param = list(title = expression(percent~of~live~cells),
                                          just=c('left', 'bottom')))
draw(hm, annotation_legend_side='top', padding = unit(c(2, 20, 2, 2), "mm"))
##dev.copy2pdf(file='figures/FAUST_heatmap_March2020.pdf', width=8, height=6)

