# Script to evaluate replicate measurements in QC runs

#----------------------------------------#
# ------------ configuration ------------#
#----------------------------------------#

library(MetaboAnalystR)
library(dplyr)
library(ggplot2)
library(cowplot)
library(limma)
library(ggrepel)
library(irr)
library(ggpubr)

wd <- '~/Documents/R/onecarbon/'
setwd(wd)

# initialize object and normalize
mSet <-InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, 
                      'data/20190710_Cell_Screening_1168_Final_QC_runs.csv', 
                      format='colu', 'disc')
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet) # use default half-min
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "SumNorm", "LogNorm", 
                      "AutoNorm", 
                      ratio=FALSE, 
                      ratioNum=20)

# get normalized data
ndat <- data.frame(t(mSet$dataSet$norm))

# SE function
std <- function(x) sd(x)/sqrt(length(x))

# calculate mean and SE of replicates
plot_df <- data.frame(ice_mean = apply(ndat[, grepl('Ice1e6', colnames(ndat))], 1, mean),
                      ice_se = apply(ndat[, grepl('Ice1e6', colnames(ndat))], 1, std),
                      norm_mean = apply(ndat[, grepl('Sep1e6', colnames(ndat))], 1, mean),
                      norm_se = apply(ndat[, grepl('Sep1e6', colnames(ndat))], 1, std)
                      )

# which names to show?
plot_df$name_keep <- ifelse(grepl('Adeno|NAD|Spermidine|S-Adeno|Inosi|Glutathione|Succin|Hypoxan|Citric|D-Glu|L-Lac|Folic|L-Glu|L-Meth', 
                                  rownames(plot_df)), 
                                  rownames(plot_df), 
                                  ''
                            )
# color vector
plot_df$col <- ifelse(plot_df$name_keep=='', rgb(0, 0, 0.8, 0.5), rgb(0.5, 0, 0, 0.5) )


#----------- stats for figure --------------#
# calculate ICCs

# for samples on ice
idat <- ndat[, grepl('Ice1e6', colnames(ndat))]
icc(idat, model='o')

# for samples not on ice -- column separated
mdat <- ndat[, grepl('Sep1e6', colnames(ndat))]
icc(mdat, model='o')

# Pearson correlation of means
cor.test(plot_df$ice_mean, plot_df$norm_mean)

#------------ plot -----------------#

ggplot(plot_df, aes(x = ice_mean, y = norm_mean)) + 
  theme_cowplot() + geom_abline(slope=1, intercept=0, color='light gray') +
  geom_errorbarh(aes(xmin = ice_mean - ice_se, xmax = ice_mean + ice_se), color = rgb(0, 0, 0.8, 0.2)) + 
  geom_errorbar(aes(ymin = norm_mean - norm_se, ymax = norm_mean + norm_se), color = rgb(0, 0, 0.8, 0.2) ) +
  geom_point(color=plot_df$col) +
  scale_color_manual('', values = point_cols) +
  geom_smooth(method='lm', fill = NA, color = rgb(0,0,1,0.7)) +
  ylab(expression(Column~separation~log[10]~IC~(ICC==0.76))) +
  xlab(expression(Iced~cells~log[10]~IC~(ICC==0.87))) +
  annotate('text', x = -1.5, y = 1.3, label = expression(r==0.77)) +
  geom_text_repel(aes(label = name_keep), size = 3, force = 20) 
  
dev.copy2pdf(file='figures/prerun_QC_scatter_plot_revised.pdf', height=5, width=5)





