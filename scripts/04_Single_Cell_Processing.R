# Single cell analysis
# Processed with Cell Ranger 3.1

# Preprocessing in script - work from saved object

library(data.table)
library(dplyr)
library(Matrix)
library(DropletUtils)
library(scran)
library(biomaRt)
library(scater)
library(vegan)
library(PCAtools)
library(SNPRelate)
library(cowplot)
library(ggplot2)

# # set wd
# wd <- '~/Documents/R/onecarbon/'
# setwd(wd)
# 
# theme_set(theme_cowplot())
# 
# # read data - ascites and tumor samples
# asc <- read10xCounts('data/filtered_feature_bc_matrix_ascites/')
# tum <- read10xCounts('data/filtered_feature_bc_matrix_tumor/')
# 
# # merge
# sce_all <- cbind(asc, tum)
# sce_all <- calculateQCMetrics(sce_all)
# 
# # get HGNCs
# sce_all <- getBMFeatureAnnos(sce_all, 
#                              filters = "ensembl_gene_id",
#                              attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
#                              dataset = "hsapiens_gene_ensembl")
# 
# # rename with HGNC symbol
# colnames(sce_all) <- colData(sce_all)$Barcode
# rownames(sce_all) <- rowData(sce_all)$Symbol
# sce_all$Sample <- ifelse(grepl('tumor', sce_all$Sample), 'tumor', 'ascites')
# 
# # define MT genes
# mt_genes <- which(rowData(sce_all)$chromosome_name == "MT")
# 
# # general ribo genes: not all
# ribo_genes <- grepl("^RP[LS]", rowData(sce_all)$Symbol)
# ## OR: drop ribosomal genes by list: not much difference
# # ribos <- fread('data/ribosomal_genes_Dec3.txt')
# # sce_all <- sce_all[!(rownames(sce_all) %in% ribos$`Approved symbol`), ]
# 
# # set controls
# feature_ctrls <- list(mito = rownames(sce_all)[mt_genes],
#                       ribo = rownames(sce_all)[ribo_genes])
# 
# # QC metrics
# sce_all <- calculateQCMetrics(sce_all, feature_controls = feature_ctrls)
# 
# # Matching Donors:
# # use cellsnp -> vireo -> denconvolution
# adonor <- fread('data/vireo_out_ascites/donor_ids.tsv') %>% setkey(cell)
# cc <- colnames(sce_all[, sce_all$Sample=='ascites'])
# adonor <- adonor[cc, ]
# 
# tdonor <- fread('data/vireo_out_tumor/donor_ids.tsv') %>% setkey(cell)
# cc <- colnames(sce_all[, sce_all$Sample=='tumor'])
# tdonor <- tdonor[cc,]
# 
# ad <- paste0(adonor$donor_id, '_a')
# td <- paste0(tdonor$donor_id, '_t')
# 
# # set donor IDs
# colData(sce_all)$donor <- c(ad, td)
# 
# # set hard thresholds for inclusion
# # high MT threshold for tumor cells
# mito_thresh <- 15
# feat_thresh <- 500
# 
# ## Alternatively: drop cells based on MADs
# ## not done here
# # libsize.drop <- isOutlier(sce_all$total_counts, nmads=3, type="lower", log=TRUE)
# # feature.drop <- isOutlier(sce_all$total_features_by_counts, nmads=3, type="lower")
# # mito.drop <- isOutlier(sce_all$pct_counts_mito, nmads=3, type="higher")
# # drop
# # sce_all <- sce_all[, !(mito.drop)]
# 
# # now drop based on donor decon
# sce_all <- sce_all[, !(sce_all$donor=='doublet_a' | sce_all$donor=='doublet_t' |
#                        sce_all$donor=='unassigned_a' | sce_all$donor=='unassigned_t')]
# sce_all <- sce_all[, sce_all$total_features_by_counts > feat_thresh]
# sce_all <- sce_all[, sce_all$pct_counts_mito < mito_thresh]
# 
# # cluster for sizefactor
# clusters <- quickCluster(sce_all)
# sce_all <- scran::computeSumFactors(sce_all, clusters=clusters)
# sce_all <- logNormCounts(sce_all)
# 
# # model variances for HVGs
# # dec <- modelGeneVar(sce_all, block=sce_all$donor)
# # plot(dec$mean, dec$total)
# # top_hvgs <- getTopHVGs(dec, prop=0.1)
# 
# # run dim reduction
# # sce_all <- runPCA(sce_all, subset_row=top_hvgs)
# # sce_all <- runTSNE(sce_all, subset_row=top_hvgs)
# # sce_all <- runUMAP(sce_all, subset_row=top_hvgs)
# 
# # # ---------------- merge donors across samples -------------- #
# # # link donors via genotyping VCFs
# # # VCFs not included here
# # snpgdsVCF2GDS('data/vireo_out_tumor/GT_donors.vireo.vcf.gz', 
# #               'data/vireo_out_tumor/GT_donors.vireo.vcf.gds')
# # 
# # snpgdsVCF2GDS('data/vireo_out_ascites/GT_donors.vireo.vcf.gz', 
# #               'data/vireo_out_ascites/GT_donors.vireo.vcf.gds')
# # 
# # 
# # # rename samples to allow genotype intersections
# # fn <- 'data/vireo_out_ascites/GT_donors.vireo.vcf.gds'
# # 
# # f <- snpgdsOpen(fn, readonly=FALSE)
# # new_samp_id <- paste0("asc", 0:4 )
# # add.gdsn(f, "sample.id", new_samp_id, replace=TRUE)
# # closefn.gds(f)
# # 
# # # merge GDS
# # snpgdsCombineGeno(c('data/vireo_out_tumor/GT_donors.vireo.vcf.gds', 
# #                     'data/vireo_out_ascites/GT_donors.vireo.vcf.gds'), 
# #                   out.fn = 'data/GT_donors.vireo_merged.vcf.gds', 
# #                   same.strand = TRUE)
# # 
# # # open
# # gmerge <- snpgdsOpen('~/Documents/R/scrna/data/GT_donors.vireo_merged.vcf.gds')
# # 
# # # calc IBS
# # ibs <- snpgdsIBS(gmerge)
# # 
# # # get sample names and rename IBS matrix
# # gname <- read.gdsn(index.gdsn(gmerge, "sample.id"))
# # dimnames(ibs$ibs) <- list(gname, gname)
# # 
# # # look at relationships
# # # tumor 3 does not align well with with either ascites sample...
# # # Exclude unaligned for now
# # # pheatmap::pheatmap(ibs$ibs)
# # 
# # # manually annotate matched patients
# # pat <- sce_all$donor
# # samp <- dplyr::recode(pat, 
# #                       'donor3_a'='donorA_a', 
# #                       'donor0_t'='donorA_t', 
# #                       'donor0_a'='donorB_a', 
# #                       'donor1_t'='donorB_t',
# #                       'donor2_a'='donorC_a', 
# #                       'donor2_t'='donorC_t'
# #                       )
# # pat <- gsub('_.$', '', samp)
# 
# # update SCE
# sce_all$tsamp <- samp
# sce_all$pat <- pat
# 
# # subset
# sce_match <- sce_all[ , sce_all$pat=='donorA' | sce_all$pat=='donorB'| sce_all$pat=='donorC']
# 
# # renormalize
# clusters <- quickCluster(sce_match)
# sce_match <- scran::computeSumFactors(sce_match, clusters=clusters)
# sce_match <- logNormCounts(sce_match)
# 
# # run dim reds - all
# sce_match <- runPCA(sce_match)
# sce_match <- runTSNE(sce_match)
# sce_match <- runUMAP(sce_match)
# #saveRDS(sce_match, file='data/Ascites_Tumor_Cells_Matched_for_Metabo_Feb2020.rds')
# 
# # add clustering
# dec <- modelGeneVar(sce_match, block=sce_match$pat)
# # plot(dec$mean, dec$total)
# top_hvgs <- getTopHVGs(dec, prop=0.05)
# 
# sc_snn <- buildSNNGraph(sce_match[top_hvgs, ], k=40, type='jaccard', d=50)
# ccs <- igraph::cluster_louvain(sc_snn)
# sce_match$cl_cluster <- as.factor(igraph::membership(ccs))

#table(sce_match$cl_cluster)
#saveRDS(sce_match, file='data/Ascites_Tumor_Cells_Matched_for_Metabo_with_Cluster_Feb2020.rds')


