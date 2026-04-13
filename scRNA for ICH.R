# setwd(choose.dir())
setwd("C://Users//qingm//Desktop//单细胞")
rm(list = ls())

#### library包 ####
library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CCA)
library(cluster)
library(cowplot)
library(monocle3)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
library(decontX)
library(scDblFinder)
library(DoubletFinder)
library(clustree)
library(SingleR)
library(clusterProfiler)
library(org.Rn.eg.db)
library(ggsci)
library(slingshot)
library(RColorBrewer)
library(tradeSeq)
library(CellChat)
library(NMF)
library(future)
library(future.apply)
library(ggalluvial)

#
#### 安装 ####
# chooseBioCmirror()
# chooseCRANmirror()

## install.packages(c("remotes", "RANN", "RColorBrewer", "cowplot", "Seurat", "SeuratObject", "fitdistrplus", "future", "future.apply", "ggplot2", "ggrepel", "ggridges", "httr", "ica", "irlba", "jsonlite", "leiden", "lmtest", "matrixStats", "miniUI", "patchwork", "rlang", "RMTstat", "R.utils", "RobustRankAggreg", "roxygen2", "stringr", "tibble", "tidyr", "tidytree", "VAM", "harmony", "hdf5r", "assertthat", "randomcoloR", "NMF", "AnnoProbe", "clustermole", "fastDummies", "qs", "igraph", "scCustomize", "paletteer", "miscTools", "SCINA", "phylogram"))

## BiocManager::install(c("DOSE", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "limma", "edgeR", "maftools", "ComplexHeatmap", "preprocessCore", "DESeq2", "TCGAbiolinks", "sva", "GenomicFeatures", "AUCell", "BiocParallel", "ComplexHeatmap", "decoupleR", "fgsea", "ggtree", "GSEABase", "GSVA", "Nebulosa", "scde", "singscore", "SummarizedExperiment", "viper", "SparseMatrixStats", "celldex", "tricycle", "GeneOverlap", "TOAST", "infercnv", "metapod", "bluster", "scran", "celda", "glmGamPoi", "SingleR", "wk", "units", "s2", "classInt", "warp", "biglm", "sf", "spData", "slider", "ResidualMatrix", "speedglm", "spdep", "rsample", "pscl", "grr", "batchelor", "scDblFinder"))

# usethis::create_github_token()
# usethis::edit_r_environ()
# devtools::install_github("GfellerLab/EPIC", build_vignettes = TRUE, force = TRUE)
# devtools::install_github("Moonerss/CIBERSORT", force = TRUE)
# devtools::install_github('dviraran/xCell', force = TRUE)
# devtools::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source", force = TRUE)
# devtools::install_github("cansysbio/ConsensusTME", force = TRUE)
# devtools::install_github("cit-bioinfo/mMCP-counter", force = TRUE)
# devtools::install_github("omnideconv/immunedeconv", force = TRUE)
# devtools::install_github("Hy4m/linkET", force = TRUE)
# devtools::install_github("chuiqin/irGSEA", force = T)
# devtools::install_github("jinworks/CellChat", force = T)
# devtools::install_github("jokergoo/circlize", force = T)
# devtools::install_github("navinlabcode/copykat", force = T)
# devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', force = T)
# devtools::install_github('immunogenomics/presto', force = T)
# devtools::install_github("campbio/decontX", force = T)
# devtools::install_github(repo = "powellgenomicslab/scPred", ref = "9f407b7436f40d44224a5976a94cc6815c6e837f", force = T)
# devtools::install_github(repo = 'genecell/COSGR', force = T)
# devtools::install_github("JerryZhang-1222/starTracer", force = T)
# devtools::install_github('cole-trapnell-lab/monocle3', force = T)
# devtools::install_github('smorabit/hdWGCNA', ref = 'dev', force = T)
# devtools::install_github('xuranw/Music', force = T)
# devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org", force = T)
# devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org", force = T)
# devtools::install_github("jackbibby1/SCPA", force = T)


#### 1. 10X数据读取 ####
data_dir <- paste0(getwd(), "/GSE282740")
samples = list.files(data_dir)
dataGSE <- gsub("(.*?)\\_.*", "\\1", samples)
lapply(levels(factor(dataGSE)), function(id){
  dir.create(paste0(data_dir, "/", id))  
})
lapply(levels(factor(dataGSE)), function(id){
  files0 <- list.files(path = paste0(data_dir, "/"), pattern = paste0(id, ".*barcodes.*"), full.names = TRUE)
  file.rename(from = files0, to = paste0(data_dir, "/", id, "/barcodes.tsv.gz"))
  files1 <- list.files(path = paste0(data_dir, "/"), pattern = paste0(id, ".*features.*"), full.names = TRUE)
  file.rename(from = files1, to = paste0(data_dir, "/", id, "/features.tsv.gz"))
  files2 <- list.files(path = paste0(data_dir, "/"), pattern = paste0(id, ".*matrix.*"), full.names = TRUE)
  file.rename(from = files2, to = paste0(data_dir, "/", id, "/matrix.mtx.gz"))
})
samples = list.files(data_dir)
dir = file.path(data_dir, samples)
scRNAdata <- Read10X(data.dir = dir)
scRNA <- CreateSeuratObject(counts = scRNAdata, 
                            min.cells = 3, 
                            min.features = 200)

scRNAidens = mapvalues(Idents(scRNA), from = levels(Idents(scRNA)), to = samples)
Idents(scRNA) = scRNAidens
scRNA$Type = Idents(scRNA)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^Mt-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP")
scRNA[["percent.hb"]] <- PercentageFeatureSet(scRNA, pattern = "^Hba|^Hbb")
qsave(scRNA, "scRNA.qs", nthreads = detectCores())


#### 2. 单细胞QC ####
scRNA = qread("scRNA.qs", nthreads = detectCores())
pdf(file = "1.1 QC前小提琴图.pdf", width = 13, height = 8)
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
dev.off()
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt") + RotatedAxis()
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + RotatedAxis()
pdf(file = "1.2 QC前指标相关性.pdf", width = 20, height = 8)
plot2 + plot3 + plot_layout(ncol = 2)
dev.off()

mask1 <- scRNA$nFeature_RNA > 200 & scRNA$nFeature_RNA < 8000 
mask2 <- scRNA$percent.mt < 15
mask3 <- scRNA$percent.rb < 40
mask4 <- scRNA$percent.hb < 5
mask <- mask1 & mask2 & mask3 & mask4
scRNA_out <- scRNA[,mask]
pdf(file = "1.3 QC后小提琴图.pdf", width = 13, height = 8)
VlnPlot(scRNA_out, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
dev.off()
plot2 <- FeatureScatter(scRNA_out, feature1 = "nCount_RNA", feature2 = "percent.mt") + RotatedAxis()
plot3 <- FeatureScatter(scRNA_out, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + RotatedAxis()
pdf(file = "1.4 QC后指标相关性.pdf", width = 20, height = 8)
plot2 + plot3 + plot_layout(ncol = 2)
dev.off()

scRNA <- scRNA_out
qsave(scRNA, "scRNA数据质控后.qs", nthreads = detectCores())


#### 3. 单细胞标准化 ####
scRNA = qread("scRNA数据质控后.qs", nthreads = detectCores())
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
pdf(file = "2.1.1 高变基因.pdf", width = 7, height = 6)
plot2
dev.off()
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- Seurat::RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- Seurat::RunTSNE(scRNA, dims = 1:20)
scRNA <- Seurat::RunUMAP(scRNA, dims = 1:20)
pdf(file = "2.1.2 tSNE图.pdf", width = 7.5, height = 5.5)
DimPlot(scRNA, reduction = "tsne", pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "2.1.3 PCA图.pdf", width = 7.5, height = 5.5)
DimPlot(scRNA, reduction = "pca", pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "2.1.4 UMAP图.pdf", width = 7.5, height = 5.5)
DimPlot(scRNA, reduction = "umap", pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "2.1.5 每个样本tSNE图.pdf", width = 12, height = 4)
do_DimPlot(sample = scRNA,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2), split.by = "Type", pt.size = 0.01
) + 
  labs(color = "combine")
dev.off()

counts <- scRNA@assays$RNA@layers$counts
decontxcounts <- decontX(counts)
scRNA$Contamination = decontxcounts$contamination
scRNAde = scRNA[,scRNA$Contamination < 0.3]
scRNAde <- RunHarmony(scRNAde, group.by.vars = "Type")
scRNAde <- FindVariableFeatures(scRNAde, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNAde), 10)
plot1 <- VariableFeaturePlot(scRNAde)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
pdf(file = "2.2.1 去批次后高变基因.pdf", width = 7, height = 6)
plot2
dev.off()
pdf(file = "2.2.2 Harmony去批次图.pdf", width = 7.5, height = 5.5)
DimPlot(scRNAde, reduction = "harmony", pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
scRNAde <- Seurat::RunTSNE(scRNAde, dims = 1:20, reduction = 'harmony')
scRNAde <- Seurat::RunUMAP(scRNAde, dims = 1:20, reduction = 'harmony')
pdf(file = "2.2.3 去批次后tSNE图.pdf", width = 7.5, height = 5.5)
DimPlot(scRNAde, reduction = "tsne", pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "2.2.4 去批次后UMAP图.pdf", width = 7.5, height = 5.5)
DimPlot(scRNAde, reduction = "umap", pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "2.2.5 去批次后每个样本tSNE图.pdf", width = 12, height = 7.5)
do_DimPlot(sample = scRNAde,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2), split.by = "Type", pt.size = 0.01
) + 
  labs(color = "combine")
dev.off()

pdf(file = "2.3 前15个PC可视化.pdf", width = 7.5, height = 9)
DimHeatmap(scRNAde, dims = 1:15, cells = 1000, balanced = T)
dev.off()
scRNAde <- JackStraw(scRNAde, num.replicate = 100)
scRNAde <- ScoreJackStraw(scRNAde, dims = 1:20)
pdf(file = "2.4 jackstrawplot.pdf", width = 7.5, height = 5.5)
JackStrawPlot(scRNAde, dims = 1:20)
dev.off()
pdf(file = "2.5 ElbowPlot.pdf", width = 5, height = 4)
ElbowPlot(scRNAde, ndims = 30)
dev.off()
qsave(scRNAde, "scRNA标准化后.qs", nthreads = detectCores())


#### 4. 双细胞去除 ####
scRNA = qread("scRNA标准化后.qs", nthreads = detectCores())
## DoubletFinder ##
# scRNAsp <- SplitObject(scRNA, split.by = "Type")
# scRNAdf <- lapply(scRNAsp, function(scRNAdfs){
#   sweeplist <- paramSweep(scRNAdfs, PCs = 1:20, sct = T)
#   sweepstats <- summarizeSweep(sweeplist, GT = F)
#   sweepstats[order(sweepstats$BCreal),]
#   bcmvn <- find.pK(sweepstats)
#   pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
#   homotypic.prop <- modelHomotypic(scRNAdfs$Type)
#   nExp_poi <- round(0.05 * nrow(scRNAdfs@meta.data))
#   nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))
#   scRNAdfs <- doubletFinder(scRNAdfs, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
#                             nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
#   pss <- "DF\\.classifications"
#   colnames(scRNAdfs@meta.data) <- ifelse(grepl(pss, colnames(scRNAdfs@meta.data)), "DF", colnames(scRNAdfs@meta.data))
#   return(scRNAdfs)
# })
# scRNA <- Reduce(function(x, y) merge(x, y, all=T), scRNAdf)
# scRNAdfs <- subset(scRNAdfs, subset = (DF == "Singlet"))
# qsave(scRNA, "scRNA剔除双细胞后2.qs", nthreads = detectCores())

## scDblFinder ##
scRNAdbin <- as.SingleCellExperiment(scRNA)
scRNAdb <- scDblFinder(scRNAdbin, sample = "Type")
scRNA <- as.Seurat(scRNAdb)
pdf("3.1 双细胞分布图.pdf")
DimPlot(scRNA, group.by = "scDblFinder.class", raster = F)
dev.off()
scRNA <- subset(scRNA, subset = (scDblFinder.class == "singlet"))
qsave(scRNA, "scRNA剔除双细胞后.qs", nthreads = detectCores())


#### 5. 双细胞后标准化 ####
scRNA = qread("scRNA剔除双细胞后.qs", nthreads = detectCores())
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- Seurat::RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- Seurat::RunTSNE(scRNA, dims = 1:20)
scRNA <- Seurat::RunUMAP(scRNA, dims = 1:20)
scRNA <- RunHarmony(scRNA, group.by.vars = "Type")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- Seurat::RunTSNE(scRNA, dims = 1:20, reduction = 'harmony')
scRNA <- Seurat::RunUMAP(scRNA, dims = 1:20, reduction = 'harmony')
scRNA <- Seurat::RunPCA(scRNA, features = VariableFeatures(object = scRNA))
qsave(scRNA, "scRNA常规标准化后.qs", nthreads = detectCores())


#### 6. 单细胞分群聚类  ####
scRNA = qread("scRNA常规标准化后.qs", nthreads = detectCores())
scRNAPC = 25
scRNA = FindNeighbors(scRNA, dims = 1:scRNAPC, reduction = "harmony")
for (res in seq(0, 0.1, by = 0.01)){
  scRNA = FindClusters(scRNA, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
apply(scRNA@meta.data[,grep("RNA_snn_res", colnames(scRNA@meta.data))], 2, table)
p2_tree = clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
pdf(file = "4.1 0至0.1节点分辨率.pdf", width = 12, height = 10)
p2_tree
dev.off()
scRNA = FindNeighbors(scRNA, dims = 1:scRNAPC, reduction = "harmony")
scRNA <- FindClusters(scRNA, resolution = 0.06)

scRNA <- RunUMAP(scRNA, dims = 1:scRNAPC, reduction = 'harmony')
scRNA <- RunTSNE(scRNA, dims = 1:scRNAPC, reduction = 'harmony')
pdf(file = "4.2 选定PC降维UMAP图.pdf", width = 6.5, height = 5.5)
DimPlot(scRNA, group.by = "RNA_snn_res.0.06", reduction = "umap", label = T, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "4.3 选定PC降维tSNE图.pdf", width = 6.5, height = 5.5)
DimPlot(scRNA, group.by = "RNA_snn_res.0.06", reduction = "tsne", label = T, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
qsave(scRNA, "scRNA分群聚类后.qs", nthreads = detectCores())


#### 7. 单细胞注释  ####
seurat = qread("scRNA分群聚类后.qs", nthreads = detectCores())
## SingleR
refdata <- celldex::MouseRNAseqData()
testdata <- GetAssayData(seurat)
cellpred <- SingleR(test = testdata,
                    ref = refdata,
                    labels = refdata$label.main)
celltype = data.frame(ClusterID = rownames(cellpred),
                      celltype = cellpred$labels,
                      stringsAsFactors = F)
write.csv(celltype, "5.1 SingleR细胞注释.csv")
pdf(file = "5.2 SingleR热图.pdf", width = 7.5, height = 5.5)
plotScoreHeatmap(cellpred, clusters = rownames(celltype), order.by = "cluster")
dev.off()

newLabels = cellpred$labels
seuratidens = mapvalues(Idents(seurat), from = Idents(seurat), to = newLabels)
Idents(seurat) = seuratidens
seurat$cellType = Idents(seurat)
pdf(file = "5.3 注释UMAP.pdf", width = 7.5, height = 5.5)
DimPlot(seurat, group.by = "cellType", reduction = "umap", label = F, label.box = F, repel = T, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "5.4 注释tSNE.pdf", width = 7.5, height = 5.5)
DimPlot(seurat, group.by = "cellType", reduction = "tsne", label = F, label.box = F, repel = T, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
qsave(seurat, "SingleR注释.qs", nthreads = detectCores())

seurat = qread("SingleR注释.qs", nthreads = detectCores())
# 谷氨酸能（例如，Slc17a7和Slc17a6）和泛GABA能（例如，Gadl、Gad2和Slc32a1）标志
seuratgenes = c("Nrg1", "Nrg3", "Psap", "Grn", "Spp1", "Igf1")
pdf(file = "5.5 基因分布UMAP.pdf", width = 12, height = 12)
FeaturePlot(seurat, features = seuratgenes, reduction = "umap", cols = c("grey", "red"), min.cutoff = 0, max.cutoff = 1, ncol = 2, pt.size = 0.01)
dev.off()
pdf(file = "5.6 基因分布tSNE.pdf", width = 12, height = 12)
FeaturePlot(seurat, features = seuratgenes, reduction = "tsne", cols = c("grey", "red"), min.cutoff = 0, max.cutoff = 1, ncol = 2, pt.size = 0.01)
dev.off()


#### 8. 计算细胞比例 ####
scRNA = qread("SingleR注释.qs", nthreads = detectCores())
cellnum <- table(scRNA$cellType, scRNA$Type)
cell.prop <- as.data.frame(prop.table(cellnum))
colnames(cell.prop) <- c("Celltype", "Group", "Proportion")
group.map <- c(GSM8649646 = "Sham", GSM8649647 = "ICH")
cell.prop$Group <- group.map[as.character(cell.prop$Group)]
cell.prop$Group <- factor(cell.prop$Group, levels = c("Sham", "ICH"))
cols <- pal_npg("nrc")(10)
p.bar <- ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = cols) + 
  ggtitle("Cell proportion") + 
  theme_bw() + 
  theme(axis.ticks.length = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(title = NULL))
pdf("6.1 细胞比例.pdf", width = 4, height = 4)
p.bar
dev.off()

#### 9. 计算marker基因 ####
scRNA = qread("SingleR注释.qs", nthreads = detectCores())
# FindAllMarkers算法
scRNA.markers <- FindAllMarkers(scRNA, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(scRNA.markers, file = "7.1 每个细胞Marker基因.csv")

# top50 <- scRNA.markers %>%
#   group_by(cluster) %>%
#   slice_max(avg_log2FC, n = 50)
# write.csv(top50, "top50.csv", row.names = FALSE)

top5scRNA.markers <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
col <- c(ggsci::pal_npg()(9), ggsci::pal_jco()(9), ggsci::pal_jama()(7), ggsci::pal_nejm()(8))
pdf(file = "7.2 每个细胞前5基因热图.pdf", width = 22, height = 16)
DoHeatmap(scRNA, features = top5scRNA.markers$gene, 
          group.colors = col) + 
  ggsci::scale_color_npg() + 
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033',
                       name = 'Z-score')
dev.off()


#### 10. KEGG/GO分析 ####
scRNA = qread("SingleR注释.qs", nthreads = detectCores())
sample.markers <- FindAllMarkers(scRNA, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top200 <- sample.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
top200_gene <- unstack(top200, gene ~ cluster)
top200_entrez <- lapply(X = top200_gene, FUN = function(x){
  x <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Rn.eg.db")
  x <- x[,-1]
})

# GO
all_cc <- compareCluster(top200_entrez, fun = 'enrichGO', ont = 'CC', OrgDb = 'org.Rn.eg.db')
write.csv(as.data.frame(all_cc), "8.1 GO分析_cc.csv")
pdf('8.2 GO分析_cc.pdf', width = 15, height = 15)
dotplot(all_cc, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

all_mf <- compareCluster(top200_entrez, fun = 'enrichGO', ont = 'MF', OrgDb = 'org.Rn.eg.db')
write.csv(as.data.frame(all_mf), "8.3 GO分析_mf.csv")
pdf('8.4 GO分析_mf.pdf', width = 15, height = 15)
dotplot(all_mf, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

all_bp <- compareCluster(top200_entrez, fun = 'enrichGO', ont = 'BP', OrgDb = 'org.Rn.eg.db')
write.csv(as.data.frame(all_bp), "8.5 GO分析_bp.csv")
pdf('8.6 GO分析_bp.pdf', width = 15, height = 15)
dotplot(all_bp, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

# 批量对每个细胞GO分析
lapply(1:length(top200_entrez), function(i){
  all_go_cell <- enrichGO(gene = top200_entrez[[i]],
                          OrgDb = 'org.Rn.eg.db',
                          ont = 'ALL',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  all_go_cell <- setReadable(all_go_cell, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(all_go_cell), paste0("8. ", names(top200_entrez[i]), " GO分析.csv"))
  plotGO <- dotplot(all_go_cell, split = "ONTOLOGY", label_format = 100) + facet_grid(ONTOLOGY~., scale = "free")
  pdf(paste0('8. ', names(top200_entrez[i]), ' GO分析.pdf'), width = 16, height = 8)
  print(plotGO)
  dev.off()
})

# KEGG
all_kegg <- compareCluster(top200_entrez, fun = 'enrichKEGG', pvalueCutoff = 0.05, organism = "rno")
write.csv(as.data.frame(all_kegg), "8.7 KEGG分析.csv")
pdf('8.8 KEGG分析.pdf', width = 15, height = 15)
dotplot(all_kegg, showCategory = 5, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

# 批量对每个细胞KEGG分析
lapply(1:length(top200_entrez), function(i){
  kegg_cell <- enrichKEGG(gene = top200_entrez[[i]],
                          organism = "rno",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1)
  kegg_cell <- setReadable(kegg_cell, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(kegg_cell), paste0("8. ", names(top200_entrez[i]), " KEGG分析.csv"))
  plotKEGG <- dotplot(kegg_cell, label_format = 100)
  pdf(paste0('8. ', names(top200_entrez[i]), ' KEGG分析.pdf'), width = 12, height = 6)
  print(plotKEGG)
  dev.off()
})


#### 11. 亚群分群聚类 ####
scRNA = qread("SingleR注释.qs", nthreads = detectCores())
scRNAf <- SplitObject(scRNA, split.by = "cellType")

# Oligo
subset <- scRNAf$Oligodendrocytes
subset <- NormalizeData(subset, normalization.method = "LogNormalize", scale.factor = 10000)
subset <- FindVariableFeatures(subset, selection.method = "vst", nfeatures = 2000)
subsetall.genes <- rownames(subset)
subset <- ScaleData(subset, features = subsetall.genes)
subset <- Seurat::RunPCA(subset, features = VariableFeatures(object = subset))
subset <- Seurat::RunTSNE(subset, dims = 1:20)
subset <- Seurat::RunUMAP(subset, dims = 1:20)
subset = FindNeighbors(subset, dims = 1:20, reduction = "harmony")
for (res in seq(0.1, 1.0, by = 0.1)){
  subset = FindClusters(subset, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
apply(subset@meta.data[,grep("RNA_snn_res", colnames(subset@meta.data))], 2, table)
p2_tree = clustree(subset@meta.data, prefix = "RNA_snn_res.")
pdf(file = "9.1.1 亚群分辨率_Oligo.pdf", width = 12, height = 12)
p2_tree
dev.off()
subset <- FindClusters(subset, resolution = 0.1)

# newLabels <- c("Oligo_0", "Oligo_1", "Oligo_2", "Oligo_3")
pdf(file = "9.1.2 亚群聚类UMAP_Oligo.pdf", width = 6.5, height = 5.5)
DimPlot(subset, group.by = "RNA_snn_res.0.1", reduction = "umap", label = F, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "9.1.3 亚群聚类tSNE_Oligo.pdf", width = 6.5, height = 5.5)
DimPlot(subset, group.by = "RNA_snn_res.0.1", reduction = "tsne", label = F, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
qsave(subset, "亚群聚类_Oligo.qs", nthreads = detectCores())

#注释
genes <- list("POPC" = c("Mki67", "Top2a"),
              "OPC" = c("Pdgfra", "Neu4"),
              "NFOL" = c("Fyn", "Itpr2"),
              "MFOL" = c("Ccp110", "Opalin"),
              "MOL" = c("Mog", "Mal"),
              "DOL" = c("Serpina3n", "C4b"))
pdf(file = "9.1.4 亚群手动注释_Oligo.pdf", width = 10, height = 8)
do_DotPlot(sample = subset, features = genes, dot.scale = 12, legend.framewidth = 2, font.size = 10)
dev.off()
ann.ids <- c("MOL", #0
             "MFOL", #1
             "NFOL", #2
             "OPC" #3
)
subsetidens = mapvalues(Idents(subset), from = levels(Idents(subset)), to = ann.ids)
Idents(subset) = subsetidens
subset$cellType = Idents(subset)
pdf(file = "9.1.5 亚群手动注释后_Oligo.pdf", width = 10, height = 8)
do_DotPlot(sample = subset, features = genes, dot.scale = 12, legend.framewidth = 2, font.size = 10)
dev.off()

#注释可视化
pdf(file = "9.1.6 亚群注释UMAP_Oligo.pdf", width = 7.5, height = 5.5)
DimPlot(subset, group.by = "cellType", reduction = "umap", label = F, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
pdf(file = "9.1.7 亚群注释tSNE_Oligo.pdf", width = 7.5, height = 5.5)
DimPlot(subset, group.by = "cellType", reduction = "tsne", label = F, label.size = 3.5, pt.size = 0.01) + theme_classic() + theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid", legend.position = "right"))
dev.off()
qsave(subset, "亚群注释_Oligo.qs", nthreads = detectCores())

#细胞比例
subset = qread("亚群注释_Oligo.qs", nthreads = detectCores())
cellnum <- table(subset$cellType, subset$Type)
cell.prop <- as.data.frame(prop.table(cellnum))
colnames(cell.prop) <- c("Celltype", "Group", "Proportion")
group.map <- c(GSM8649646 = "Sham", GSM8649647 = "ICH")
cell.prop$Group <- group.map[as.character(cell.prop$Group)]
cell.prop$Group <- factor(cell.prop$Group, levels = c("Sham", "ICH"))
cols <- pal_npg("nrc")(10)
p.bar <- ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = cols) + 
  ggtitle("Cell proportion") + 
  theme_bw() + 
  theme(axis.ticks.length = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(title = NULL))
pdf("9.1.8 亚群细胞比例_Oligo.pdf", width = 4, height = 4)
p.bar
dev.off()

#计算marker
subset.markers <- FindAllMarkers(subset, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top5subset.markers <- subset.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
col <- c(ggsci::pal_npg()(9), ggsci::pal_jco()(9), ggsci::pal_jama()(7), ggsci::pal_nejm()(8))
pdf(file = "9.1.9 亚群前5基因热图_Oligo.pdf", width = 22, height = 16)
DoHeatmap(subset, features = top5subset.markers$gene, 
          group.colors = col) + 
  ggsci::scale_color_npg() + 
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033',
                       name = 'Z-score')
dev.off()

# 高变基因分布
subsetgenes = c("Erbb4", "Gpr37", "Sort1", "Itgav", "Itgb1", "Igf1r")
pdf(file = "9.1.10 基因分布UMAP.pdf", width = 12, height = 12)
FeaturePlot(subset, features = subsetgenes, reduction = "umap", cols = c("grey", "red"), min.cutoff = 0, max.cutoff = 1, ncol = 2, pt.size = 0.01)
dev.off()
pdf(file = "9.1.11 基因分布tSNE.pdf", width = 12, height = 12)
FeaturePlot(subset, features = subsetgenes, reduction = "tsne", cols = c("grey", "red"), min.cutoff = 0, max.cutoff = 1, ncol = 2, pt.size = 0.01)
dev.off()


#### 12. 亚群KEGG/GO ####
scRNA = qread("亚群注释_Oligo.qs", nthreads = detectCores())
sample.markers <- FindAllMarkers(scRNA, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top200 <- sample.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
top200_gene <- unstack(top200, gene ~ cluster)
top200_entrez <- lapply(X = top200_gene, FUN = function(x){
  x <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Rn.eg.db")
  x <- x[,-1]
})

# GO
all_cc <- compareCluster(top200_entrez, fun = 'enrichGO', ont = 'CC', OrgDb = 'org.Rn.eg.db')
write.csv(as.data.frame(all_cc), "10.1 GO分析_cc_Oligo.csv")
pdf('10.2 GO分析_cc_Oligo.pdf', width = 15, height = 15)
dotplot(all_cc, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

all_mf <- compareCluster(top200_entrez, fun = 'enrichGO', ont = 'MF', OrgDb = 'org.Rn.eg.db')
write.csv(as.data.frame(all_mf), "10.3 GO分析_mf_Oligo.csv")
pdf('10.4 GO分析_mf_Oligo.pdf', width = 15, height = 15)
dotplot(all_mf, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

all_bp <- compareCluster(top200_entrez, fun = 'enrichGO', ont = 'BP', OrgDb = 'org.Rn.eg.db')
write.csv(as.data.frame(all_bp), "10.5 GO分析_bp_Oligo.csv")
pdf('10.6 GO分析_bp_Oligo.pdf', width = 15, height = 15)
dotplot(all_bp, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

# 批量对每个细胞GO分析
lapply(1:length(top200_entrez), function(i){
  all_go_cell <- enrichGO(gene = top200_entrez[[i]],
                          OrgDb = 'org.Rn.eg.db',
                          ont = 'ALL',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  all_go_cell <- setReadable(all_go_cell, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(all_go_cell), paste0("10. ", names(top200_entrez[i]), " GO分析.csv"))
  plotGO <- dotplot(all_go_cell, split = "ONTOLOGY", label_format = 100) + facet_grid(ONTOLOGY~., scale = "free")
  pdf(paste0('10. ', names(top200_entrez[i]), ' GO分析.pdf'), width = 16, height = 8)
  print(plotGO)
  dev.off()
})

# KEGG
all_kegg <- compareCluster(top200_entrez, fun = 'enrichKEGG', pvalueCutoff = 0.05, organism = "rno")
write.csv(as.data.frame(all_kegg), "10.7 KEGG分析_Oligo.csv")
pdf('10.8 KEGG分析_Oligo.pdf', width = 15, height = 15)
dotplot(all_kegg, showCategory = 5, includeAll = F, label_format = 100) + theme(axis.text = element_text(angle = 30, hjust = 1))
dev.off()

# 批量对每个细胞KEGG分析
lapply(1:length(top200_entrez), function(i){
  kegg_cell <- enrichKEGG(gene = top200_entrez[[i]],
                          organism = "rno",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1)
  kegg_cell <- setReadable(kegg_cell, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(kegg_cell), paste0("10. ", names(top200_entrez[i]), " KEGG分析.csv"))
  plotKEGG <- dotplot(kegg_cell, label_format = 100)
  pdf(paste0('10. ', names(top200_entrez[i]), ' KEGG分析.pdf'), width = 12, height = 6)
  print(plotKEGG)
  dev.off()
})


#### 13. 拟时序分析 ####
# Slingshot
scRNA = qread("亚群注释_Oligo.qs", nthreads = detectCores())
sce <- as.SingleCellExperiment(scRNA)
sce <- slingshot(data = sce, clusterLabels = 'cellType', reducedDim = 'UMAP', start.clus = NULL, end.clus = NULL)
celltypes <- unique(sce$cellType)
nature_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#FFB000", "#A20056")
pal <- setNames(nature_colors[1:length(celltypes)], celltypes)

pdf("11.1 Slingshot细胞轨迹_Oligo.pdf", width = 7 ,height = 6)
plot(reducedDims(sce)$UMAP, col = pal[sce$cellType], cex = 0.5, pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
legend("topright", legend = names(pal), col = pal, pch = 16, cex = 0.6, box.lwd = 0.5)
title("Slingshot Trajectory Inference", cex.main = 1.2, font.main = 2)
dev.off()

pdf("11.2 Slingshot拟时序_1_Oligo.pdf", width = 7, height = 6)
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks = 100)]
plotcol[is.na(plotcol)] <- "#CCCCCC"
plot(reducedDims(sce)$UMAP, col = plotcol, pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 2, col = brewer.pal(6, "Set2"))
legend("topright", legend = paste0("lineage", 1:2), col = unique(brewer.pal(4, "Set2")),
       inset = 0.8, pch = 16, cex = 0.6, box.lwd = 0.5)
title("Slingshot Pseudotime_1", cex.main = 1.2, font.main = 2)
dev.off()

pdf("11.3 Slingshot拟时序_2_Oligo.pdf", width = 7, height = 6)
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_2, breaks = 100)]
plotcol[is.na(plotcol)] <- "#CCCCCC"
plot(reducedDims(sce)$UMAP, col = plotcol, pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 2, col = brewer.pal(6, "Set2"))
legend("topright", legend = paste0("lineage", 1:2), col = unique(brewer.pal(4, "Set2")),
       inset = 0.8, pch = 16, cex = 0.6, box.lwd = 0.5)
title("Slingshot Pseudotime_2", cex.main = 1.2, font.main = 2)
dev.off()

diff.wilcox <- FindAllMarkers(scRNA, only.pos = T, min.pct = 0.1, thresh.use = 0.1)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val < 0.05)
diff.genes <- unique(subset(all.markers, p_val_adj < 0.001)$gene)
pseudotimeED <-slingPseudotime(sce, na = F)
cellWeightsED <- slingCurveWeights(sce)
counts <- scRNA@assays$RNA$counts
counts <- counts[diff.genes, ]
sce_slingshot <- fitGAM(counts = as.matrix(counts), pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T)
qsave(sce_slingshot, "Slingshot拟时序_Oligo.qs", nthreads = detectCores())

sce_slingshot = qread("Slingshot拟时序_Oligo.qs", nthreads = detectCores())
gene_dynamic <- list()
genes_plot <- c("Erbb4", "Gpr37", "Sort1", "Itgb1", "Igf1r")
pdf("11.4 拟时序基因表达.pdf", width = 10, height = 8)
for(i in 1:length(genes_plot)){
  p = plotSmoothers(sce_slingshot, assays(sce_slingshot)$counts, gene = genes_plot[i], 
                    alpha = 0.6, border = T, lwd = 2) + ggtitle(genes_plot[i])
  gene_dynamic[[i]] <- p
}
Seurat::CombinePlots(gene_dynamic, ncol = 2)
dev.off()

assoRes <- associationTest(sce_slingshot)
startRes <- startVsEndTest(sce_slingshot)

gene_dynamic <- list()
genes_plot <- c("Erbb4", "Gpr37", "Sort1", "Itgb1", "Igf1r")
pdf("11.5 每条lineage拟时序基因表达.pdf", width = 10, height = 10)
for(i in 1:length(genes_plot)){
  lineage_plots <- list()
  for(lineage_num in 1:2){
    p <- plotSmoothers(sce_slingshot, assays(sce_slingshot)$counts, gene = genes_plot[i], 
                       lineage = lineage_num, alpha = 0.6, border = T, lwd = 2) + 
      ggtitle(paste0(genes_plot[i], " - Lineage ", lineage_num)) + theme(legend.position = "none")
    lineage_plots[[lineage_num]] <- p
  }
  gene_dynamic[[i]] <- cowplot::plot_grid(plotlist = lineage_plots, ncol = 2)
}
cowplot::plot_grid(plotlist = gene_dynamic, ncol = 1)
dev.off()


#### 14. CellChat细胞通讯 ####
scRNA = qread("SingleR注释.qs", nthreads = detectCores())
cellchat <- createCellChat(scRNA@assays$RNA@data)
meta <- data.frame(cellType = scRNA$cellType, row.names = Cells(scRNA))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
cellchat <- setIdent(cellchat, ident.use = "cellType")
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

pdf("12.1.1 细胞通讯网络_细胞交互次数.pdf")
par(mfrow = c(1,1), xpd = T)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
dev.off()
pdf("12.1.2 细胞通讯网络_细胞交互权重.pdf")
par(mfrow = c(1, 1), xpd = T)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
dev.off()
mat <- cellchat@net$weight
pdf("12.1.3 不同细胞通讯网络权重.pdf")
par(mfrow = c(3, 3), xpd = T, mar = c(1, 1, 4, 1))
for (i in 1: nrow(mat)){
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, vertex.label.cex = 0.8, weight.scale = T, label.edge = F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

## NRG
cellchat@netP$pathways
  pathways.show <- "NRG"
  pdf("12.2.1 重要通路通讯_贝壳图_NRG.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  pdf("12.2.2 重要通路通讯_弦图_NRG.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  pdf("12.2.3 通讯热图_NRG.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
  
  pdf("12.3.1 配体受体贡献度_NRG.pdf")
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
  LR.show <- pairLR.CXCL[5, ]
  pdf("12.3.2 单个配体受体通讯网络_贝壳图_NRG.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  dev.off()
  pdf("12.3.3 单个配体受体通讯网络_弦图_NRG.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  dev.off()

## PSAP
cellchat@netP$pathways
  pathways.show <- "PSAP"
  pdf("12.2.1 重要通路通讯_贝壳图_PSAP.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  pdf("12.2.2 重要通路通讯_弦图_PSAP.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  pdf("12.2.3 通讯热图_PSAP.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
  
  pdf("12.3.1 配体受体贡献度_PSAP.pdf")
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
  LR.show <- pairLR.CXCL[1, ]
  pdf("12.3.2 单个配体受体通讯网络_贝壳图_PSAP.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  dev.off()
  pdf("12.3.3 单个配体受体通讯网络_弦图_PSAP.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  dev.off()

## GRN
cellchat@netP$pathways
  pathways.show <- "GRN"
  pdf("12.2.1 重要通路通讯_贝壳图_GRN.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  pdf("12.2.2 重要通路通讯_弦图_GRN.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  pdf("12.2.3 通讯热图_GRN.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
  
  pdf("12.3.1 配体受体贡献度_GRN.pdf")
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
  LR.show <- pairLR.CXCL[1, ]
  pdf("12.3.2 单个配体受体通讯网络_贝壳图_GRN.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  dev.off()
  pdf("12.3.3 单个配体受体通讯网络_弦图_GRN.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  dev.off()

## SPP1
cellchat@netP$pathways
  pathways.show <- "SPP1"
  pdf("12.2.1 重要通路通讯_贝壳图_SPP1.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  pdf("12.2.2 重要通路通讯_弦图_SPP1.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  pdf("12.2.3 通讯热图_SPP1.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
  
  pdf("12.3.1 配体受体贡献度_SPP1.pdf")
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
  LR.show <- pairLR.CXCL[2, ]
  pdf("12.3.2 单个配体受体通讯网络_贝壳图_SPP1.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  dev.off()
  pdf("12.3.3 单个配体受体通讯网络_弦图_SPP1.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  dev.off()

## IGF
cellchat@netP$pathways
  pathways.show <- "IGF"
  pdf("12.2.1 重要通路通讯_贝壳图_IGF.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  pdf("12.2.2 重要通路通讯_弦图_IGF.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  pdf("12.2.3 通讯热图_IGF.pdf")
  par(mfrow = c(1, 1), xpd = T)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
  
  pdf("12.3.1 配体受体贡献度_IGF.pdf")
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
  LR.show <- pairLR.CXCL[1, ]
  pdf("12.3.2 单个配体受体通讯网络_贝壳图_IGF.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  dev.off()
  pdf("12.3.3 单个配体受体通讯网络_弦图_IGF.pdf")
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  dev.off()

## Neurons
levels(cellchat@idents)
  pdf("12.4.1 多配体受体介导气泡图_Neu.pdf")
  netVisual_bubble(cellchat, sources.use = 1, remove.isolate = F)
  dev.off()
  pdf("12.4.2 多配体受体介导弦图_Neu.pdf")
  netVisual_chord_gene(cellchat, sources.use = 1, lab.cex = 0.25, legend.pos.y = 10)
  dev.off()

## Astrocytes
levels(cellchat@idents)
  pdf("12.4.1 多配体受体介导气泡图_Astro.pdf")
  netVisual_bubble(cellchat, sources.use = 5, remove.isolate = F)
  dev.off()
  pdf("12.4.2 多配体受体介导弦图_Astro.pdf")
  netVisual_chord_gene(cellchat, sources.use = 5, lab.cex = 0.25, legend.pos.y = 10)
  dev.off()

## Microglia
levels(cellchat@idents)
  pdf("12.4.1 多配体受体介导气泡图_Micro.pdf")
  netVisual_bubble(cellchat, sources.use = 2, remove.isolate = F)
  dev.off()
  pdf("12.4.2 多配体受体介导弦图_Micro.pdf")
  netVisual_chord_gene(cellchat, sources.use = 2, lab.cex = 0.25, legend.pos.y = 10)
  dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
## NRG
  cellchat@netP$pathways
  pathways.show <- "NRG"
  pdf("12.5.1 中心性得分_NRG.pdf")
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  pdf("12.5.2 发送者和接收者散点图_NRG.pdf")
  netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
  dev.off()

## PSAP
  cellchat@netP$pathways
  pathways.show <- "PSAP"
  pdf("12.5.1 中心性得分_PSAP.pdf")
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  pdf("12.5.2 发送者和接收者散点图_PSAP.pdf")
  netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
  dev.off()

## GRN
  cellchat@netP$pathways
  pathways.show <- "GRN"
  pdf("12.5.1 中心性得分_GRN.pdf")
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  pdf("12.5.2 发送者和接收者散点图_GRN.pdf")
  netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
  dev.off()

## SPP1
  cellchat@netP$pathways
  pathways.show <- "SPP1"
  pdf("12.5.1 中心性得分_SPP1.pdf")
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  pdf("12.5.2 发送者和接收者散点图_SPP1.pdf")
  netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
  dev.off()

## IGF
  cellchat@netP$pathways
  pathways.show <- "IGF"
  pdf("12.5.1 中心性得分_IGF.pdf")
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  pdf("12.5.2 发送者和接收者散点图_IGF.pdf")
  netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
  dev.off()

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("12.5.3 信号通路在细胞通讯中的作用.pdf", width = 12, height = 6)
ht1 + ht2
dev.off()

pdf("12.6.1 selectK推断_输出.pdf")
selectK(cellchat, pattern = "outgoing")
dev.off()
nPatterns = 9
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
pdf("12.6.2 通讯模式热图_输出.pdf", width = 8, height = 5)
identifyCommunicationPatterns(cellchat, pattern = "outgoing", font.size = 6, k = nPatterns)
dev.off()
pdf("12.6.3 通讯模式桑基图_输出.pdf")
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
pdf("12.6.4 通讯模式点图_输出.pdf")
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

pdf("12.7.1 selectK推断_输入.pdf")
selectK(cellchat, pattern = "incoming")
dev.off()
nPatterns = 8
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
pdf("12.7.2 通讯模式热图_输入.pdf", width = 8, height = 5)
identifyCommunicationPatterns(cellchat, pattern = "incoming", font.size = 6, k = nPatterns)
dev.off()
pdf("12.7.3 通讯模式桑基图_输入.pdf")
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
pdf("12.7.4 通讯模式点图_输入.pdf")
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional")
pdf("12.8.1 二维分类图.pdf")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
plan(sequential)   # 强制单线程
options(future.globals.maxSize = Inf)
cellchat <- netClustering(cellchat, type = "structural")
pdf("12.8.2 二维信号组图.pdf")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()
pdf("12.8.3 分组二维信号组图.pdf")
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()
qsave(cellchat, file = "cellchat.qs", nthreads = detectCores())
cellchat = qread("cellchat.qs", nthreads = detectCores())
