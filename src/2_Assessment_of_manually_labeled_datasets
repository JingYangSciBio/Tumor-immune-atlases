#####################################################################
#### 2023.05.18 knnown single-cell data from as positive test
#####################################################################

rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
"%ni%" <- Negate("%in%")

# load("/Users/jingyang/Dropbox/work/cancer-immu for biological discoveries/code/scRNA-seq_from_Qi_object.Rdata")
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
author_embedding <- read.table(paste0(workdir,"2.5_mapping_known_dataset_Bi_embedding_from_authhor.txt"),header=T,sep="\t")
author_embedding <- author_embedding[-1,]
rownames(author_embedding) <- author_embedding[,1]
author_embedding <- author_embedding[,-1]
author_embedding_subset <- author_embedding[rownames(author_embedding) %in% rownames(Embeddings(object)),]
test1 <- as.numeric(author_embedding_subset[,1])
test2 <- as.numeric(author_embedding_subset[,2])
author_embedding_subset_final <- matrix(c(test1,test2),ncol=2)
rownames(author_embedding_subset_final) <- rownames(author_embedding_subset)
object$cell_types_from_author <- object$FinalCellType
object[["umap_from_author"]] <- CreateDimReducObject(embeddings = author_embedding_subset_final,assay='RNA',key="umap_from_author_")
DimPlot(object,reduction = "umap_from_author",group.by = "cell_types_from_author")
saveRDS(object,file=paste0(workdir,"2.5_mapping_known_dataset_Bi.rds"))
object_T_NK <- subset(object,subset=cell_types_from_author %in% c("41BB-Hi CD8+ T cell","41BB-Lo CD8+ T cell","Cycling CD8+ T cell",
                                                                  "Effector T-Helper","FGFBP2- NK","FGFBP2+ NK","Memory T-Helper",
                                                                  "MX1-Hi CD8+ T cell","NKT","T-Reg"))

object_T_NK <- FindVariableFeatures(object_T_NK,selection.method = "vst",nfeatures = 2000);
var.genes <- VariableFeatures(object_T_NK);
object_T_NK <- ScaleData(object_T_NK)
object_T_NK <- RunPCA(object_T_NK, features = var.genes)
object_T_NK <- FindNeighbors(object_T_NK, dims = 1:30)
object_T_NK <- FindClusters(object_T_NK, resolution = 0.5)
object_T_NK <- RunUMAP(object_T_NK, reduction = "pca", dims = 1:30, return.model=T)
medianexp<-t(apply(GetAssayData(object_T_NK,slot="data"),1,function(x){tapply(x,object_T_NK@active.ident,mean)}))
if (nrow(medianexp)==1) medianexp<-t(medianexp)
markerfile<-"/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/PanglaoDB_markers_21_Oct_2019.tsv"
add_Subtype=FALSE
species <- "Hs"
source("/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/markerCode_filter.R")
celltype_predictmethod="cta"
predict_celltype<-ORA_celltype(medianexp,cellType,method=celltype_predictmethod)
new.cluster.ids<-paste(levels(object_T_NK@active.ident),rownames(predict_celltype$predict_result),sep="_")
names(new.cluster.ids) <- levels(object_T_NK)
object_T_NK <- RenameIdents(object_T_NK, new.cluster.ids)
object_T_NK$cellcluster_ann <- Idents(object_T_NK)
saveRDS(object_T_NK,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_20250506.rds"))
object_T_NK_markers <- FindAllMarkers(object_T_NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
saveRDS(object_T_NK_markers,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_markers_our_calculation_20250506.rds"))
Idents(object_T_NK) <- object_T_NK$cell_types_from_author
object_T_NK_markers <- FindAllMarkers(object_T_NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
saveRDS(object_T_NK_markers,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_markers_from_author_20250506.rds"))
# 
# #### check Bi dataset
Bi_T_NK_markers_from_authors <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_markers_from_author.rds"))
Bi_T_NK_markers_our_calculation <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_markers_our_calculation.rds"))
# 
unique(Bi_T_NK_markers_from_authors$cluster)
Bi_T_NK_markers_from_authors[Bi_T_NK_markers_from_authors$cluster=="41BB-Hi CD8+ T cell","gene"][1:40]
Bi_T_NK <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK.rds"))
DimPlot(Bi_T_NK,group.by="cell_types_from_author",label = T)
DimPlot(Bi_T_NK,label = T)
DimPlot(Bi_T_NK,group.by = "CellType",label = T)
Bi_T_NK_markers_our_calculation %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC) -> top10_our_calculation
DoHeatmap(Bi_T_NK,features=c(top10_our_calculation$gene,"CD4","CD8A","CD8B"),slot = "data") + NoLegend()
# 
Bi_T_NK_markers_from_authors %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC) -> top10_from_authors
DoHeatmap(Bi_T_NK,features=top10_from_authors$gene,group.by = "cell_types_from_author") + NoLegend()
# 
# #### use pbmc 10K as known dataset
# pbmc10k <- Read10X(data.dir = paste0(workdir,"2.5_pbmc_10k_v3_filtered_feature_bc_matrix/"))
# pbmc10k <- CreateSeuratObject(counts = pbmc10k, project = "pbmc10k", min.cells = 3, min.features = 200)
# pbmc10k[["percent.mt"]] <- PercentageFeatureSet(pbmc10k, pattern = "^MT-")
# VlnPlot(pbmc10k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(pbmc10k, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(pbmc10k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# pbmc10k <- NormalizeData(pbmc10k)
# pbmc10k <- FindVariableFeatures(pbmc10k, selection.method = "vst", nfeatures = 2000)
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc10k), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc10k)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# all.genes <- rownames(pbmc10k)
# pbmc10k <- ScaleData(pbmc10k, features = all.genes)
# pbmc10k <- RunPCA(pbmc10k, features = VariableFeatures(object = pbmc10k))
# pbmc10k <- FindNeighbors(pbmc10k, dims = 1:10)
# pbmc10k <- FindClusters(pbmc10k, resolution = 0.3)
# pbmc10k <- RunUMAP(pbmc10k, dims = 1:10)
# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes
# pbmc10k <- CellCycleScoring(pbmc10k, s.features = s.genes, g2m.features = g2m.genes)
# medianexp<-t(apply(GetAssayData(pbmc10k,slot="data"),1,function(x){tapply(x,pbmc10k@active.ident,mean)}))
# if (nrow(medianexp)==1) medianexp<-t(medianexp)
# markerfile<-"/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/PanglaoDB_markers_21_Oct_2019.tsv"
# add_Subtype=FALSE
# species <- "Hs"
# source("/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/markerCode_filter.R")
# celltype_predictmethod="cta"
# predict_celltype<-ORA_celltype(medianexp,cellType,method=celltype_predictmethod)
# new.cluster.ids<-paste(levels(pbmc10k@active.ident),rownames(predict_celltype$predict_result),sep="_")
# names(new.cluster.ids) <- levels(pbmc10k)
# pbmc10k <- RenameIdents(pbmc10k, new.cluster.ids)
# pbmc10k$cellcluster_ann <- Idents(pbmc10k)
# pbmc10k_markers <- FindAllMarkers(pbmc10k,only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0)
# 
# hpca.ref <- celldex::HumanPrimaryCellAtlasData()
# dice.ref <- celldex::DatabaseImmuneCellExpressionData()
# monaco.ref <- celldex::MonacoImmuneData()
# sce <- as.SingleCellExperiment(pbmc10k)
# hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
# monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
# monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
# pbmc10k$hpca.main   <- hpca.main$pruned.labels
# pbmc10k$dice.main   <- dice.main$pruned.labels
# pbmc10k$monaco.main <- monaco.main$pruned.labels
# pbmc10k$hpca.fine   <- hpca.fine$pruned.labels
# pbmc10k$dice.fine   <- dice.fine$pruned.labels
# pbmc10k$monaco.fine <- monaco.fine$pruned.labels
# 
# DimPlot(pbmc10k,group.by = "cellcluster_ann",label = T)
# DimPlot(pbmc10k,group.by = "hpca.fine",label = T)
# DimPlot(pbmc10k,group.by = "hpca.main",label = T)
# saveRDS(pbmc10k,file=paste0(workdir,"2.5_pbmc_10k.rds"))
# saveRDS(pbmc10k_markers,file=paste0(workdir,"2.5_pbmc_10k_markers.rds"))
# 
# 
# #### T and NK cells
# pbmc10k_T_NK <- subset(pbmc10k, idents = c("1_T cells","6_NK cells","2_T memory cells","3_T cells"))
# pbmc10k_T_NK <- FindVariableFeatures(pbmc10k_T_NK, selection.method = "vst", nfeatures = 2000)
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc10k_T_NK), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc10k_T_NK)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# all.genes <- rownames(pbmc10k_T_NK)
# pbmc10k_T_NK <- ScaleData(pbmc10k_T_NK, features = all.genes)
# pbmc10k_T_NK <- RunPCA(pbmc10k_T_NK, features = VariableFeatures(object = pbmc10k_T_NK))
# pbmc10k_T_NK <- FindNeighbors(pbmc10k_T_NK, dims = 1:10)
# pbmc10k_T_NK <- FindClusters(pbmc10k_T_NK, resolution = 0.3)
# pbmc10k_T_NK <- RunUMAP(pbmc10k_T_NK, dims = 1:10)
# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes
# pbmc10k_T_NK <- CellCycleScoring(pbmc10k_T_NK, s.features = s.genes, g2m.features = g2m.genes)
# medianexp<-t(apply(GetAssayData(pbmc10k_T_NK,slot="data"),1,function(x){tapply(x,pbmc10k_T_NK@active.ident,mean)}))
# if (nrow(medianexp)==1) medianexp<-t(medianexp)
# markerfile<-"/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/PanglaoDB_markers_21_Oct_2019.tsv"
# add_Subtype=FALSE
# species <- "Hs"
# source("/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/markerCode_filter.R")
# celltype_predictmethod="cta"
# predict_celltype<-ORA_celltype(medianexp,cellType,method=celltype_predictmethod)
# new.cluster.ids<-paste(levels(pbmc10k_T_NK@active.ident),rownames(predict_celltype$predict_result),sep="_")
# names(new.cluster.ids) <- levels(pbmc10k_T_NK)
# pbmc10k_T_NK <- RenameIdents(pbmc10k_T_NK, new.cluster.ids)
# pbmc10k_T_NK$cellcluster_ann <- Idents(pbmc10k_T_NK)
# pbmc10k_T_NK_markers <- FindAllMarkers(pbmc10k_T_NK,only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0)
# 
# hpca.ref <- celldex::HumanPrimaryCellAtlasData()
# dice.ref <- celldex::DatabaseImmuneCellExpressionData()
# monaco.ref <- celldex::MonacoImmuneData()
# sce <- as.SingleCellExperiment(pbmc10k_T_NK)
# hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
# monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
# monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
# pbmc10k_T_NK$hpca.main   <- hpca.main$pruned.labels
# pbmc10k_T_NK$dice.main   <- dice.main$pruned.labels
# pbmc10k_T_NK$monaco.main <- monaco.main$pruned.labels
# pbmc10k_T_NK$hpca.fine   <- hpca.fine$pruned.labels
# pbmc10k_T_NK$dice.fine   <- dice.fine$pruned.labels
# pbmc10k_T_NK$monaco.fine <- monaco.fine$pruned.labels
# 
# DimPlot(pbmc10k_T_NK,group.by = "cellcluster_ann",label = T)
# DimPlot(pbmc10k_T_NK,group.by = "hpca.fine",label = T)
# DimPlot(pbmc10k_T_NK,group.by = "hpca.main",label = T)
# saveRDS(pbmc10k_T_NK,file=paste0(workdir,"2.5_pbmc_10k_T_NK.rds"))
# saveRDS(pbmc10k_T_NK_markers,file=paste0(workdir,"2.5_pbmc_10k_T_NK_markers.rds"))
# 
# 
# ##### check pbmc10k
# DimPlot(pbmc10k,group.by = "hpca.main",label = T)




# #### use pbmc 3k as known dataset
# pbmc3k <- pbmc3k.SeuratData::pbmc3k
# pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")
# VlnPlot(pbmc3k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# pbmc3k <- NormalizeData(pbmc3k)
# pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 2000)
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc3k), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc3k)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# all.genes <- rownames(pbmc3k)
# pbmc3k <- ScaleData(pbmc3k, features = all.genes)
# pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(object = pbmc3k))
# pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)
# pbmc3k <- FindClusters(pbmc3k, resolution = 0.5)
# pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)
# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes
# pbmc3k <- CellCycleScoring(pbmc3k, s.features = s.genes, g2m.features = g2m.genes)
# medianexp<-t(apply(GetAssayData(pbmc3k,slot="data"),1,function(x){tapply(x,pbmc3k@active.ident,mean)}))
# if (nrow(medianexp)==1) medianexp<-t(medianexp)
# markerfile<-"/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/PanglaoDB_markers_21_Oct_2019.tsv"
# add_Subtype=FALSE
# species <- "Hs"
# source("/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/markerCode_filter.R")
# celltype_predictmethod="cta"
# predict_celltype<-ORA_celltype(medianexp,cellType,method=celltype_predictmethod)
# new.cluster.ids<-paste(levels(pbmc3k@active.ident),rownames(predict_celltype$predict_result),sep="_")
# names(new.cluster.ids) <- levels(pbmc3k)
# pbmc3k <- RenameIdents(pbmc3k, new.cluster.ids)
# pbmc3k$cellcluster_ann <- Idents(pbmc3k)
# pbmc3k_markers <- FindAllMarkers(pbmc3k,only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0)
# 
# hpca.ref <- celldex::HumanPrimaryCellAtlasData()
# dice.ref <- celldex::DatabaseImmuneCellExpressionData()
# monaco.ref <- celldex::MonacoImmuneData()
# sce <- as.SingleCellExperiment(pbmc3k)
# hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
# monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
# monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
# pbmc3k$hpca.main   <- hpca.main$pruned.labels
# pbmc3k$dice.main   <- dice.main$pruned.labels
# pbmc3k$monaco.main <- monaco.main$pruned.labels
# pbmc3k$hpca.fine   <- hpca.fine$pruned.labels
# pbmc3k$dice.fine   <- dice.fine$pruned.labels
# pbmc3k$monaco.fine <- monaco.fine$pruned.labels
# 
# DimPlot(pbmc3k,group.by = "cellcluster_ann",label = T)
# DimPlot(pbmc3k,group.by = "hpca.fine",label = T)
# DimPlot(pbmc3k,group.by = "hpca.main",label = T)
# saveRDS(pbmc3k,file=paste0(workdir,"2.5_pbmc_3k.rds"))
# saveRDS(pbmc3k_markers,file=paste0(workdir,"2.5_pbmc_3k_markers.rds"))
# 
# 
# #### T and NK cells
# pbmc3k_T_NK <- subset(pbmc3k, idents = c("0_T memory cells","2_T memory cells","4_T cells","6_Gamma delta T cells"))
# pbmc3k_T_NK <- FindVariableFeatures(pbmc3k_T_NK, selection.method = "vst", nfeatures = 2000)
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc3k_T_NK), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc3k_T_NK)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# all.genes <- rownames(pbmc3k_T_NK)
# pbmc3k_T_NK <- ScaleData(pbmc3k_T_NK, features = all.genes)
# pbmc3k_T_NK <- RunPCA(pbmc3k_T_NK, features = VariableFeatures(object = pbmc3k_T_NK))
# pbmc3k_T_NK <- FindNeighbors(pbmc3k_T_NK, dims = 1:10)
# pbmc3k_T_NK <- FindClusters(pbmc3k_T_NK, resolution = 0.5)
# pbmc3k_T_NK <- RunUMAP(pbmc3k_T_NK, dims = 1:10)
# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes
# pbmc3k_T_NK <- CellCycleScoring(pbmc3k_T_NK, s.features = s.genes, g2m.features = g2m.genes)
# medianexp<-t(apply(GetAssayData(pbmc3k_T_NK,slot="data"),1,function(x){tapply(x,pbmc3k_T_NK@active.ident,mean)}))
# if (nrow(medianexp)==1) medianexp<-t(medianexp)
# markerfile<-"/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/PanglaoDB_markers_21_Oct_2019.tsv"
# add_Subtype=FALSE
# species <- "Hs"
# source("/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/subscripts/markerCode_filter.R")
# celltype_predictmethod="cta"
# predict_celltype<-ORA_celltype(medianexp,cellType,method=celltype_predictmethod)
# new.cluster.ids<-paste(levels(pbmc3k_T_NK@active.ident),rownames(predict_celltype$predict_result),sep="_")
# names(new.cluster.ids) <- levels(pbmc3k_T_NK)
# pbmc3k_T_NK <- RenameIdents(pbmc3k_T_NK, new.cluster.ids)
# pbmc3k_T_NK$cellcluster_ann <- Idents(pbmc3k_T_NK)
# pbmc3k_T_NK_markers <- FindAllMarkers(pbmc3k_T_NK,only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0)
# 
# hpca.ref <- celldex::HumanPrimaryCellAtlasData()
# dice.ref <- celldex::DatabaseImmuneCellExpressionData()
# monaco.ref <- celldex::MonacoImmuneData()
# sce <- as.SingleCellExperiment(pbmc3k_T_NK)
# hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
# monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
# monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
# pbmc3k_T_NK$hpca.main   <- hpca.main$pruned.labels
# pbmc3k_T_NK$dice.main   <- dice.main$pruned.labels
# pbmc3k_T_NK$monaco.main <- monaco.main$pruned.labels
# pbmc3k_T_NK$hpca.fine   <- hpca.fine$pruned.labels
# pbmc3k_T_NK$dice.fine   <- dice.fine$pruned.labels
# pbmc3k_T_NK$monaco.fine <- monaco.fine$pruned.labels
# 
# DimPlot(pbmc3k_T_NK,group.by = "cellcluster_ann",label = T)
# DimPlot(pbmc3k_T_NK,group.by = "hpca.fine",label = T)
# DimPlot(pbmc3k_T_NK,group.by = "hpca.main",label = T)
# saveRDS(pbmc3k_T_NK,file=paste0(workdir,"2.5_pbmc_3k_T_NK.rds"))
# saveRDS(pbmc3k_T_NK_markers,file=paste0(workdir,"2.5_pbmc_3k_T_NK_markers.rds"))
# 
# 
# ##### check pbmc3k
# DimPlot(pbmc3k,group.by = "hpca.main",label = T)



# #####  mapping to atlases using pbmc10k, not Bi dataset. Because Bi dataset is not good
# pbmc10k_T_NK_cellsbubian <- readRDS(paste0(workdir,"2.5_pbmc_10k_T_NK.rds"))
# 
# #### 1. TICA
# 
# atlas_TICA_T_NK_cells <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells.rds"))
# pbmc10k_T_NK_cells <- pbmc10k_T_NK_cellsbubian
# anchors_map_pbmc10k_T_NK_cells_to_TICA_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_TICA_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc10k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc10k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc10k_T_NK_cells_to_TICA_T_NK_cells,
#   query = pbmc10k_T_NK_cells,
#   reference = atlas_TICA_T_NK_cells,
#   refdata = list(
#     to_TICA_T_NK_cells = "cell_types_from_author",
#     to_TICA_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_TICA_T_NK_cells = "integrated"#"RNA"
#   ),
#   new.reduction.name = "to_TICA_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc10k_T_NK_cells_to_TICA_T_NK_cells")
# rm("atlas_TICA_T_NK_cells")
# pbmc10k_T_NK_cells[["to_TICA_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc10k_T_NK_cells,reduction = "ref.umap"),
#                                                                         key = "to_TICA_T_NK_cells_",
#                                                                         assay = DefaultAssay(pbmc10k_T_NK_cells))
# saveRDS(pbmc10k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_TICA_T_NK_cells_integrated.rds"))
# rm("pbmc10k_T_NK_cells")
# 
# #### pancancer T cells CD8 CD4
# pbmc10k_T_NK_cells <- pbmc10k_T_NK_cellsbubian
# atlas_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_filled.rds"))
# anchors_map_pbmc10k_T_NK_cells_to_pancancer_T_cells_CD8_CD4 <- FindTransferAnchors(
#   reference = atlas_pancancer_T_cells_CD8_CD4,#atlas_pancancer_T_cells_CD8,
#   query = pbmc10k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc10k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc10k_T_NK_cells_to_pancancer_T_cells_CD8_CD4,
#   query = pbmc10k_T_NK_cells,
#   reference = atlas_pancancer_T_cells_CD8_CD4,
#   refdata = list(
#     to_pancancer_T_cells_CD8_CD4 = "cell_types_from_author", #"cell_types_from_author_simplized"
#     to_pancancer_T_cells_CD8_CD4_our_calculation = "cellcluster_ann",
#     predicted_to_pancancer_T_cells_CD8_CD4 = "RNA"
#   ),
#   new.reduction.name = "to_pancancer_T_cells_CD8_CD4_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc10k_T_NK_cells_to_pancancer_T_cells_CD8_CD4")
# rm("atlas_pancancer_T_cells_CD8_CD4")
# pbmc10k_T_NK_cells[["to_pancancer_T_cells_CD8_CD4_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc10k_T_NK_cells,reduction = "ref.umap"),
#                                                                                   key = "to_pancancer_T_cells_CD8_CD4_", 
#                                                                                   assay = DefaultAssay(pbmc10k_T_NK_cells))
# saveRDS(pbmc10k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# rm("pbmc10k_T_NK_cells")
# 
# # 3. map to pancancer_blueprint_T_NK_cells
# atlas_pancancer_blueprint_T_NK_cells <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells.rds"))
# DefaultAssay(atlas_pancancer_blueprint_T_NK_cells) <- "integrated"
# pbmc10k_T_NK_cells <- pbmc10k_T_NK_cellsbubian
# anchors_map_pbmc10k_T_NK_cells_to_pancancer_blueprint_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_pancancer_blueprint_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc10k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc10k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc10k_T_NK_cells_to_pancancer_blueprint_T_NK_cells,
#   query = pbmc10k_T_NK_cells,
#   reference = atlas_pancancer_blueprint_T_NK_cells,
#   refdata = list(
#     to_pancancer_blueprint_T_NK_cells = "cell_types_from_author",
#     to_pancancer_blueprint_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_pancancer_blueprint_T_NK_cells = "integrated"#"RNA"
#   ),
#   new.reduction.name = "to_pancancer_blueprint_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc10k_T_NK_cells_to_pancancer_blueprint_T_NK_cells")
# rm("atlas_pancancer_blueprint_T_NK_cells")
# pbmc10k_T_NK_cells[["to_pancancer_blueprint_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc10k_T_NK_cells,reduction = "ref.umap"),
#                                                                                        key = "to_pancancer_blueprint_T_NK_cells_",
#                                                                                        assay = DefaultAssay(pbmc10k_T_NK_cells))
# saveRDS(pbmc10k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated.rds"))
# rm("pbmc10k_T_NK_cells")
# 
# # 4. map to ABTC_T_NK_cells
# pbmc10k_T_NK_cells <- pbmc10k_T_NK_cellsbubian
# atlas_ABTC_T_NK_cells <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells.rds"))
# anchors_map_pbmc10k_T_NK_cells_to_ABTC_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_ABTC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc10k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc10k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc10k_T_NK_cells_to_ABTC_T_NK_cells,
#   query = pbmc10k_T_NK_cells,
#   reference = atlas_ABTC_T_NK_cells,
#   refdata = list(
#     to_ABTC_T_NK_cells = "cell_types_from_author",
#     to_ABTC_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_ABTC_T_NK_cells = "RNA"
#   ),
#   new.reduction.name = "to_ABTC_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc10k_T_NK_cells_to_ABTC_T_NK_cells")
# rm("atlas_ABTC_T_NK_cells")
# pbmc10k_T_NK_cells[["to_ABTC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc10k_T_NK_cells,reduction = "ref.umap"),
#                                                                         key = "to_ABTC_T_NK_cells_",
#                                                                         assay = DefaultAssay(pbmc10k_T_NK_cells))
# saveRDS(pbmc10k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# rm("pbmc10k_T_NK_cells")
# 
# # 5. map to HCC_T_NK_cells
# pbmc10k_T_NK_cells <- pbmc10k_T_NK_cellsbubian
# atlas_HCC_T_NK_cells <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells.rds"))
# DefaultAssay(atlas_HCC_T_NK_cells) <- "integrated"
# anchors_map_pbmc10k_T_NK_cells_to_HCC_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_HCC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc10k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc10k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc10k_T_NK_cells_to_HCC_T_NK_cells,
#   query = pbmc10k_T_NK_cells,
#   reference = atlas_HCC_T_NK_cells,
#   refdata = list(
#     to_HCC_T_NK_cells = "cell_types_from_author",
#     to_HCC_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_HCC_T_NK_cells = "RNA"
#   ),
#   new.reduction.name = "to_HCC_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc10k_T_NK_cells_to_HCC_T_NK_cells")
# rm("atlas_HCC_T_NK_cells")
# pbmc10k_T_NK_cells[["to_HCC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc10k_T_NK_cells,reduction = "ref.umap"),
#                                                                        key = "to_HCC_T_NK_cells_",
#                                                                        assay = DefaultAssay(pbmc10k_T_NK_cells))
# saveRDS(pbmc10k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_HCC_T_NK_cells.rds"))
# rm("pbmc10k_T_NK_cells")
# 
# # 6. map to STAD_T_NK_cells
# pbmc10k_T_NK_cells <- pbmc10k_T_NK_cellsbubian
# atlas_STAD_T_NK_cells <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells.rds"))
# anchors_map_pbmc10k_T_NK_cells_to_STAD_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_STAD_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc10k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc10k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc10k_T_NK_cells_to_STAD_T_NK_cells,
#   query = pbmc10k_T_NK_cells,
#   reference = atlas_STAD_T_NK_cells,
#   refdata = list(
#     to_STAD_T_NK_cells = "cell_types_from_author",
#     to_STAD_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_STAD_T_NK_cells = "RNA"
#   ),
#   new.reduction.name = "to_STAD_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc10k_T_NK_cells_to_STAD_T_NK_cells")
# rm("atlas_STAD_T_NK_cells")
# pbmc10k_T_NK_cells[["to_STAD_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc10k_T_NK_cells,reduction = "ref.umap"),
#                                                                         key = "to_STAD_T_NK_cells_",
#                                                                         assay = DefaultAssay(pbmc10k_T_NK_cells))
# saveRDS(pbmc10k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_STAD_T_NK_cells.rds"))
# rm("pbmc10k_T_NK_cells")


#####  mapping to atlases using pbmc3k, not Bi dataset. Because Bi dataset is not good
# pbmc3k_T_NK_cellsbubian <- readRDS(paste0(workdir,"2.5_pbmc_3k_T_NK.rds"))
# 
# #### 1. TICA
# 
# atlas_TICA_T_NK_cells <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells.rds"))
# pbmc3k_T_NK_cells <- pbmc3k_T_NK_cellsbubian
# anchors_map_pbmc3k_T_NK_cells_to_TICA_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_TICA_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc3k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc3k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc3k_T_NK_cells_to_TICA_T_NK_cells,
#   query = pbmc3k_T_NK_cells,
#   reference = atlas_TICA_T_NK_cells,
#   refdata = list(
#     to_TICA_T_NK_cells = "cell_types_from_author",
#     to_TICA_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_TICA_T_NK_cells = "integrated"#"RNA"
#   ),
#   new.reduction.name = "to_TICA_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc3k_T_NK_cells_to_TICA_T_NK_cells")
# rm("atlas_TICA_T_NK_cells")
# pbmc3k_T_NK_cells[["to_TICA_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc3k_T_NK_cells,reduction = "ref.umap"),
#                                                                        key = "to_TICA_T_NK_cells_",
#                                                                        assay = DefaultAssay(pbmc3k_T_NK_cells))
# saveRDS(pbmc3k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_TICA_T_NK_cells_integrated.rds"))
# rm("pbmc3k_T_NK_cells")
# 
# #### pancancer T cells CD8 CD4
# pbmc3k_T_NK_cells <- pbmc3k_T_NK_cellsbubian
# atlas_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_filled.rds"))
# anchors_map_pbmc3k_T_NK_cells_to_pancancer_T_cells_CD8_CD4 <- FindTransferAnchors(
#   reference = atlas_pancancer_T_cells_CD8_CD4,#atlas_pancancer_T_cells_CD8,
#   query = pbmc3k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc3k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc3k_T_NK_cells_to_pancancer_T_cells_CD8_CD4,
#   query = pbmc3k_T_NK_cells,
#   reference = atlas_pancancer_T_cells_CD8_CD4,
#   refdata = list(
#     to_pancancer_T_cells_CD8_CD4 = "cell_types_from_author", #"cell_types_from_author_simplized"
#     to_pancancer_T_cells_CD8_CD4_our_calculation = "cellcluster_ann",
#     predicted_to_pancancer_T_cells_CD8_CD4 = "RNA"
#   ),
#   new.reduction.name = "to_pancancer_T_cells_CD8_CD4_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc3k_T_NK_cells_to_pancancer_T_cells_CD8_CD4")
# rm("atlas_pancancer_T_cells_CD8_CD4")
# pbmc3k_T_NK_cells[["to_pancancer_T_cells_CD8_CD4_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc3k_T_NK_cells,reduction = "ref.umap"),
#                                                                                  key = "to_pancancer_T_cells_CD8_CD4_", 
#                                                                                  assay = DefaultAssay(pbmc3k_T_NK_cells))
# saveRDS(pbmc3k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# rm("pbmc3k_T_NK_cells")
# 
# # 3. map to pancancer_blueprint_T_NK_cells
# atlas_pancancer_blueprint_T_NK_cells <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells.rds"))
# DefaultAssay(atlas_pancancer_blueprint_T_NK_cells) <- "integrated"
# pbmc3k_T_NK_cells <- pbmc3k_T_NK_cellsbubian
# anchors_map_pbmc3k_T_NK_cells_to_pancancer_blueprint_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_pancancer_blueprint_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc3k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc3k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc3k_T_NK_cells_to_pancancer_blueprint_T_NK_cells,
#   query = pbmc3k_T_NK_cells,
#   reference = atlas_pancancer_blueprint_T_NK_cells,
#   refdata = list(
#     to_pancancer_blueprint_T_NK_cells = "cell_types_from_author",
#     to_pancancer_blueprint_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_pancancer_blueprint_T_NK_cells = "integrated"#"RNA"
#   ),
#   new.reduction.name = "to_pancancer_blueprint_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc3k_T_NK_cells_to_pancancer_blueprint_T_NK_cells")
# rm("atlas_pancancer_blueprint_T_NK_cells")
# pbmc3k_T_NK_cells[["to_pancancer_blueprint_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc3k_T_NK_cells,reduction = "ref.umap"),
#                                                                                       key = "to_pancancer_blueprint_T_NK_cells_",
#                                                                                       assay = DefaultAssay(pbmc3k_T_NK_cells))
# saveRDS(pbmc3k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated.rds"))
# rm("pbmc3k_T_NK_cells")
# 
# # 4. map to ABTC_T_NK_cells
# pbmc3k_T_NK_cells <- pbmc3k_T_NK_cellsbubian
# atlas_ABTC_T_NK_cells <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells.rds"))
# anchors_map_pbmc3k_T_NK_cells_to_ABTC_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_ABTC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc3k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc3k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc3k_T_NK_cells_to_ABTC_T_NK_cells,
#   query = pbmc3k_T_NK_cells,
#   reference = atlas_ABTC_T_NK_cells,
#   refdata = list(
#     to_ABTC_T_NK_cells = "cell_types_from_author",
#     to_ABTC_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_ABTC_T_NK_cells = "RNA"
#   ),
#   new.reduction.name = "to_ABTC_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc3k_T_NK_cells_to_ABTC_T_NK_cells")
# rm("atlas_ABTC_T_NK_cells")
# pbmc3k_T_NK_cells[["to_ABTC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc3k_T_NK_cells,reduction = "ref.umap"),
#                                                                        key = "to_ABTC_T_NK_cells_",
#                                                                        assay = DefaultAssay(pbmc3k_T_NK_cells))
# saveRDS(pbmc3k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# rm("pbmc3k_T_NK_cells")
# 
# # 5. map to HCC_T_NK_cells
# pbmc3k_T_NK_cells <- pbmc3k_T_NK_cellsbubian
# atlas_HCC_T_NK_cells <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells.rds"))
# DefaultAssay(atlas_HCC_T_NK_cells) <- "integrated"
# anchors_map_pbmc3k_T_NK_cells_to_HCC_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_HCC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc3k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc3k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc3k_T_NK_cells_to_HCC_T_NK_cells,
#   query = pbmc3k_T_NK_cells,
#   reference = atlas_HCC_T_NK_cells,
#   refdata = list(
#     to_HCC_T_NK_cells = "cell_types_from_author",
#     to_HCC_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_HCC_T_NK_cells = "RNA"
#   ),
#   new.reduction.name = "to_HCC_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc3k_T_NK_cells_to_HCC_T_NK_cells")
# rm("atlas_HCC_T_NK_cells")
# pbmc3k_T_NK_cells[["to_HCC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc3k_T_NK_cells,reduction = "ref.umap"),
#                                                                       key = "to_HCC_T_NK_cells_",
#                                                                       assay = DefaultAssay(pbmc3k_T_NK_cells))
# saveRDS(pbmc3k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_HCC_T_NK_cells.rds"))
# rm("pbmc3k_T_NK_cells")
# 
# # 6. map to STAD_T_NK_cells
# pbmc3k_T_NK_cells <- pbmc3k_T_NK_cellsbubian
# atlas_STAD_T_NK_cells <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells.rds"))
# anchors_map_pbmc3k_T_NK_cells_to_STAD_T_NK_cells <- FindTransferAnchors(
#   reference = atlas_STAD_T_NK_cells,#atlas_pancancer_T_cells_CD8,
#   query = pbmc3k_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
#   # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
#   # normalization.method = "LogNormalize",
#   reference.reduction = "pca",
#   # k.filter = NA,
#   dims = 1:30 #reference dim
# )
# pbmc3k_T_NK_cells <- MapQuery(
#   anchorset = anchors_map_pbmc3k_T_NK_cells_to_STAD_T_NK_cells,
#   query = pbmc3k_T_NK_cells,
#   reference = atlas_STAD_T_NK_cells,
#   refdata = list(
#     to_STAD_T_NK_cells = "cell_types_from_author",
#     to_STAD_T_NK_cells_our_calculation = "cellcluster_ann",
#     predicted_to_STAD_T_NK_cells = "RNA"
#   ),
#   new.reduction.name = "to_STAD_T_NK_cells_pca",
#   reference.reduction = "pca",
#   reduction.model = "umap"
# )
# rm("anchors_map_pbmc3k_T_NK_cells_to_STAD_T_NK_cells")
# rm("atlas_STAD_T_NK_cells")
# pbmc3k_T_NK_cells[["to_STAD_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(pbmc3k_T_NK_cells,reduction = "ref.umap"),
#                                                                        key = "to_STAD_T_NK_cells_",
#                                                                        assay = DefaultAssay(pbmc3k_T_NK_cells))
# saveRDS(pbmc3k_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_STAD_T_NK_cells.rds"))
# rm("pbmc3k_T_NK_cells")



#######################################
#### 2023.05.24 tidy up data
#######################################

#### pbmc10k
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

pbmc10k_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_TICA_T_NK_cells.rds"))
cell_summary_TICA <- data.frame(pbmc10k_T_NK_cells_TICA$cellcluster_ann,paste0("TICA: ",pbmc10k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells))
rm("pbmc10k_T_NK_cells_TICA")
pbmc10k_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
cell_summary_panT <- data.frame(pbmc10k_T_NK_cells_panT$cellcluster_ann,paste0("pan- T: ",pbmc10k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4))
rm("pbmc10k_T_NK_cells_panT")
pbmc10k_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
cell_summary_panblueprint <- data.frame(pbmc10k_T_NK_cells_panblueprint$cellcluster_ann,paste0("pan- blueprint: ",pbmc10k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))
rm("pbmc10k_T_NK_cells_panblueprint")
pbmc10k_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_ABTC_T_NK_cells.rds"))
cell_summary_ABTC <- data.frame(pbmc10k_T_NK_cells_ABTC$cellcluster_ann,paste0("ABTC: ",pbmc10k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells))
rm("pbmc10k_T_NK_cells_ABTC")
pbmc10k_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_HCC_T_NK_cells.rds"))
cell_summary_HCC <- data.frame(pbmc10k_T_NK_cells_HCC$cellcluster_ann,paste0("HCC: ",pbmc10k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells))
rm("pbmc10k_T_NK_cells_HCC")
pbmc10k_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc10k_T_NK_cells_to_STAD_T_NK_cells.rds"))
cell_summary_STAD <- data.frame(pbmc10k_T_NK_cells_STAD$cellcluster_ann,paste0("STAD: ",pbmc10k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells),pbmc10k_T_NK_cells_STAD$hpca.fine)
rm("pbmc10k_T_NK_cells_STAD")

cell_summary <- cbind(cell_summary_TICA,to_panT=cell_summary_panT[,2],
                      to_panblueprint=cell_summary_panblueprint[,2],to_ABTC=cell_summary_ABTC[,2],to_HCC=cell_summary_HCC[,2],
                      to_STAD=cell_summary_STAD[,2],hpca_fine=cell_summary_STAD[,3])
colnames(cell_summary)[1:2] <-c("cellcluster_ann","to_TICA")
cell_summary[1:3,]
cell_summary[,"cellcluster_ann"] <- as.character(cell_summary[,"cellcluster_ann"])
cell_summary[,"hpca_fine"] <- as.character(cell_summary[,"hpca_fine"])


add_num_annotation <- c(unique(cell_summary$to_TICA),unique(cell_summary$to_panT),unique(cell_summary$to_panblueprint),
                        unique(cell_summary$to_ABTC),unique(cell_summary$to_HCC),unique(cell_summary$to_STAD))
add_num_annotation <- as.character(new_cluster_names[add_num_annotation])
test <- 1:49
names(test) <- add_num_annotation

index <- c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD")
# for T_cell:CD4+_effector_memory
output <- NULL

for (i in 1:length(index)) {
    # cell_summary_subset <- subset(cell_summary,hpca_fine %in% c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
    # cell_summary_subset$hpca_fine <- "CD4+ T"
    cell_summary_subset <- subset(cell_summary,hpca_fine %in% c("T_cell:CD8+","T_cell:CD8+_Central_memory","T_cell:CD8+_effector_memory","T_cell:CD8+_effector_memory_RA"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
    cell_summary_subset$hpca_fine <- "CD8+ T"#"CD4+ T"
    # cell_summary_subset <- subset(cell_summary,hpca_fine %in% c("NK_cell","NK_cell:IL2"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
    # cell_summary_subset$hpca_fine <- "NK"#"CD4+ T"
  output_pre <- as.data.frame(table(cell_summary_subset[,c(index[i],"hpca_fine")]))
  output_pre <- data.frame(group=index[i],output_pre)
  output_pre[,2] <- as.character(output_pre[,2])
  output_pre$Freq <- as.numeric(output_pre$Freq)
  # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
  colnames(output_pre)[2] <- "fill"
  output_pre$fill <- as.character(new_cluster_names[output_pre$fill])
  output_pre <- data.frame(output_pre,populations=output_pre$fill)
  output_pre$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre$fill),": ",2)))[,2])
  output_pre$Freq <- as.numeric(output_pre$Freq)
  output_pre <- output_pre[order(output_pre$Freq),]
  output_pre <- output_pre %>%
    group_by(group) %>%
    mutate(label_y = cumsum(Freq)-0.5*Freq)
  output_pre <- output_pre[order(-output_pre$Freq),]
  output <- rbind(output,output_pre)
}

output <- as.data.frame(output)
# output <- rbind(output,c("RF","CD4+ T","CD4+ T",sum(output_pre$Freq),"RF: CD4+ T",sum(output_pre$Freq)/2))
output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
# output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2)) #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
output$fill <- factor(output$fill,levels = c(output$fill))
output$populations <- factor(output$populations,levels = c(output$populations))
output$group <- factor(output$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
output$Freq <- as.numeric(output$Freq)
output$label_y <- as.numeric(output$label_y)


cols <- rainbow(length(c(unique(cell_summary$to_TICA),unique(cell_summary$to_panT),unique(cell_summary$to_panblueprint),unique(cell_summary$to_ABTC),unique(cell_summary$to_HCC),unique(cell_summary$to_STAD))))
names(cols) <- as.character(new_cluster_names[c(unique(cell_summary$to_TICA),unique(cell_summary$to_panT),unique(cell_summary$to_panblueprint),unique(cell_summary$to_ABTC),unique(cell_summary$to_HCC),unique(cell_summary$to_STAD))])
cols <- c("RF: CD4+ T" = "#000000","RF: CD8+ T" = "#000000","RF: NK" = "#000000",cols)

output$fill <- as.character(test[as.character(output$populations)])
output[output$Freq < 30,"fill"] <- "" #10

ggplot(output,aes(x=group,y=Freq,fill=populations)) + 
  geom_bar(stat="identity",colour="black") + 
  scale_fill_manual(values=cols) +
  geom_text(aes(y=label_y,label=fill),size =8) +
  ylab("cell number") + xlab ("") +
  theme(axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18,color = "black"),
        axis.text.y = element_text(size = 20,color = "black"))#,
        legend.text=element_text(size=20)) #+
  # guides(color = guide_legend(override.aes = list(size=14), ncol=1) )#+ #labs(y="")
  theme_classic()



# output <- as.data.frame(output)
#   # output <- rbind(output,c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   output <- rbind(output,c("original","CD8+ T","CD8+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2)) #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output,c("original","NK cell","NK",sum(output_pre$Freq),sum(output_pre$Freq)/2)) #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
# output$fill <- factor(output$fill,levels = c(output$fill))
# output$group <- factor(output$group,levels = c("original","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
# output$Freq <- as.numeric(output$Freq)
# output$label_y <- as.numeric(output$label_y)
# 
# ggplot(output,aes(x=group,y=Freq,fill=fill)) + 
#   geom_bar(stat="identity",colour="black") + 
#   geom_text(aes(y=label_y,label=fill))



#### pbmc3k
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

pbmc3k_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_TICA_T_NK_cells.rds"))
cell_summary_TICA <- data.frame(pbmc3k_T_NK_cells_TICA$cell_types_from_pbmc3k,paste0("TICA: ",pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells))
rm("pbmc3k_T_NK_cells_TICA")
pbmc3k_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
cell_summary_panT <- data.frame(pbmc3k_T_NK_cells_panT$cell_types_from_pbmc3k,paste0("pan- T: ",pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4))
rm("pbmc3k_T_NK_cells_panT")
pbmc3k_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
cell_summary_panblueprint <- data.frame(pbmc3k_T_NK_cells_panblueprint$cell_types_from_pbmc3k,paste0("pan- blueprint: ",pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))
rm("pbmc3k_T_NK_cells_panblueprint")
pbmc3k_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_ABTC_T_NK_cells.rds"))
cell_summary_ABTC <- data.frame(pbmc3k_T_NK_cells_ABTC$cell_types_from_pbmc3k,paste0("ABTC: ",pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells))
rm("pbmc3k_T_NK_cells_ABTC")
pbmc3k_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_HCC_T_NK_cells.rds"))
cell_summary_HCC <- data.frame(pbmc3k_T_NK_cells_HCC$cell_types_from_pbmc3k,paste0("HCC: ",pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells))
rm("pbmc3k_T_NK_cells_HCC")
pbmc3k_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_STAD_T_NK_cells.rds"))
cell_summary_STAD <- data.frame(pbmc3k_T_NK_cells_STAD$cell_types_from_pbmc3k,paste0("STAD: ",pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells),pbmc3k_T_NK_cells_STAD$hpca.fine)
rm("pbmc3k_T_NK_cells_STAD")

cell_summary <- cbind(cell_summary_TICA,to_panT=cell_summary_panT[,2],
                      to_panblueprint=cell_summary_panblueprint[,2],to_ABTC=cell_summary_ABTC[,2],to_HCC=cell_summary_HCC[,2],
                      to_STAD=cell_summary_STAD[,2],hpca_fine=cell_summary_STAD[,3])
colnames(cell_summary)[1:2] <-c("cell_types_from_pbmc3k","to_TICA")
cell_summary[1:3,]
cell_summary[,"cell_types_from_pbmc3k"] <- as.character(cell_summary[,"cell_types_from_pbmc3k"])
cell_summary[,"hpca_fine"] <- as.character(cell_summary[,"hpca_fine"])

unique(cell_summary$cell_types_from_pbmc3k)
# [1] "0_T memory cells"      "2_T memory cells"      "6_Gamma delta T cells" "4_T cells"     

index <- c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD")
# for T_cell:CD4+_effector_memory
output <- NULL
for (i in 1:length(index)) {
  # cell_summary_subset <- subset(cell_summary,cell_types_from_pbmc3k %in% c("0_T memory cells","2_T memory cells")) #
  # cell_summary_subset$hpca_fine <- "CD4+ T"
  # cell_summary_subset <- subset(cell_summary,cell_types_from_pbmc3k %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
  # cell_summary_subset$hpca_fine <- "CD8+ T"#"CD4+ T"
  cell_summary_subset <- subset(cell_summary,cell_types_from_pbmc3k %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
  cell_summary_subset$hpca_fine <- "NK"#"CD4+ T"
  output_pre <- as.data.frame(table(cell_summary_subset[,c(index[i],"hpca_fine")]))
  output_pre <- data.frame(group=index[i],output_pre)
  output_pre[,2] <- as.character(output_pre[,2])
  output_pre$Freq <- as.numeric(output_pre$Freq)
  # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
  colnames(output_pre)[2] <- "fill"
  output_pre$fill <- as.character(new_cluster_names[output_pre$fill])
  output_pre <- data.frame(output_pre,populations=output_pre$fill)
  output_pre$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre$fill),": ",2)))[,2])
  output_pre$Freq <- as.numeric(output_pre$Freq)
  output_pre <- output_pre[order(output_pre$Freq),]
  output_pre <- output_pre %>%
    group_by(group) %>%
    mutate(label_y = cumsum(Freq)-0.5*Freq)
  output_pre <- output_pre[order(-output_pre$Freq),]
  output <- rbind(output,output_pre)
}
output <- as.data.frame(output)
# output <- rbind(output,c("RF","CD4+ T","CD4+ T",sum(output_pre$Freq),"RF: CD4+ T",sum(output_pre$Freq)/2))
# output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2)) #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
output$fill <- factor(output$fill,levels = c(output$fill))
output$populations <- factor(output$populations,levels = c(output$populations))
output$group <- factor(output$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
output$Freq <- as.numeric(output$Freq)
output$label_y <- as.numeric(output$label_y)
output[output$Freq < 3,"fill"] <- "" #10

cols <- rainbow(length(c(unique(cell_summary$to_TICA),unique(cell_summary$to_panT),unique(cell_summary$to_panblueprint),unique(cell_summary$to_ABTC),unique(cell_summary$to_HCC),unique(cell_summary$to_STAD))))
names(cols) <- as.character(new_cluster_names[c(unique(cell_summary$to_TICA),unique(cell_summary$to_panT),unique(cell_summary$to_panblueprint),unique(cell_summary$to_ABTC),unique(cell_summary$to_HCC),unique(cell_summary$to_STAD))])
cols <- c("RF: CD4+ T" = "#000000","RF: CD8+ T" = "#000000","RF: NK" = "#000000",cols)

ggplot(output,aes(x=group,y=Freq,fill=populations)) + 
  geom_bar(stat="identity",colour="black") + 
  scale_fill_manual(values=cols) +
  # geom_text(aes(y=label_y,label=fill),size =3) +
  theme_classic()



###########################################
#### 2023.06.22 new plots for pbmc3k
###########################################
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
library(scales)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
new_cluster_names <- c("TICA: Recently activated CD4 T cells", "TICA: Naive-memory CD4 T cells", "TICA: Transitional memory CD4 T cells","TICA: Cytotoxic CD8 T cells",
                       "TICA: Effector memory CD8 T cells", "TICA: Th17 cells", "TICA: NK", "TICA: Terminally exhausted CD8 T cells",
                       "TICA: Naive T cells", "TICA: Regulatory T cells", "TICA: T helper cells", "TICA: Proliferative T cells",
                       "TICA: Pre-exhausted CD8 T cells",
                       
                       
                       "pan- T: CD8.c07.Temra.CX3CR1", "pan- T: CD8.c09.Tk.KIR2DL4", "pan- T: CD8.c05.Tem.CXCR5" ,
                       "pan- T: CD8.c15.ISG.IFIT1", "pan- T: CD8.c01.Tn.MAL", "pan- T: CD8.c17.Tm.NME1", "pan- T: CD8.c08.Tk.TYROBP" ,
                       "pan- T: CD8.c12.Tex.CXCL13", "pan- T: CD8.c02.Tm.IL7R", "pan- T: CD8.c10.Trm.ZNF683", "pan- T: CD8.c03.Tm.RPS12" ,
                       "pan- T: CD8.c14.Tex.TCF7", "pan- T: CD8.c06.Tem.GZMK", "pan- T: CD8.c11.Tex.PDCD1", "pan- T: CD8.c16.MAIT.SLC4A10" ,
                       "pan- T: CD8.c04.Tm.CD52", "pan- T: CD8.c13.Tex.myl12a", "pan- T: CD4.c20.Treg.TNFRSF9", "pan- T: CD4.c18.Treg.RTKN2" ,
                       "pan- T: CD4.c06.Tm.ANXA1", "pan- T: CD4.c22.ISG.IFIT1", "pan- T: CD4.c01.Tn.TCF7", "pan- T: CD4.c12.Tem.GZMK" ,
                       "pan- T: CD4.c11.Tm.GZMA", "pan- T: CD4.c13.Temra.CX3CR1", "pan- T: CD4.c16.Tfh.CXCR5", "pan- T: CD4.c07.Tm.ANXA2" ,
                       "pan- T: CD4.c17.TfhTh1.CXCL13", "pan- T: CD4.c19.Treg.S1PR1", "pan- T: CD4.c10.Tm.CAPG", "pan- T: CD4.c24.Mix.NME2" ,
                       "pan- T: CD4.c21.Treg.OAS1", "pan- T: CD4.c14.Th17.SLC4A10", "pan- T: CD4.c09.Tm.CCL5", "pan- T: CD4.c02.Tn.PASK" ,
                       "pan- T: CD4.c03.Tn.ADSL", "pan- T: CD4.c05.Tm.TNF", "pan- T: CD4.c08.Tm.CREM", "pan- T: CD4.c04.Tn.il7r" ,
                       "pan- T: CD4.c15.Th17.IL23R", "pan- T: CD4.c23.Mix.NME1",
                       
                       
                       "pan- blueprint: CD8 memory T", "pan- blueprint: CD8 exhausted cytotoxic T",
                       "pan- blueprint: CD8 pre-effector T",   "pan- blueprint: CD4 memory/effector T",
                       "pan- blueprint: NK_XCL1",   "pan- blueprint: CD4 regulatory T",
                       "pan- blueprint: CD4 naive T",   "pan- blueprint: CD8 effector T",
                       "pan- blueprint: CD4 exhausted effector T", "pan- blueprint: Cytotoxic NK",
                       
                       "ABTC: Naive CD4 T cell",     "ABTC: Exhausted CD8 T cell",
                       "ABTC: Memory CD8 T cell",    "ABTC: Treg" ,               
                       "ABTC: DNT",                  "ABTC: Effector CD8 T cell" ,
                       "ABTC: NCAM1nFCGR3Ap NK",     "ABTC: Naive CD8 T cell" ,   
                       "ABTC: NKT and NK cell",      "ABTC: NCAM1pFCGR3Ap NK"   ,"ABTC: NCAM1pFCGR3An NK",
                       
                       "HCC: 13_CD8_effector memory T", "HCC: 2_NK", "HCC: 2_CD8_MAIT",  "HCC: 20_CD4_naive T",
                       "HCC: 9_CD4_regulatory T",  "HCC: 8_NK",  "HCC: 35_CD4_central memory T", "HCC: 18_NK",
                       "HCC: 11_CD8_inter-states T", "HCC: 14_CD8_TR_memory T", "HCC: 32_CD4_other T", "HCC: 1_CD8_cytotoxic T" ,
                       "HCC: 28_NK", "HCC: 51_CD8_other T",
                       
                       "STAD: LNK1", "STAD: LNK3", "STAD: LNK4", "STAD: CD8 naive T" ,
                       "STAD: CD4 naive T cell",  "STAD: LNK2", "STAD: CD4 helper T",  "STAD: Regulatory T" ,
                       "STAD: CD8 effector T cell",  "STAD: Proliferative T")
# names(new_cluster_names) <- c("TICA: Recently.activated.CD4.T.cells", "TICA: Naive.memory.CD4.T.cells", "TICA: Transitional.memory.CD4.T.cells", "TICA: Cytotoxic.CD8.T.cells",
#                               "TICA: Effector.memory.CD8.T.cells", "TICA: Th17.cells", "TICA: NK", "TICA: Terminally.exhausted.CD8.T.cells",
#                               "TICA: Naive.T.cells", "TICA: Regulatory.T.cells", "TICA: T.helper.cells", "TICA: Proliferative.T.cells" ,
#                               "TICA: Pre.exhausted.CD8.T.cells",
#                               
#                               "pan- T: CD8.c07.Temra.CX3CR1", "pan- T: CD8.c09.Tk.KIR2DL4", "pan- T: CD8.c05.Tem.CXCR5" ,
#                               "pan- T: CD8.c15.ISG.IFIT1", "pan- T: CD8.c01.Tn.MAL", "pan- T: CD8.c17.Tm.NME1", "pan- T: CD8.c08.Tk.TYROBP" ,
#                               "pan- T: CD8.c12.Tex.CXCL13", "pan- T: CD8.c02.Tm.IL7R", "pan- T: CD8.c10.Trm.ZNF683", "pan- T: CD8.c03.Tm.RPS12" ,
#                               "pan- T: CD8.c14.Tex.TCF7", "pan- T: CD8.c06.Tem.GZMK", "pan- T: CD8.c11.Tex.PDCD1", "pan- T: CD8.c16.MAIT.SLC4A10" ,
#                               "pan- T: CD8.c04.Tm.CD52", "pan- T: CD8.c13.Tex.myl12a", "pan- T: CD4.c20.Treg.TNFRSF9", "pan- T: CD4.c18.Treg.RTKN2" ,
#                               "pan- T: CD4.c06.Tm.ANXA1", "pan- T: CD4.c22.ISG.IFIT1", "pan- T: CD4.c01.Tn.TCF7", "pan- T: CD4.c12.Tem.GZMK" ,
#                               "pan- T: CD4.c11.Tm.GZMA", "pan- T: CD4.c13.Temra.CX3CR1", "pan- T: CD4.c16.Tfh.CXCR5", "pan- T: CD4.c07.Tm.ANXA2" ,
#                               "pan- T: CD4.c17.TfhTh1.CXCL13", "pan- T: CD4.c19.Treg.S1PR1", "pan- T: CD4.c10.Tm.CAPG", "pan- T: CD4.c24.Mix.NME2" ,
#                               "pan- T: CD4.c21.Treg.OAS1", "pan- T: CD4.c14.Th17.SLC4A10", "pan- T: CD4.c09.Tm.CCL5", "pan- T: CD4.c02.Tn.PASK" ,
#                               "pan- T: CD4.c03.Tn.ADSL", "pan- T: CD4.c05.Tm.TNF", "pan- T: CD4.c08.Tm.CREM", "pan- T: CD4.c04.Tn.il7r" ,
#                               "pan- T: CD4.c15.Th17.IL23R", "pan- T: CD4.c23.Mix.NME1",
#                               
#                               "pan- blueprint: C3_CD8_ZNF683", "pan- blueprint: C1_CD8_HAVCR2" ,
#                               "pan- blueprint: C2_CD8_GZMK", "pan- blueprint: C6_CD4_GZMA", "pan- blueprint: C10_NK_XCL1", "pan- blueprint: C8_CD4_FOXP3",
#                               "pan- blueprint: C5_CD4_CCR7", "pan- blueprint: C4_CD8_CX3CR1", "pan- blueprint: C7_CD4_CXCL13", "pan- blueprint: C9_NK_FGFBP2" ,
#                               
#                               
#                               "ABTC: Naive.CD4.T.cell", "ABTC: Exhausted.CD8.T.cell", "ABTC: Memory.CD8.T.cell", "ABTC: Treg" ,
#                               "ABTC: DNT", "ABTC: Effector.CD8.T.cell", "ABTC: NCAM1.FCGR3A..NK", "ABTC: Naive.CD8.T.cell" ,
#                               "ABTC: NKT.and.NK.cell", "ABTC: NCAM1.FCGR3A..NK.1", "ABTC: NCAM1.FCGR3A..NK.2", "HCC: 13_T.NK" ,
#                               
#                               
#                               "HCC: 25_T.NK", "HCC: 2_T.NK", "HCC: 20_T.NK", "HCC: 9_T.NK" ,
#                               "HCC: 8_T.NK", "HCC: 35_T.NK", "HCC: 18_T.NK", "HCC: 11_T.NK" ,
#                               "HCC: 14_T.NK", "HCC: 32_T.NK", "HCC: 1_T.NK", "HCC: 28_T.NK" ,
#                               "HCC: 51_T.NK",
#                               
#                               "STAD: LNK1", "STAD: LNK3", "STAD: LNK4" ,
#                               "STAD: LT2", "STAD: LT3", "STAD: LNK2", "STAD: LT4" ,
#                               "STAD: LT5", "STAD: LT1", "STAD: LPT")
names(new_cluster_names) <- c("TICA: Recently activated CD4 T cells",   "TICA: Naive-memory CD4 T cells",         "TICA: Transitional memory CD4 T cells",
                              "TICA: Cytotoxic CD8 T cells",            "TICA: Effector memory CD8 T cells",      "TICA: Th17 cells",
                              "TICA: NK",                               "TICA: Terminally exhausted CD8 T cells", "TICA: Naive T cells",
                              "TICA: Regulatory T cells",               "TICA: T helper cells",                   "TICA: Proliferative T cells",
                              "TICA: Pre-exhausted CD8 T cells",
                              
                              "pan- T: CD8.c07.Temra.CX3CR1",   "pan- T: CD8.c09.Tk.KIR2DL4",
                              "pan- T: CD8.c05.Tem.CXCR5",      "pan- T: CD8.c15.ISG.IFIT1",      "pan- T: CD8.c01.Tn.MAL",
                              "pan- T: CD8.c17.Tm.NME1",        "pan- T: CD8.c08.Tk.TYROBP",      "pan- T: CD8.c12.Tex.CXCL13",
                              "pan- T: CD8.c02.Tm.IL7R",        "pan- T: CD8.c10.Trm.ZNF683",     "pan- T: CD8.c03.Tm.RPS12",
                              "pan- T: CD8.c14.Tex.TCF7",       "pan- T: CD8.c06.Tem.GZMK",       "pan- T: CD8.c11.Tex.PDCD1",
                              "pan- T: CD8.c16.MAIT.SLC4A10",   "pan- T: CD8.c04.Tm.CD52",        "pan- T: CD8.c13.Tex.myl12a",
                              "pan- T: CD4.c20.Treg.TNFRSF9",   "pan- T: CD4.c18.Treg.RTKN2",     "pan- T: CD4.c06.Tm.ANXA1",
                              "pan- T: CD4.c22.ISG.IFIT1",      "pan- T: CD4.c01.Tn.TCF7",        "pan- T: CD4.c12.Tem.GZMK",
                              "pan- T: CD4.c11.Tm.GZMA",        "pan- T: CD4.c13.Temra.CX3CR1",   "pan- T: CD4.c16.Tfh.CXCR5",
                              "pan- T: CD4.c07.Tm.ANXA2",       "pan- T: CD4.c17.TfhTh1.CXCL13",  "pan- T: CD4.c19.Treg.S1PR1",
                              "pan- T: CD4.c10.Tm.CAPG",        "pan- T: CD4.c24.Mix.NME2",       "pan- T: CD4.c21.Treg.OAS1",
                              "pan- T: CD4.c14.Th17.SLC4A10",   "pan- T: CD4.c09.Tm.CCL5",        "pan- T: CD4.c02.Tn.PASK",
                              "pan- T: CD4.c03.Tn.ADSL",        "pan- T: CD4.c05.Tm.TNF",         "pan- T: CD4.c08.Tm.CREM",
                              "pan- T: CD4.c04.Tn.il7r",        "pan- T: CD4.c15.Th17.IL23R",     "pan- T: CD4.c23.Mix.NME1",
                              
                              
                              "pan- blueprint: C3_CD8_ZNF683",     "pan- blueprint: C1_CD8_HAVCR2",     "pan- blueprint: C2_CD8_GZMK",
                              "pan- blueprint: C6_CD4_GZMA",       "pan- blueprint: C10_NK_XCL1",       "pan- blueprint: C8_CD4_FOXP3",
                              "pan- blueprint: C5_CD4_CCR7",       "pan- blueprint: C4_CD8_CX3CR1",     "pan- blueprint: C7_CD4_CXCL13",
                              "pan- blueprint: C9_NK_FGFBP2",
                              
                              "ABTC: Naive CD4 T cell",                 "ABTC: Exhausted CD8 T cell",
                              "ABTC: Memory CD8 T cell",                "ABTC: Treg",                             "ABTC: DNT",
                              "ABTC: Effector CD8 T cell",              "ABTC: NCAM1-FCGR3A+ NK",                 "ABTC: Naive CD8 T cell",
                              "ABTC: NKT and NK cell",                  "ABTC: NCAM1+FCGR3A+ NK",                 "ABTC: NCAM1+FCGR3A- NK",
                              
                              
                              "HCC: 13_T/NK",                           "HCC: 25_T/NK",                           "HCC: 2_T/NK",
                              "HCC: 20_T/NK",                           "HCC: 9_T/NK",                            "HCC: 8_T/NK",
                              "HCC: 35_T/NK",                           "HCC: 18_T/NK",                           "HCC: 11_T/NK",
                              "HCC: 14_T/NK",                           "HCC: 32_T/NK",                           "HCC: 1_T/NK",
                              "HCC: 28_T/NK",                           "HCC: 51_T/NK",                           "STAD: LNK1",
                              
                              
                              "STAD: LNK3",                             "STAD: LNK4",                             "STAD: LT2",
                              "STAD: LT3",                              "STAD: LNK2",                             "STAD: LT4",
                              "STAD: LT5",                              "STAD: LT1",                              "STAD: LPT")

pbmc3k_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_TICA_T_NK_cells.rds"))
pbmc3k_T_NK_cells_TICA <- subset(pbmc3k_T_NK_cells_TICA,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
cell_summary_TICA <- data.frame(pbmc3k_T_NK_cells_TICA$cell_types_from_pbmc3k,paste0("TICA: ",pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells))
# rm("pbmc3k_T_NK_cells_TICA")
pbmc3k_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
pbmc3k_T_NK_cells_panT <- subset(pbmc3k_T_NK_cells_panT,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
cell_summary_panT <- data.frame(pbmc3k_T_NK_cells_panT$cell_types_from_pbmc3k,paste0("pan- T: ",pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4))
# rm("pbmc3k_T_NK_cells_panT")
pbmc3k_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
pbmc3k_T_NK_cells_panblueprint <- subset(pbmc3k_T_NK_cells_panblueprint,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
cell_summary_panblueprint <- data.frame(pbmc3k_T_NK_cells_panblueprint$cell_types_from_pbmc3k,paste0("pan- blueprint: ",pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))
# rm("pbmc3k_T_NK_cells_panblueprint")
pbmc3k_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_ABTC_T_NK_cells.rds"))
pbmc3k_T_NK_cells_ABTC <- subset(pbmc3k_T_NK_cells_ABTC,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
cell_summary_ABTC <- data.frame(pbmc3k_T_NK_cells_ABTC$cell_types_from_pbmc3k,paste0("ABTC: ",pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells))
# rm("pbmc3k_T_NK_cells_ABTC")
pbmc3k_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_HCC_T_NK_cells.rds"))
pbmc3k_T_NK_cells_HCC <- subset(pbmc3k_T_NK_cells_HCC,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
cell_summary_HCC <- data.frame(pbmc3k_T_NK_cells_HCC$cell_types_from_pbmc3k,paste0("HCC: ",pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells))
# rm("pbmc3k_T_NK_cells_HCC")
pbmc3k_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_STAD_T_NK_cells.rds"))
pbmc3k_T_NK_cells_STAD <- subset(pbmc3k_T_NK_cells_STAD,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
cell_summary_STAD <- data.frame(pbmc3k_T_NK_cells_STAD$cell_types_from_pbmc3k,paste0("STAD: ",pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells),pbmc3k_T_NK_cells_STAD$hpca.fine)
# rm("pbmc3k_T_NK_cells_STAD")

cell_summary <- cbind(cell_summary_TICA,to_panT=cell_summary_panT[,2],
                      to_panblueprint=cell_summary_panblueprint[,2],to_ABTC=cell_summary_ABTC[,2],to_HCC=cell_summary_HCC[,2],
                      to_STAD=cell_summary_STAD[,2],hpca_fine=cell_summary_STAD[,3])
colnames(cell_summary)[1:2] <-c("cell_types_from_pbmc3k","to_TICA")
cell_summary[1:3,]
cell_summary[,"cell_types_from_pbmc3k"] <- as.character(cell_summary[,"cell_types_from_pbmc3k"])
cell_summary[,"hpca_fine"] <- as.character(cell_summary[,"hpca_fine"])

unique(cell_summary$cell_types_from_pbmc3k)
# [1] "0_T memory cells"      "2_T memory cells"      "6_Gamma delta T cells" "4_T cells"     

index <- c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD")
index2 <- c("TICA: ","pan- T: ","pan- blueprint: ","ABTC: ","HCC: ","STAD: ")
# for T_cell:CD4+_effector_memory
output1 <- NULL;output2 <- NULL;output3<-NULL
pp <- NULL
for (i in 1:length(index)) {
  #### for p1
  
  cell_summary_subset1 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("0_T memory cells","2_T memory cells")) #
  cell_summary_subset1$hpca_fine <- "CD4+ T"
  # cell_summary_subset2 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
  # cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
  # cell_summary_subset3 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
  # cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
  cell_summary_subset <- cell_summary_subset1 #rbind(cell_summary_subset1,cell_summary_subset2,cell_summary_subset3)
  output_pre1 <- as.data.frame(table(cell_summary_subset1[,c(index[i],"hpca_fine")]))
  output_pre1 <- data.frame(group=index[i],output_pre1)
  output_pre1[,2] <- as.character(output_pre1[,2])
  output_pre1$Freq <- as.numeric(output_pre1$Freq)
  # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
  colnames(output_pre1)[2] <- "fill"
  output_pre1$fill <- as.character(new_cluster_names[output_pre1$fill])
  output_pre1 <- data.frame(output_pre1,populations=output_pre1$fill)
  output_pre1$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre1$fill),": ",2)))[,2])
  output_pre1$Freq <- as.numeric(output_pre1$Freq)
  output_pre1 <- output_pre1[order(output_pre1$Freq),]
  output_pre1 <- output_pre1 %>%
    group_by(group) %>%
    mutate(label_y = cumsum(Freq)-0.5*Freq)
  output_pre1 <- output_pre1[order(-output_pre1$Freq),]
  output1 <- rbind(output1,output_pre1)

output1 <- as.data.frame(output1)
output1 <- rbind(output1,c("RF","CD4+ T","CD4+ T",sum(output_pre1$Freq),"RF: CD4+ T",sum(output_pre1$Freq)/2))#;output1 <- output
# output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
# output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2));output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
# output <- rbind(output1,output2,output3)
output1$fill <- factor(output1$fill,levels = c(output1$fill))
output1$populations <- factor(output1$populations,levels = c(output1$populations))
output1$group <- factor(output1$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
output1$Freq <- as.numeric(output1$Freq)
output1$label_y <- as.numeric(output1$label_y)
output1[output1$Freq < 3,"fill"] <- "" #10
if (i ==1) {
  cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
  names(cols) <- paste0(index2[i],sort(unique(pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
}
if (i ==2) {
cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
names(cols) <- paste0(index2[i],sort(unique(pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
}
  if (i ==3) {
cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
names(cols) <- sort(as.character(new_cluster_names[paste0("pan- blueprint: ",unique(pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))]))
}
    if (i ==4) {
cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
}
      if (i ==5) {
cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
}
        if (i ==6) {
cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
names(cols) <- sort(as.character(new_cluster_names[paste0("STAD: ",unique(pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells))]))
}



cols <- c("RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
p1 <- ggplot(output1,aes(x=group,y=Freq,fill=populations)) + 
  geom_bar(stat="identity",colour="black") + 
  scale_fill_manual(values=cols) +
  # geom_text(aes(y=label_y,label=fill),size =10) +
  theme_classic()


##### for p2

# cell_summary_subset2 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("0_T memory cells","2_T memory cells")) #
# cell_summary_subset2$hpca_fine <- "CD4+ T"
cell_summary_subset2 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
# cell_summary_subset3 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
# cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
# cell_summary_subset <- cell_summary_subset2 #rbind(cell_summary_subset2,cell_summary_subset2,cell_summary_subset3)
output_pre2 <- as.data.frame(table(cell_summary_subset2[,c(index[i],"hpca_fine")]))
output_pre2 <- data.frame(group=index[i],output_pre2)
output_pre2[,2] <- as.character(output_pre2[,2])
output_pre2$Freq <- as.numeric(output_pre2$Freq)
# output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
colnames(output_pre2)[2] <- "fill"
output_pre2$fill <- as.character(new_cluster_names[output_pre2$fill])
output_pre2 <- data.frame(output_pre2,populations=output_pre2$fill)
output_pre2$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre2$fill),": ",2)))[,2])
output_pre2$Freq <- as.numeric(output_pre2$Freq)
output_pre2 <- output_pre2[order(output_pre2$Freq),]
output_pre2 <- output_pre2 %>%
  group_by(group) %>%
  mutate(label_y = cumsum(Freq)-0.5*Freq)
output_pre2 <- output_pre2[order(-output_pre2$Freq),]
output2 <- rbind(output2,output_pre2)

output2 <- as.data.frame(output2)
# output2 <- rbind(output2,c("RF","CD4+ T","CD4+ T",sum(output_pre2$Freq),"RF: CD4+ T",sum(output_pre2$Freq)/2))#;output2 <- output
output2 <- rbind(output2,c("RF","CD8+ T","CD8+ T",sum(output_pre2$Freq),"RF: CD8+ T",sum(output_pre2$Freq)/2))# ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
# output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2));output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
# output <- rbind(output2,output2,output3)
output2$fill <- factor(output2$fill,levels = c(output2$fill))
output2$populations <- factor(output2$populations,levels = c(output2$populations))
output2$group <- factor(output2$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
output2$Freq <- as.numeric(output2$Freq)
output2$label_y <- as.numeric(output2$label_y)
output2[output2$Freq < 3,"fill"] <- "" #20
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
# names(cols) <- paste0(index2[i],sort(unique(pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
# names(cols) <- paste0(index2[i],sort(unique(pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("pan- blueprint: ",sort(unique(pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))])
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("STAD: ",sort(unique(pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))])



cols <- c("RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
p2 <- ggplot(output2,aes(x=group,y=Freq,fill=populations)) + 
  geom_bar(stat="identity",colour="black") + 
  scale_fill_manual(values=cols) +
  # geom_text(aes(y=label_y,label=fill),size =20) +
  theme_classic()


#### for p3

# cell_summary_subset3 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("0_T memory cells","2_T memory cells")) #
# cell_summary_subset3$hpca_fine <- "CD4+ T"
# cell_summary_subset2 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
# cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
cell_summary_subset3 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
# cell_summary_subset <- cell_summary_subset3 #rbind(cell_summary_subset3,cell_summary_subset2,cell_summary_subset3)
output_pre3 <- as.data.frame(table(cell_summary_subset3[,c(index[i],"hpca_fine")]))
output_pre3 <- data.frame(group=index[i],output_pre3)
output_pre3[,2] <- as.character(output_pre3[,2])
output_pre3$Freq <- as.numeric(output_pre3$Freq)
# output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
colnames(output_pre3)[2] <- "fill"
output_pre3$fill <- as.character(new_cluster_names[output_pre3$fill])
output_pre3 <- data.frame(output_pre3,populations=output_pre3$fill)
output_pre3$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre3$fill),": ",2)))[,2])
output_pre3$Freq <- as.numeric(output_pre3$Freq)
output_pre3 <- output_pre3[order(output_pre3$Freq),]
output_pre3 <- output_pre3 %>%
  group_by(group) %>%
  mutate(label_y = cumsum(Freq)-0.5*Freq)
output_pre3 <- output_pre3[order(-output_pre3$Freq),]
output3 <- rbind(output3,output_pre3)

output3 <- as.data.frame(output3)
# output3 <- rbind(output3,c("RF","CD4+ T","CD4+ T",sum(output_pre3$Freq),"RF: CD4+ T",sum(output_pre3$Freq)/2))#;output3 <- output
# output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
output3 <- rbind(output3,c("RF","NK cell","NK",sum(output_pre3$Freq),"RF: NK",sum(output_pre3$Freq)/2))#;output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
# output <- rbind(output3,output2,output3)
output3$fill <- factor(output3$fill,levels = c(output3$fill))
output3$populations <- factor(output3$populations,levels = c(output3$populations))
output3$group <- factor(output3$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
output3$Freq <- as.numeric(output3$Freq)
output3$label_y <- as.numeric(output3$label_y)
output3[output3$Freq < 3,"fill"] <- "" #30
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
# names(cols) <- paste0(index2[i],sort(unique(pbmc3k_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
# names(cols) <- paste0(index2[i],sort(unique(pbmc3k_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("pan- blueprint: ",sort(unique(pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))])
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(pbmc3k_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
# cols <- hue_pal()(length(unique(pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
# names(cols) <- as.character(new_cluster_names[paste0("STAD: ",sort(unique(pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))])


cols <- c("RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
p3 <- ggplot(output3,aes(x=group,y=Freq,fill=populations)) + 
  geom_bar(stat="identity",colour="black") + 
  scale_fill_manual(values=cols) +
  # geom_text(aes(y=label_y,label=fill),size =30) +
  theme_classic()


# CombinePlots(plots=list(p1,p2,p3),nrow=1,legend="right")

p1 <- p1 + theme(axis.text.x = element_text(size = 28,color = "black"),axis.text.y = element_text(size = 28,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="")
p2 <- p2 + theme(axis.text.x = element_text(size = 28,color = "black"),axis.text.y = element_text(size = 28,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="")
p3 <- p3 + theme(axis.text.x = element_text(size = 28,color = "black"),axis.text.y = element_text(size = 28,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="")
CombinePlots(plots=list(p1,p2,p3),nrow=1,legend="none")
# pdf(paste0(workdir,"test.pdf"))
# pp <- list(pp,CombinePlots(plots=list(p1,p2,p3),nrow=1,legend="right"))
# dev.off()
}


########  check the agreement

cell_summary_subset1 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("0_T memory cells","2_T memory cells")) #
cell_summary_subset1$hpca_fine <- "CD4+ T"
cell_summary_subset2 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
cell_summary_subset3 <- subset(cell_summary,cell_types_from_pbmc3k %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
cell_summary_subset123 <- rbind(cell_summary_subset1,cell_summary_subset2,cell_summary_subset3)
cell_summary_subset123$to_HCC <- as.character(new_cluster_names[cell_summary_subset123$to_HCC])
cell_summary_subset123$to_STAD <- as.character(new_cluster_names[cell_summary_subset123$to_STAD])
table(cell_summary_subset123[,c("to_TICA","hpca_fine")])

# hpca_fine
# to_TICA                                CD4+ T CD8+ T  NK
# TICA: Cytotoxic CD8 T cells               9     60   3
# TICA: Naive T cells                     450     24   0
# TICA: Naive-memory CD4 T cells          591     79   2
# TICA: NK                                  0     49 130
# TICA: Pre-exhausted CD8 T cells          24      0   2
# TICA: Proliferative T cells               1      1  11
# TICA: Recently activated CD4 T cells      6      2   0
# TICA: Regulatory T cells                  2      0   0
# TICA: T helper cells                      4      0   0  ## 60+450+591+130+6+2+4 / 1450 = 0.8572414

# hpca_fine
# to_panT                        CD4+ T CD8+ T   NK
# pan- T: CD4.c03.Tn.ADSL        1047     71    1
# pan- T: CD4.c10.Tm.CAPG           5     47    0
# pan- T: CD4.c24.Mix.NME2          2      0    0
# pan- T: CD8.c07.Temra.CX3CR1      0     80  139
# pan- T: CD8.c10.Trm.ZNF683        3      0    8
# pan- T: CD8.c17.Tm.NME1          30     17    0## 1047+5+2+80+17 / 1450 = 0.7937931

# hpca_fine
# to_panblueprint                 CD4+ T CD8+ T  NK
# pan- blueprint: C1_CD8_HAVCR2      0      0   6
# pan- blueprint: C10_NK_XCL1        0     13   4
# pan- blueprint: C2_CD8_GZMK        6     83   0
# pan- blueprint: C3_CD8_ZNF683      7     11   0
# pan- blueprint: C4_CD8_CX3CR1      0     87   4
# pan- blueprint: C5_CD4_CCR7      942     11   0
# pan- blueprint: C6_CD4_GZMA      123      6   2
# pan- blueprint: C8_CD4_FOXP3       9      0   0
# pan- blueprint: C9_NK_FGFBP2       0      4 132 ## 83+11+87+942+123+9+132 / 1450 = 0.9565517

# hpca_fine
# to_ABTC                      CD4+ T CD8+ T  NK
# ABTC: DNT                      46      8   0
# ABTC: Exhausted CD8 T cell      0      0   6
# ABTC: Memory CD8 T cell         9    168   3
# ABTC: Naive CD4 T cell        991     15   0
# ABTC: Naive CD8 T cell         35      4   0
# ABTC: NCAM1-FCGR3A+ NK          1      9 138
# ABTC: NCAM1+FCGR3A- NK          0     11   1
# ABTC: Treg                      5      0   0 ## 168+991+4+138+1+5 / 1450 = 0.9013793


# hpca_fine
# to_HCC                          CD4+ T CD8+ T  NK
# HCC: 1_CD8_cytotoxic T             0      0   1
# HCC: 11_CD8_inter-states T        41     60   0
# HCC: 13_CD8_effector memory T     11     83   3
# HCC: 14_CD8_TR_memory T            0      2   0
# HCC: 2_CD8_MAIT                  345     26   0
# HCC: 20_CD4_naive T              680     11   0
# HCC: 28_NK                         1     32 135
# HCC: 32_CD4_other T                0      0   8
# HCC: 8_NK                          0      0   1
# HCC: 9_CD4_regulatory T            9      1   0 ## 60+83+2+26+680+135+1+9 / 1450 = 0.6868966

# hpca_fine
# to_STAD                     CD4+ T CD8+ T  NK
# STAD: CD4 helper T            16      0   0
# STAD: CD4 naive T cell       975     20   0
# STAD: CD8 effector T cell      0      1   0
# STAD: CD8 naive T             31     25   0
# STAD: LNK1                     0     96 139
# STAD: LNK2                    21     69   0
# STAD: Proliferative T          1      2   9
# STAD: Regulatory T            43      2   0 ## 16+975+1+25+139+43 / 1450 = 0.8268966



# output$fill <- as.character(new_cluster_names[output$fill])
# output$fill <- as.character(t(as.data.frame(str_split(as.character(output$fill),"__",2)))[,2])
# output$Freq <- as.numeric(output$Freq)
# output <- output %>%
#   group_by(group) %>%
#   mutate(Freq=sort(Freq))
# 
# 
# output <- output %>%
#   group_by(group) %>%
#   mutate(label_y = cumsum(Freq)-0.5*Freq)
# output %>%
#   group_by(group) %>%
#   mutate(Freq=-sort(Freq))


##########################################
#### 2023.06.22 plot umpas
###########################################
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
changenames <- c("CD4 T","CD4 T","CD8 T","NK")
names(changenames) <- c("0_T memory cells","2_T memory cells","4_T cells","6_Gamma delta T cells")
new_cluster_names <- c("Recently activated CD4 T cells", "Naive-memory CD4 T cells", "Transitional memory CD4 T cells","Cytotoxic CD8 T cells",
                       "Effector memory CD8 T cells", "Th17 cells", "NK", "Terminally exhausted CD8 T cells",
                       "Naive T cells", "Regulatory T cells", "T helper cells", "Proliferative T cells",
                       "Pre-exhausted CD8 T cells",
                       
                       
                       "CD8.c07.Temra.CX3CR1", "CD8.c09.Tk.KIR2DL4", "CD8.c05.Tem.CXCR5" ,
                       "CD8.c15.ISG.IFIT1", "CD8.c01.Tn.MAL", "CD8.c17.Tm.NME1", "CD8.c08.Tk.TYROBP" ,
                       "CD8.c12.Tex.CXCL13", "CD8.c02.Tm.IL7R", "CD8.c10.Trm.ZNF683", "CD8.c03.Tm.RPS12" ,
                       "CD8.c14.Tex.TCF7", "CD8.c06.Tem.GZMK", "CD8.c11.Tex.PDCD1", "CD8.c16.MAIT.SLC4A10" ,
                       "CD8.c04.Tm.CD52", "CD8.c13.Tex.myl12a", "CD4.c20.Treg.TNFRSF9", "CD4.c18.Treg.RTKN2" ,
                       "CD4.c06.Tm.ANXA1", "CD4.c22.ISG.IFIT1", "CD4.c01.Tn.TCF7", "CD4.c12.Tem.GZMK" ,
                       "CD4.c11.Tm.GZMA", "CD4.c13.Temra.CX3CR1", "CD4.c16.Tfh.CXCR5", "CD4.c07.Tm.ANXA2" ,
                       "CD4.c17.TfhTh1.CXCL13", "CD4.c19.Treg.S1PR1", "CD4.c10.Tm.CAPG", "CD4.c24.Mix.NME2" ,
                       "CD4.c21.Treg.OAS1", "CD4.c14.Th17.SLC4A10", "CD4.c09.Tm.CCL5", "CD4.c02.Tn.PASK" ,
                       "CD4.c03.Tn.ADSL", "CD4.c05.Tm.TNF", "CD4.c08.Tm.CREM", "CD4.c04.Tn.il7r" ,
                       "CD4.c15.Th17.IL23R", "CD4.c23.Mix.NME1",
                       
                       
                       "CD8 memory T", "CD8 exhausted cytotoxic T",
                       "CD8 pre-effector T",   "CD4 memory/effector T",
                       "NK_XCL1",   "CD4 regulatory T",
                       "CD4 naive T",   "CD8 effector T",
                       "CD4 exhausted effector T", "Cytotoxic NK",
                       
                       "Naive CD4 T cell",     "Exhausted CD8 T cell",
                       "Memory CD8 T cell",    "Treg" ,               
                       "DNT",                  "Effector CD8 T cell" ,
                       "NCAM1nFCGR3Ap NK",     "Naive CD8 T cell" ,   
                       "NKT and NK cell",      "NCAM1pFCGR3Ap NK"   ,"NCAM1pFCGR3An NK",
                       
                       "13_CD8_effector memory T", "2_NK", "2_CD8_MAIT",  "20_CD4_naive T",
                       "9_CD4_regulatory T",  "8_NK",  "35_CD4_central memory T", "18_NK",
                       "11_CD8_inter-states T", "14_CD8_TR_memory T", "32_CD4_other T", "1_CD8_cytotoxic T" ,
                       "28_NK", "51_CD8_other T",
                       
                       "LNK1", "LNK3", "LNK4", "CD8 naive T" ,
                       "CD4 naive T cell",  "LNK2", "CD4 helper T",  "Regulatory T" ,
                       "CD8 effector T cell",  "Proliferative T")

names(new_cluster_names) <- c("Recently activated CD4 T cells",   "Naive-memory CD4 T cells",         "Transitional memory CD4 T cells",
                              "Cytotoxic CD8 T cells",            "Effector memory CD8 T cells",      "Th17 cells",
                              "NK",                               "Terminally exhausted CD8 T cells", "Naive T cells",
                              "Regulatory T cells",               "T helper cells",                   "Proliferative T cells",
                              "Pre-exhausted CD8 T cells",
                              
                              "CD8.c07.Temra.CX3CR1",   "CD8.c09.Tk.KIR2DL4",
                              "CD8.c05.Tem.CXCR5",      "CD8.c15.ISG.IFIT1",      "CD8.c01.Tn.MAL",
                              "CD8.c17.Tm.NME1",        "CD8.c08.Tk.TYROBP",      "CD8.c12.Tex.CXCL13",
                              "CD8.c02.Tm.IL7R",        "CD8.c10.Trm.ZNF683",     "CD8.c03.Tm.RPS12",
                              "CD8.c14.Tex.TCF7",       "CD8.c06.Tem.GZMK",       "CD8.c11.Tex.PDCD1",
                              "CD8.c16.MAIT.SLC4A10",   "CD8.c04.Tm.CD52",        "CD8.c13.Tex.myl12a",
                              "CD4.c20.Treg.TNFRSF9",   "CD4.c18.Treg.RTKN2",     "CD4.c06.Tm.ANXA1",
                              "CD4.c22.ISG.IFIT1",      "CD4.c01.Tn.TCF7",        "CD4.c12.Tem.GZMK",
                              "CD4.c11.Tm.GZMA",        "CD4.c13.Temra.CX3CR1",   "CD4.c16.Tfh.CXCR5",
                              "CD4.c07.Tm.ANXA2",       "CD4.c17.TfhTh1.CXCL13",  "CD4.c19.Treg.S1PR1",
                              "CD4.c10.Tm.CAPG",        "CD4.c24.Mix.NME2",       "CD4.c21.Treg.OAS1",
                              "CD4.c14.Th17.SLC4A10",   "CD4.c09.Tm.CCL5",        "CD4.c02.Tn.PASK",
                              "CD4.c03.Tn.ADSL",        "CD4.c05.Tm.TNF",         "CD4.c08.Tm.CREM",
                              "CD4.c04.Tn.il7r",        "CD4.c15.Th17.IL23R",     "CD4.c23.Mix.NME1",
                              
                              
                              "C3_CD8_ZNF683",     "C1_CD8_HAVCR2",     "C2_CD8_GZMK",
                              "C6_CD4_GZMA",       "C10_NK_XCL1",       "C8_CD4_FOXP3",
                              "C5_CD4_CCR7",       "C4_CD8_CX3CR1",     "C7_CD4_CXCL13",
                              "C9_NK_FGFBP2",
                              
                              "Naive CD4 T cell",                 "Exhausted CD8 T cell",
                              "Memory CD8 T cell",                "Treg",                             "DNT",
                              "Effector CD8 T cell",              "NCAM1-FCGR3A+ NK",                 "Naive CD8 T cell",
                              "NKT and NK cell",                  "NCAM1+FCGR3A+ NK",                 "NCAM1+FCGR3A- NK",
                              
                              
                              "13_T/NK",                           "25_T/NK",                           "2_T/NK",
                              "20_T/NK",                           "9_T/NK",                            "8_T/NK",
                              "35_T/NK",                           "18_T/NK",                           "11_T/NK",
                              "14_T/NK",                           "32_T/NK",                           "1_T/NK",
                              "28_T/NK",                           "51_T/NK",                           "LNK1",
                              
                              
                              "LNK3",                             "LNK4",                             "LT2",
                              "LT3",                              "LNK2",                             "LT4",
                              "LT5",                              "LT1",                              "LPT")

pbmc3k_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_TICA_T_NK_cells.rds"))
pbmc3k_T_NK_cells_TICA <- subset(pbmc3k_T_NK_cells_TICA,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
pbmc3k_T_NK_cells_TICA$cell_types_from_pbmc3k <- as.character(changenames[as.character(pbmc3k_T_NK_cells_TICA$cell_types_from_pbmc3k)])
DimPlot(pbmc3k_T_NK_cells_TICA,group.by = "cell_types_from_pbmc3k",label = T,pt.size=2,label.size=12,cols = c("green","blue","red"))
DimPlot(pbmc3k_T_NK_cells_TICA,group.by = "predicted.to_TICA_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=24)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

pbmc3k_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
pbmc3k_T_NK_cells_panT <- subset(pbmc3k_T_NK_cells_panT,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
DimPlot(pbmc3k_T_NK_cells_panT,group.by = "predicted.to_pancancer_T_cells_CD8_CD4",pt.size=2) + theme(legend.text=element_text(size=24)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

pbmc3k_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
pbmc3k_T_NK_cells_panblueprint <- subset(pbmc3k_T_NK_cells_panblueprint,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells <- as.character(new_cluster_names[as.character(pbmc3k_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)])
DimPlot(pbmc3k_T_NK_cells_panblueprint,group.by = "predicted.to_pancancer_blueprint_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=24)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

pbmc3k_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_ABTC_T_NK_cells.rds"))
pbmc3k_T_NK_cells_panblueprint <- subset(pbmc3k_T_NK_cells_panblueprint,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
DimPlot(pbmc3k_T_NK_cells_ABTC,group.by = "predicted.to_ABTC_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=24)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

pbmc3k_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_HCC_T_NK_cells.rds"))
pbmc3k_T_NK_cells_HCC <- subset(pbmc3k_T_NK_cells_HCC,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells <- as.character(new_cluster_names[as.character(pbmc3k_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)])
DimPlot(pbmc3k_T_NK_cells_HCC,group.by = "predicted.to_HCC_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=24)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

pbmc3k_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_pbmc3k_T_NK_cells_to_STAD_T_NK_cells.rds"))
pbmc3k_T_NK_cells_STAD <- subset(pbmc3k_T_NK_cells_STAD,ident = c("0_T memory cells","1_T memory cells","3_NK cells","4_Gamma delta T cells"))
pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells <- as.character(new_cluster_names[as.character(pbmc3k_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)])
DimPlot(pbmc3k_T_NK_cells_STAD,group.by = "predicted.to_STAD_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=24)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )



########################################
#### 2024.02.21 Jing Yang
#### using bi dataset as known sc dataset
########################################
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
"%ni%" <- Negate("%in%")

# load("/Users/jingyang/Dropbox/work/cancer-immu for biological discoveries/code/scRNA-seq_from_Qi_object.Rdata")
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
Bi_T_NK <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_20250506.rds"))
DimPlot(Bi_T_NK,group.by="cell_types_from_author",label = T)
DimPlot(Bi_T_NK,label = T)
DimPlot(Bi_T_NK,group.by = "CellType",label = T)

Bi_T_NK <- subset(Bi_T_NK, idents = c("3_T cells"), invert = TRUE)
DimPlot(Bi_T_NK,label = T)
FeaturePlot(Bi_T_NK,features = c("CD4","CD8A","CD8B","FOXP3","TNFRSF4","KLRC1","STMN1","MKI67","PDCD1","HAVCR2","IFIT3","GZMB","PRF1","GZMA","GZMH","ZNF683","NCR1","NCAM1","FCGR3A"))

# new.cluster.ids <- c("0_CD8+ T cells","1_CD8+ T cells","2_T reg",,"4_Proliferative T cells","5_NK cells","6_T helpers","7_CD8+ T cells","8_CD8+ T cells")
new.cluster.ids <- c("0_CD8+ T cells","1_CD8+ T cells","2_T reg","3_CD4+CD8+ T cells","4_CD4+ T cells","5_NK cells","6_Proliferative T cells","7_CD8+ T cells")
names(new.cluster.ids) <- levels(Bi_T_NK)
Bi_T_NK <- RenameIdents(Bi_T_NK, new.cluster.ids)
Bi_T_NK$cellcluster_ann <- Idents(Bi_T_NK)
DimPlot(Bi_T_NK,label = T)


### mapping to atlases using Bi dataset 
Bi_T_NK_cellsbubian <- Bi_T_NK

#### 1. TICA

atlas_TICA_T_NK_cells <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells.rds"))
atlas_TICA_T_NK_cells <- subset(atlas_TICA_T_NK_cells,downsample=1000)
Bi_T_NK_cells <- Bi_T_NK_cellsbubian
anchors_map_Bi_T_NK_cells_to_TICA_T_NK_cells <- FindTransferAnchors(
  reference = atlas_TICA_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = Bi_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
Bi_T_NK_cells <- MapQuery(
  anchorset = anchors_map_Bi_T_NK_cells_to_TICA_T_NK_cells,
  query = Bi_T_NK_cells,
  reference = atlas_TICA_T_NK_cells,
  refdata = list(
    to_TICA_T_NK_cells = "cell_types_from_author",
    to_TICA_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_TICA_T_NK_cells = "integrated"#"RNA"
  ),
  new.reduction.name = "to_TICA_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_Bi_T_NK_cells_to_TICA_T_NK_cells")
rm("atlas_TICA_T_NK_cells")
Bi_T_NK_cells[["to_TICA_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(Bi_T_NK_cells,reduction = "ref.umap"),
                                                                   key = "to_TICA_T_NK_cells_",
                                                                   assay = DefaultAssay(Bi_T_NK_cells))
saveRDS(Bi_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_integrated_20250506.rds"))
rm("Bi_T_NK_cells")

#### pancancer T cells CD8 CD4
Bi_T_NK_cells <- Bi_T_NK_cellsbubian
atlas_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_filled.rds"))
anchors_map_Bi_T_NK_cells_to_pancancer_T_cells_CD8_CD4 <- FindTransferAnchors(
  reference = atlas_pancancer_T_cells_CD8_CD4,#atlas_pancancer_T_cells_CD8,
  query = Bi_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
Bi_T_NK_cells <- MapQuery(
  anchorset = anchors_map_Bi_T_NK_cells_to_pancancer_T_cells_CD8_CD4,
  query = Bi_T_NK_cells,
  reference = atlas_pancancer_T_cells_CD8_CD4,
  refdata = list(
    to_pancancer_T_cells_CD8_CD4 = "cell_types_from_author", #"cell_types_from_author_simplized"
    to_pancancer_T_cells_CD8_CD4_our_calculation = "cellcluster_ann",
    predicted_to_pancancer_T_cells_CD8_CD4 = "RNA"
  ),
  new.reduction.name = "to_pancancer_T_cells_CD8_CD4_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_Bi_T_NK_cells_to_pancancer_T_cells_CD8_CD4")
rm("atlas_pancancer_T_cells_CD8_CD4")
Bi_T_NK_cells[["to_pancancer_T_cells_CD8_CD4_umap"]] <- CreateDimReducObject(embeddings = Embeddings(Bi_T_NK_cells,reduction = "ref.umap"),
                                                                             key = "to_pancancer_T_cells_CD8_CD4_",
                                                                             assay = DefaultAssay(Bi_T_NK_cells))
saveRDS(Bi_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_20250506.rds"))
rm("Bi_T_NK_cells")

# 3. map to pancancer_blueprint_T_NK_cells
atlas_pancancer_blueprint_T_NK_cells <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells.rds"))
DefaultAssay(atlas_pancancer_blueprint_T_NK_cells) <- "integrated"
Bi_T_NK_cells <- Bi_T_NK_cellsbubian
anchors_map_Bi_T_NK_cells_to_pancancer_blueprint_T_NK_cells <- FindTransferAnchors(
  reference = atlas_pancancer_blueprint_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = Bi_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
Bi_T_NK_cells <- MapQuery(
  anchorset = anchors_map_Bi_T_NK_cells_to_pancancer_blueprint_T_NK_cells,
  query = Bi_T_NK_cells,
  reference = atlas_pancancer_blueprint_T_NK_cells,
  refdata = list(
    to_pancancer_blueprint_T_NK_cells = "cell_types_from_author",
    to_pancancer_blueprint_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_pancancer_blueprint_T_NK_cells = "integrated"#"RNA"
  ),
  new.reduction.name = "to_pancancer_blueprint_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_Bi_T_NK_cells_to_pancancer_blueprint_T_NK_cells")
rm("atlas_pancancer_blueprint_T_NK_cells")
Bi_T_NK_cells[["to_pancancer_blueprint_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(Bi_T_NK_cells,reduction = "ref.umap"),
                                                                                  key = "to_pancancer_blueprint_T_NK_cells_",
                                                                                  assay = DefaultAssay(Bi_T_NK_cells))
saveRDS(Bi_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated_20250506.rds"))
rm("Bi_T_NK_cells")

# 4. map to ABTC_T_NK_cells
Bi_T_NK_cells <- Bi_T_NK_cellsbubian
atlas_ABTC_T_NK_cells <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells.rds"))
anchors_map_Bi_T_NK_cells_to_ABTC_T_NK_cells <- FindTransferAnchors(
  reference = atlas_ABTC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = Bi_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
Bi_T_NK_cells <- MapQuery(
  anchorset = anchors_map_Bi_T_NK_cells_to_ABTC_T_NK_cells,
  query = Bi_T_NK_cells,
  reference = atlas_ABTC_T_NK_cells,
  refdata = list(
    to_ABTC_T_NK_cells = "cell_types_from_author",
    to_ABTC_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_ABTC_T_NK_cells = "RNA"
  ),
  new.reduction.name = "to_ABTC_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_Bi_T_NK_cells_to_ABTC_T_NK_cells")
rm("atlas_ABTC_T_NK_cells")
Bi_T_NK_cells[["to_ABTC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(Bi_T_NK_cells,reduction = "ref.umap"),
                                                                   key = "to_ABTC_T_NK_cells_",
                                                                   assay = DefaultAssay(Bi_T_NK_cells))
saveRDS(Bi_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_ABTC_T_NK_cells_20250506.rds"))
rm("Bi_T_NK_cells")

# 5. map to HCC_T_NK_cells
Bi_T_NK_cells <- Bi_T_NK_cellsbubian
atlas_HCC_T_NK_cells <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells.rds"))
DefaultAssay(atlas_HCC_T_NK_cells) <- "integrated"
anchors_map_Bi_T_NK_cells_to_HCC_T_NK_cells <- FindTransferAnchors(
  reference = atlas_HCC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = Bi_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
Bi_T_NK_cells <- MapQuery(
  anchorset = anchors_map_Bi_T_NK_cells_to_HCC_T_NK_cells,
  query = Bi_T_NK_cells,
  reference = atlas_HCC_T_NK_cells,
  refdata = list(
    to_HCC_T_NK_cells = "cell_types_from_author",
    to_HCC_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_HCC_T_NK_cells = "RNA"
  ),
  new.reduction.name = "to_HCC_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_Bi_T_NK_cells_to_HCC_T_NK_cells")
rm("atlas_HCC_T_NK_cells")
Bi_T_NK_cells[["to_HCC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(Bi_T_NK_cells,reduction = "ref.umap"),
                                                                  key = "to_HCC_T_NK_cells_",
                                                                  assay = DefaultAssay(Bi_T_NK_cells))
saveRDS(Bi_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_HCC_T_NK_cells_20250506.rds"))
rm("Bi_T_NK_cells")

# 6. map to STAD_T_NK_cells
Bi_T_NK_cells <- Bi_T_NK_cellsbubian
atlas_STAD_T_NK_cells <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells.rds"))
anchors_map_Bi_T_NK_cells_to_STAD_T_NK_cells <- FindTransferAnchors(
  reference = atlas_STAD_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = Bi_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
Bi_T_NK_cells <- MapQuery(
  anchorset = anchors_map_Bi_T_NK_cells_to_STAD_T_NK_cells,
  query = Bi_T_NK_cells,
  reference = atlas_STAD_T_NK_cells,
  refdata = list(
    to_STAD_T_NK_cells = "cell_types_from_author",
    to_STAD_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_STAD_T_NK_cells = "RNA"
  ),
  new.reduction.name = "to_STAD_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_Bi_T_NK_cells_to_STAD_T_NK_cells")
rm("atlas_STAD_T_NK_cells")
Bi_T_NK_cells[["to_STAD_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(Bi_T_NK_cells,reduction = "ref.umap"),
                                                                   key = "to_STAD_T_NK_cells_",
                                                                   assay = DefaultAssay(Bi_T_NK_cells))
saveRDS(Bi_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_STAD_T_NK_cells_20250506.rds"))
rm("Bi_T_NK_cells")

##########################################
#### 2024.02.23 plot umaps for bi dataset
###########################################
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
changenames <- c("CD4 T","CD4 T","CD8 T","NK")
names(changenames) <- c("0_T memory cells","2_T memory cells","4_T cells","6_Gamma delta T cells")
new_cluster_names <- c("Recently activated CD4 T cells", "Naive-memory CD4 T cells", "Transitional memory CD4 T cells","Cytotoxic CD8 T cells",
                       "Effector memory CD8 T cells", "Th17 cells", "NK", "Terminally exhausted CD8 T cells",
                       "Naive T cells", "Regulatory T cells", "T helper cells", "Proliferative T cells",
                       "Pre-exhausted CD8 T cells",
                       
                       
                       "CD8.c07.Temra.CX3CR1", "CD8.c09.Tk.KIR2DL4", "CD8.c05.Tem.CXCR5" ,
                       "CD8.c15.ISG.IFIT1", "CD8.c01.Tn.MAL", "CD8.c17.Tm.NME1", "CD8.c08.Tk.TYROBP" ,
                       "CD8.c12.Tex.CXCL13", "CD8.c02.Tm.IL7R", "CD8.c10.Trm.ZNF683", "CD8.c03.Tm.RPS12" ,
                       "CD8.c14.Tex.TCF7", "CD8.c06.Tem.GZMK", "CD8.c11.Tex.PDCD1", "CD8.c16.MAIT.SLC4A10" ,
                       "CD8.c04.Tm.CD52", "CD8.c13.Tex.myl12a", "CD4.c20.Treg.TNFRSF9", "CD4.c18.Treg.RTKN2" ,
                       "CD4.c06.Tm.ANXA1", "CD4.c22.ISG.IFIT1", "CD4.c01.Tn.TCF7", "CD4.c12.Tem.GZMK" ,
                       "CD4.c11.Tm.GZMA", "CD4.c13.Temra.CX3CR1", "CD4.c16.Tfh.CXCR5", "CD4.c07.Tm.ANXA2" ,
                       "CD4.c17.TfhTh1.CXCL13", "CD4.c19.Treg.S1PR1", "CD4.c10.Tm.CAPG", "CD4.c24.Mix.NME2" ,
                       "CD4.c21.Treg.OAS1", "CD4.c14.Th17.SLC4A10", "CD4.c09.Tm.CCL5", "CD4.c02.Tn.PASK" ,
                       "CD4.c03.Tn.ADSL", "CD4.c05.Tm.TNF", "CD4.c08.Tm.CREM", "CD4.c04.Tn.il7r" ,
                       "CD4.c15.Th17.IL23R", "CD4.c23.Mix.NME1",
                       
                       
                       "CD8 memory T", "CD8 exhausted cytotoxic T",
                       "CD8 pre-effector T",   "CD4 memory/effector T",
                       "NK_XCL1",   "CD4 regulatory T",
                       "CD4 naive T",   "CD8 effector T",
                       "CD4 exhausted effector T", "Cytotoxic NK",
                       
                       "Naive CD4 T cell",     "Exhausted CD8 T cell",
                       "Memory CD8 T cell",    "Treg" ,               
                       "DNT",                  "Effector CD8 T cell" ,
                       "NCAM1nFCGR3Ap NK",     "Naive CD8 T cell" ,   
                       "NKT and NK cell",      "NCAM1pFCGR3Ap NK"   ,"NCAM1pFCGR3An NK",
                       
                       "13_CD8_effector memory T", "2_NK", "2_CD8_MAIT",  "20_CD4_naive T",
                       "9_CD4_regulatory T",  "8_NK",  "35_CD4_central memory T", "18_NK",
                       "11_CD8_inter-states T", "14_CD8_TR_memory T", "32_CD4_other T", "1_CD8_cytotoxic T" ,
                       "28_NK", "51_CD8_other T",
                       
                       "LNK1", "LNK3", "LNK4", "CD8 naive T" ,
                       "CD4 naive T cell",  "LNK2", "CD4 helper T",  "Regulatory T" ,
                       "CD8 effector T cell",  "Proliferative T")

names(new_cluster_names) <- c("Recently activated CD4 T cells",   "Naive-memory CD4 T cells",         "Transitional memory CD4 T cells",
                              "Cytotoxic CD8 T cells",            "Effector memory CD8 T cells",      "Th17 cells",
                              "NK",                               "Terminally exhausted CD8 T cells", "Naive T cells",
                              "Regulatory T cells",               "T helper cells",                   "Proliferative T cells",
                              "Pre-exhausted CD8 T cells",
                              
                              "CD8.c07.Temra.CX3CR1",   "CD8.c09.Tk.KIR2DL4",
                              "CD8.c05.Tem.CXCR5",      "CD8.c15.ISG.IFIT1",      "CD8.c01.Tn.MAL",
                              "CD8.c17.Tm.NME1",        "CD8.c08.Tk.TYROBP",      "CD8.c12.Tex.CXCL13",
                              "CD8.c02.Tm.IL7R",        "CD8.c10.Trm.ZNF683",     "CD8.c03.Tm.RPS12",
                              "CD8.c14.Tex.TCF7",       "CD8.c06.Tem.GZMK",       "CD8.c11.Tex.PDCD1",
                              "CD8.c16.MAIT.SLC4A10",   "CD8.c04.Tm.CD52",        "CD8.c13.Tex.myl12a",
                              "CD4.c20.Treg.TNFRSF9",   "CD4.c18.Treg.RTKN2",     "CD4.c06.Tm.ANXA1",
                              "CD4.c22.ISG.IFIT1",      "CD4.c01.Tn.TCF7",        "CD4.c12.Tem.GZMK",
                              "CD4.c11.Tm.GZMA",        "CD4.c13.Temra.CX3CR1",   "CD4.c16.Tfh.CXCR5",
                              "CD4.c07.Tm.ANXA2",       "CD4.c17.TfhTh1.CXCL13",  "CD4.c19.Treg.S1PR1",
                              "CD4.c10.Tm.CAPG",        "CD4.c24.Mix.NME2",       "CD4.c21.Treg.OAS1",
                              "CD4.c14.Th17.SLC4A10",   "CD4.c09.Tm.CCL5",        "CD4.c02.Tn.PASK",
                              "CD4.c03.Tn.ADSL",        "CD4.c05.Tm.TNF",         "CD4.c08.Tm.CREM",
                              "CD4.c04.Tn.il7r",        "CD4.c15.Th17.IL23R",     "CD4.c23.Mix.NME1",
                              
                              
                              "C3_CD8_ZNF683",     "C1_CD8_HAVCR2",     "C2_CD8_GZMK",
                              "C6_CD4_GZMA",       "C10_NK_XCL1",       "C8_CD4_FOXP3",
                              "C5_CD4_CCR7",       "C4_CD8_CX3CR1",     "C7_CD4_CXCL13",
                              "C9_NK_FGFBP2",
                              
                              "Naive CD4 T cell",                 "Exhausted CD8 T cell",
                              "Memory CD8 T cell",                "Treg",                             "DNT",
                              "Effector CD8 T cell",              "NCAM1-FCGR3A+ NK",                 "Naive CD8 T cell",
                              "NKT and NK cell",                  "NCAM1+FCGR3A+ NK",                 "NCAM1+FCGR3A- NK",
                              
                              
                              "13_T/NK",                           "25_T/NK",                           "2_T/NK",
                              "20_T/NK",                           "9_T/NK",                            "8_T/NK",
                              "35_T/NK",                           "18_T/NK",                           "11_T/NK",
                              "14_T/NK",                           "32_T/NK",                           "1_T/NK",
                              "28_T/NK",                           "51_T/NK",                           "LNK1",
                              
                              
                              "LNK3",                             "LNK4",                             "LT2",
                              "LT3",                              "LNK2",                             "LT4",
                              "LT5",                              "LT1",                              "LPT")

Bi_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_integrated_20250506.rds"))
Bi_T_NK_cells_TICA <- subset(Bi_T_NK_cells_TICA, cells=rownames(Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_TICA <- subset(Bi_T_NK_cells_TICA, cells=rownames(Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)])
DimPlot(Bi_T_NK_cells_TICA,reduction = "umap_from_author",group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
Bi_T_NK_cells_TICA$CellType2 <- Bi_T_NK_cells_TICA$CellType
DimPlot(Bi_T_NK_cells_TICA,reduction = "umap_from_author",group.by = "CellType2",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(Bi_T_NK_cells_TICA,reduction = "umap_from_author",group.by = "predicted.to_TICA_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

Bi_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_20250506.rds"))
Bi_T_NK_cells_panT <- subset(Bi_T_NK_cells_panT, cells=rownames(Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_panT <- subset(Bi_T_NK_cells_panT, cells=rownames(Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4 <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)])
# DimPlot(Bi_T_NK_cells_panT,reduction = "umap_from_author",group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(Bi_T_NK_cells_panT,reduction = "umap_from_author",group.by = "predicted.to_pancancer_T_cells_CD8_CD4",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

Bi_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated_20250506.rds"))
Bi_T_NK_cells_panblueprint <- subset(Bi_T_NK_cells_panblueprint, cells=rownames(Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_panblueprint <- subset(Bi_T_NK_cells_panblueprint, cells=rownames(Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)])
# DimPlot(Bi_T_NK_cells_panblueprint,reduction = "umap_from_author",group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(Bi_T_NK_cells_panblueprint,reduction = "umap_from_author",group.by = "predicted.to_pancancer_blueprint_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

Bi_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_ABTC_T_NK_cells_20250506.rds"))
Bi_T_NK_cells_ABTC <- subset(Bi_T_NK_cells_ABTC, cells=rownames(Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_ABTC <- subset(Bi_T_NK_cells_ABTC, cells=rownames(Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)])
# DimPlot(Bi_T_NK_cells_ABTC,reduction = "umap_from_author",group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(Bi_T_NK_cells_ABTC,reduction = "umap_from_author",group.by = "predicted.to_ABTC_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

Bi_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_HCC_T_NK_cells_20250506.rds"))
Bi_T_NK_cells_HCC <- subset(Bi_T_NK_cells_HCC, cells=rownames(Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_HCC <- subset(Bi_T_NK_cells_HCC, cells=rownames(Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)])
# DimPlot(Bi_T_NK_cells_HCC,reduction = "umap_from_author",group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(Bi_T_NK_cells_HCC,reduction = "umap_from_author",group.by = "predicted.to_HCC_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

Bi_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_STAD_T_NK_cells_20250506.rds"))
Bi_T_NK_cells_STAD <- subset(Bi_T_NK_cells_STAD, cells=rownames(Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_STAD <- subset(Bi_T_NK_cells_STAD, cells=rownames(Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)])
# DimPlot(Bi_T_NK_cells_STAD,reduction = "umap_from_author",group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(Bi_T_NK_cells_STAD,reduction = "umap_from_author",group.by = "predicted.to_STAD_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )






##############################
#### 2024.02.23 barplot for bi dataset
#### 2025.05.06 sankey plots
##############################

rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
library(scales)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
new_cluster_names <- c("Recently activated CD4 T cells", "Naive-memory CD4 T cells", "Transitional memory CD4 T cells","Cytotoxic CD8 T cells",
                       "Effector memory CD8 T cells", "Th17 cells", "NK", "Terminally exhausted CD8 T cells",
                       "Naive T cells", "Regulatory T cells", "T helper cells", "Proliferative T cells",
                       "Pre-exhausted CD8 T cells",
                       
                       
                       "CD8.c07.Temra.CX3CR1", "CD8.c09.Tk.KIR2DL4", "CD8.c05.Tem.CXCR5" ,
                       "CD8.c15.ISG.IFIT1", "CD8.c01.Tn.MAL", "CD8.c17.Tm.NME1", "CD8.c08.Tk.TYROBP" ,
                       "CD8.c12.Tex.CXCL13", "CD8.c02.Tm.IL7R", "CD8.c10.Trm.ZNF683", "CD8.c03.Tm.RPS12" ,
                       "CD8.c14.Tex.TCF7", "CD8.c06.Tem.GZMK", "CD8.c11.Tex.PDCD1", "CD8.c16.MAIT.SLC4A10" ,
                       "CD8.c04.Tm.CD52", "CD8.c13.Tex.myl12a", "CD4.c20.Treg.TNFRSF9", "CD4.c18.Treg.RTKN2" ,
                       "CD4.c06.Tm.ANXA1", "CD4.c22.ISG.IFIT1", "CD4.c01.Tn.TCF7", "CD4.c12.Tem.GZMK" ,
                       "CD4.c11.Tm.GZMA", "CD4.c13.Temra.CX3CR1", "CD4.c16.Tfh.CXCR5", "CD4.c07.Tm.ANXA2" ,
                       "CD4.c17.TfhTh1.CXCL13", "CD4.c19.Treg.S1PR1", "CD4.c10.Tm.CAPG", "CD4.c24.Mix.NME2" ,
                       "CD4.c21.Treg.OAS1", "CD4.c14.Th17.SLC4A10", "CD4.c09.Tm.CCL5", "CD4.c02.Tn.PASK" ,
                       "CD4.c03.Tn.ADSL", "CD4.c05.Tm.TNF", "CD4.c08.Tm.CREM", "CD4.c04.Tn.il7r" ,
                       "CD4.c15.Th17.IL23R", "CD4.c23.Mix.NME1",
                       
                       
                       "CD8 memory T", "CD8 exhausted cytotoxic T",
                       "CD8 pre-effector T",   "CD4 memory/effector T",
                       "NK_XCL1",   "CD4 regulatory T",
                       "CD4 naive T",   "CD8 effector T",
                       "CD4 exhausted effector T", "Cytotoxic NK",
                       
                       "Naive CD4 T cell",     "Exhausted CD8 T cell",
                       "Memory CD8 T cell",    "Treg" ,               
                       "DNT",                  "Effector CD8 T cell" ,
                       "NCAM1nFCGR3Ap NK",     "Naive CD8 T cell" ,   
                       "NKT and NK cell",      "NCAM1pFCGR3Ap NK"   ,"NCAM1pFCGR3An NK",
                       
                       "13_CD8_effector memory T", "2_NK", "2_CD8_MAIT",  "20_CD4_naive T",
                       "9_CD4_regulatory T",  "8_NK",  "35_CD4_central memory T", "18_NK",
                       "11_CD8_inter-states T", "14_CD8_TR_memory T", "32_CD4_other T", "1_CD8_cytotoxic T" ,
                       "28_NK", "51_CD8_other T",
                       
                       "LNK1", "LNK3", "LNK4", "CD8 naive T" ,
                       "CD4 naive T cell",  "LNK2", "CD4 helper T",  "Regulatory T" ,
                       "CD8 effector T cell",  "Proliferative T")

names(new_cluster_names) <- c("Recently activated CD4 T cells",   "Naive-memory CD4 T cells",         "Transitional memory CD4 T cells",
                              "Cytotoxic CD8 T cells",            "Effector memory CD8 T cells",      "Th17 cells",
                              "NK",                               "Terminally exhausted CD8 T cells", "Naive T cells",
                              "Regulatory T cells",               "T helper cells",                   "Proliferative T cells",
                              "Pre-exhausted CD8 T cells",
                              
                              "CD8.c07.Temra.CX3CR1",   "CD8.c09.Tk.KIR2DL4",
                              "CD8.c05.Tem.CXCR5",      "CD8.c15.ISG.IFIT1",      "CD8.c01.Tn.MAL",
                              "CD8.c17.Tm.NME1",        "CD8.c08.Tk.TYROBP",      "CD8.c12.Tex.CXCL13",
                              "CD8.c02.Tm.IL7R",        "CD8.c10.Trm.ZNF683",     "CD8.c03.Tm.RPS12",
                              "CD8.c14.Tex.TCF7",       "CD8.c06.Tem.GZMK",       "CD8.c11.Tex.PDCD1",
                              "CD8.c16.MAIT.SLC4A10",   "CD8.c04.Tm.CD52",        "CD8.c13.Tex.myl12a",
                              "CD4.c20.Treg.TNFRSF9",   "CD4.c18.Treg.RTKN2",     "CD4.c06.Tm.ANXA1",
                              "CD4.c22.ISG.IFIT1",      "CD4.c01.Tn.TCF7",        "CD4.c12.Tem.GZMK",
                              "CD4.c11.Tm.GZMA",        "CD4.c13.Temra.CX3CR1",   "CD4.c16.Tfh.CXCR5",
                              "CD4.c07.Tm.ANXA2",       "CD4.c17.TfhTh1.CXCL13",  "CD4.c19.Treg.S1PR1",
                              "CD4.c10.Tm.CAPG",        "CD4.c24.Mix.NME2",       "CD4.c21.Treg.OAS1",
                              "CD4.c14.Th17.SLC4A10",   "CD4.c09.Tm.CCL5",        "CD4.c02.Tn.PASK",
                              "CD4.c03.Tn.ADSL",        "CD4.c05.Tm.TNF",         "CD4.c08.Tm.CREM",
                              "CD4.c04.Tn.il7r",        "CD4.c15.Th17.IL23R",     "CD4.c23.Mix.NME1",
                              
                              
                              "C3_CD8_ZNF683",     "C1_CD8_HAVCR2",     "C2_CD8_GZMK",
                              "C6_CD4_GZMA",       "C10_NK_XCL1",       "C8_CD4_FOXP3",
                              "C5_CD4_CCR7",       "C4_CD8_CX3CR1",     "C7_CD4_CXCL13",
                              "C9_NK_FGFBP2",
                              
                              "Naive CD4 T cell",                 "Exhausted CD8 T cell",
                              "Memory CD8 T cell",                "Treg",                             "DNT",
                              "Effector CD8 T cell",              "NCAM1-FCGR3A+ NK",                 "Naive CD8 T cell",
                              "NKT and NK cell",                  "NCAM1+FCGR3A+ NK",                 "NCAM1+FCGR3A- NK",
                              
                              
                              "13_T/NK",                           "25_T/NK",                           "2_T/NK",
                              "20_T/NK",                           "9_T/NK",                            "8_T/NK",
                              "35_T/NK",                           "18_T/NK",                           "11_T/NK",
                              "14_T/NK",                           "32_T/NK",                           "1_T/NK",
                              "28_T/NK",                           "51_T/NK",                           "LNK1",
                              
                              
                              "LNK3",                             "LNK4",                             "LT2",
                              "LT3",                              "LNK2",                             "LT4",
                              "LT5",                              "LT1",                              "LPT")

Bi_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_integrated_20250506.rds"))
Bi_T_NK_cells_TICA <- subset(Bi_T_NK_cells_TICA, cells=rownames(Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_TICA <- subset(Bi_T_NK_cells_TICA, cells=rownames(Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_TICA@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)])
cell_summary_TICA <- data.frame(Bi_T_NK_cells_TICA$cell_types_from_author,paste0("TICA: ",Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells))
# rm("Bi_T_NK_cells_TICA")
Bi_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_20250506.rds"))
Bi_T_NK_cells_panT <- subset(Bi_T_NK_cells_panT, cells=rownames(Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_panT <- subset(Bi_T_NK_cells_panT, cells=rownames(Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_panT@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4 <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)])
cell_summary_panT <- data.frame(Bi_T_NK_cells_panT$cell_types_from_author,paste0("pan- T: ",Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4))
# rm("Bi_T_NK_cells_panT")
Bi_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated_20250506.rds"))
Bi_T_NK_cells_panblueprint <- subset(Bi_T_NK_cells_panblueprint, cells=rownames(Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_panblueprint <- subset(Bi_T_NK_cells_panblueprint, cells=rownames(Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_panblueprint@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)])
cell_summary_panblueprint <- data.frame(Bi_T_NK_cells_panblueprint$cell_types_from_author,paste0("pan- blueprint: ",Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))
# rm("Bi_T_NK_cells_panblueprint")
Bi_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_ABTC_T_NK_cells_20250506.rds"))
Bi_T_NK_cells_ABTC <- subset(Bi_T_NK_cells_ABTC, cells=rownames(Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_ABTC <- subset(Bi_T_NK_cells_ABTC, cells=rownames(Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_ABTC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)])
cell_summary_ABTC <- data.frame(Bi_T_NK_cells_ABTC$cell_types_from_author,paste0("ABTC: ",Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells))
# rm("Bi_T_NK_cells_ABTC")
Bi_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_HCC_T_NK_cells_20250506.rds"))
Bi_T_NK_cells_HCC <- subset(Bi_T_NK_cells_HCC, cells=rownames(Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_HCC <- subset(Bi_T_NK_cells_HCC, cells=rownames(Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_HCC@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)])
cell_summary_HCC <- data.frame(Bi_T_NK_cells_HCC$cell_types_from_author,paste0("HCC: ",Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells))
# rm("Bi_T_NK_cells_HCC")
Bi_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_STAD_T_NK_cells_20250506.rds"))
Bi_T_NK_cells_STAD <- subset(Bi_T_NK_cells_STAD, cells=rownames(Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_1"] >0 & Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] < 0,]),invert = T)
Bi_T_NK_cells_STAD <- subset(Bi_T_NK_cells_STAD, cells=rownames(Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings)[Bi_T_NK_cells_STAD@reductions$umap_from_author@cell.embeddings[,"umapfromauthor_2"] > 5],invert = T)
Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells <- as.character(new_cluster_names[as.character(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)])
cell_summary_STAD <- data.frame(Bi_T_NK_cells_STAD$cell_types_from_author,paste0("STAD: ",Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells))
# rm("Bi_T_NK_cells_STAD")

cell_summary <- cbind(cell_summary_TICA,to_panT=cell_summary_panT[,2],
                      to_panblueprint=cell_summary_panblueprint[,2],to_ABTC=cell_summary_ABTC[,2],to_HCC=cell_summary_HCC[,2],
                      to_STAD=cell_summary_STAD[,2])
colnames(cell_summary)[1:2] <-c("cell_types_from_Bi","to_TICA")
# cell_summary[1:3,]
# cell_summary[,"cell_types_from_Bi"] <- as.character(cell_summary[,"cell_types_from_Bi"])
# cell_summary[,"hpca_fine"] <- as.character(cell_summary[,"hpca_fine"])
unique(cell_summary$cell_types_from_Bi)
# [1] "41BB-Hi CD8+ T cell" "41BB-Lo CD8+ T cell" "T-Reg"               "NKT"                 "Effector T-Helper"   "FGFBP2- NK"         
# [7] "Cycling CD8+ T cell" "FGFBP2+ NK"          "MX1-Hi CD8+ T cell"  "Memory T-Helper"    




library(networkD3)
library(ggalluvial)
# cell_summary[,"cell_types_from_Bi"] <- as.character(cell_summary[,"cell_types_from_Bi"])
colnames_atlas <- colnames(cell_summary)
my_colors <- c(
  "Cycling CD8+ T cell"         = "#E41A1C",
  "T-Reg"     = "#377EB8",
  "Effector T-Helper"   = "#4DAF4A",
  "Memory T-Helper"      = "#FF7F00",
  "41BB-Lo CD8+ T cell"      = "#984EA3",
  "41BB-Hi CD8+ T cell"        = "#00CED1",
  "MX1-Hi CD8+ T cell" = "#A65628",
  "NKT"   = "#FFD700",
  "FGFBP2- NK" = "#999999",
  "FGFBP2+ NK" = "#66C2A5"
)
  ### For TICA
  links <- dplyr::count(cell_summary,cell_types_from_Bi,to_TICA)
  sankeyCols <- c("source", "target", "value")
  colnames(links) <- sankeyCols
  links <- links[links$value >10,]
  links$source <- factor(links$source,levels = c("Cycling CD8+ T cell","T-Reg","Effector T-Helper","Memory T-Helper",
                                                 "41BB-Lo CD8+ T cell","41BB-Hi CD8+ T cell","MX1-Hi CD8+ T cell","NKT",
                                               "FGFBP2- NK","FGFBP2+ NK"))
  links$target <- factor(links$target,levels = c("TICA: Proliferative T cells",
                                                 "TICA: Regulatory T cells","TICA: T helper cells","TICA: Naive-memory CD4 T cells","TICA: Th17 cells","TICA: Transitional memory CD4 T cells","TICA: Naive T cells",
                                                 "TICA: Cytotoxic CD8 T cells","TICA: Recently activated CD4 T cells","TICA: Effector memory CD8 T cells","TICA: Pre-exhausted CD8 T cells","TICA: Terminally exhausted CD8 T cells",
                                                 "TICA: NK" ))

  # pdf(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.pdf"),width=10)
  # png(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.png"),width=1800)
  p1 <- ggplot(data = links,
         aes(axis1 = source, axis2 = target, y = value)) +
    geom_alluvium(aes(fill = source),width = 1/10) +
    geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
    # geom_bar() +
    # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
    scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
    scale_fill_manual(values=my_colors) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
      aes(label = as.character(target)),
      stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
      aes(label = as.character(source)),
      stat = "stratum", size = 6, direction = "y", nudge_x = -0.4) +
    theme(legend.position = "none")
  # dev.off()
  
  ### For pan-T
  links <- dplyr::count(cell_summary,cell_types_from_Bi,to_panT)
  sankeyCols <- c("source", "target", "value")
  colnames(links) <- sankeyCols
  # links <- links[links$value >10,]
  links$source <- factor(links$source,levels = c("Cycling CD8+ T cell","T-Reg","Effector T-Helper","Memory T-Helper",
                                                 "41BB-Lo CD8+ T cell","41BB-Hi CD8+ T cell","MX1-Hi CD8+ T cell","NKT",
                                                 "FGFBP2- NK","FGFBP2+ NK"))
  links$target <- factor(links$target,levels = c("pan- T: CD4.c10.Tm.CAPG","pan- T: CD8.c10.Trm.ZNF683","pan- T: CD4.c03.Tn.ADSL","pan- T: CD4.c06.Tm.ANXA1","pan- T: CD8.c01.Tn.MAL",      
                                                 "pan- T: CD8.c07.Temra.CX3CR1","pan- T: CD8.c17.Tm.NME1","pan- T: CD4.c07.Tm.ANXA2","pan- T: CD4.c20.Treg.TNFRSF9" ))
  
  # pdf(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.pdf"),width=10)
  # png(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.png"),width=1800)
  p_nouse <- ggplot(data = links,
               aes(axis1 = source, axis2 = target, y = value)) +
    geom_alluvium(aes(fill = source),width = 1/10) +
    geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
    # geom_bar() +
    # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
    scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
    scale_fill_manual(values=my_colors) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
      aes(label = as.character(target)),
      stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
      aes(label = as.character(source)),
      stat = "stratum", size = 6, direction = "y", nudge_x = -0.4) +
    theme(legend.position = "none")
  # dev.off()
  
  ## For blueprint
  links <- dplyr::count(cell_summary,cell_types_from_Bi,to_panblueprint)
  sankeyCols <- c("source", "target", "value")
  colnames(links) <- sankeyCols
  links <- links[links$value >20,]
  links$source <- factor(links$source,levels = c("Cycling CD8+ T cell","T-Reg","Effector T-Helper","Memory T-Helper",
                                                 "41BB-Lo CD8+ T cell","41BB-Hi CD8+ T cell","MX1-Hi CD8+ T cell","NKT",
                                                 "FGFBP2- NK","FGFBP2+ NK"))
  links$target <- factor(links$target,levels = c("pan- blueprint: CD4 regulatory T","pan- blueprint: CD4 memory/effector T","pan- blueprint: CD4 exhausted effector T","pan- blueprint: CD4 naive T",
                                                 "pan- blueprint: CD8 pre-effector T","pan- blueprint: CD8 exhausted cytotoxic T","pan- blueprint: CD8 memory T","pan- blueprint: CD8 effector T",
                                                 "pan- blueprint: NK_XCL1","pan- blueprint: Cytotoxic NK"))
  # pdf(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_flow_20250506.pdf"),width=10)
  p2 <- ggplot(data = links,
         aes(axis1 = source, axis2 = target, y = value)) +
    geom_alluvium(aes(fill = source),width = 1/10) +
    geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
    # geom_bar() +
    # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
    scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
    scale_fill_manual(values=my_colors) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
      aes(label = as.character(target)),
      stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
      aes(label = as.character(source)),
      stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
    theme(legend.position = "none")
  # dev.off()


  ## For ABTC
  links <- dplyr::count(cell_summary,cell_types_from_Bi,to_ABTC)
  sankeyCols <- c("source", "target", "value")
  colnames(links) <- sankeyCols
  links <- links[links$value >10,]
  links$source <- factor(links$source,levels = c("Cycling CD8+ T cell","T-Reg","Effector T-Helper","Memory T-Helper",
                                                 "41BB-Lo CD8+ T cell","41BB-Hi CD8+ T cell","MX1-Hi CD8+ T cell","NKT",
                                                 "FGFBP2- NK","FGFBP2+ NK"))
  links$target <- factor(links$target,levels = c("ABTC: Treg","ABTC: Naive CD4 T cell","ABTC: DNT",
                                                 "ABTC: Exhausted CD8 T cell","ABTC: Effector CD8 T cell","ABTC: Memory CD8 T cell","ABTC: Naive CD8 T cell",
                                                 "ABTC: NCAM1pFCGR3An NK","ABTC: NCAM1pFCGR3Ap NK","ABTC: NCAM1nFCGR3Ap NK" ))
  # pdf(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_ABTC_T_NK_cells_flow_20250506.pdf"),width=10)
  p3 <- ggplot(data = links,
         aes(axis1 = source, axis2 = target, y = value)) +
    geom_alluvium(aes(fill = source),width = 1/10) +
    geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
    # geom_bar() +
    # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
    scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
    scale_fill_manual(values=my_colors) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
      aes(label = as.character(target)),
      stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
      aes(label = as.character(source)),
      stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
    theme(legend.position = "none")
  # dev.off()
  
  ## For HCC
  links <- dplyr::count(cell_summary,cell_types_from_Bi,to_HCC)
  sankeyCols <- c("source", "target", "value")
  colnames(links) <- sankeyCols
  links <- links[links$value >10,]
  links$source <- factor(links$source,levels = c("Cycling CD8+ T cell","T-Reg","Effector T-Helper","Memory T-Helper",
                                                 "41BB-Lo CD8+ T cell","41BB-Hi CD8+ T cell","MX1-Hi CD8+ T cell","NKT",
                                                 "FGFBP2- NK","FGFBP2+ NK"))
  links$target <- factor(links$target,levels = c("HCC: 32_CD4_other T","HCC: 51_CD8_other T",
                                                 "HCC: 9_CD4_regulatory T","HCC: 20_CD4_naive T","HCC: 35_CD4_central memory T",
                                                 "HCC: 28_NK","HCC: 13_CD8_effector memory T","HCC: 11_CD8_inter-states T","HCC: 1_CD8_cytotoxic T","HCC: 2_CD8_MAIT","HCC: 14_CD8_TR_memory T",
                                                 "HCC: 8_NK","HCC: 25_NK", "HCC: 18_NK"))
  # pdf(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_HCC_T_NK_cells_flow_20250506.pdf"),width=10,height=4)
  p4 <- ggplot(data = links,
         aes(axis1 = source, axis2 = target, y = value)) +
    geom_alluvium(aes(fill = source),width = 1/10) +
    geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
    # geom_bar() +
    # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
    scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
    scale_fill_manual(values=my_colors) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
      aes(label = as.character(target)),
      stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
      aes(label = as.character(source)),
      stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
    theme(legend.position = "none")
  # dev.off()

  ## For STAD
  links <- dplyr::count(cell_summary,cell_types_from_Bi,to_STAD)
  sankeyCols <- c("source", "target", "value")
  colnames(links) <- sankeyCols
  links <- links[links$value >10,]
  links$source <- factor(links$source,levels = c("Cycling CD8+ T cell","T-Reg","Effector T-Helper","Memory T-Helper",
                                                 "41BB-Lo CD8+ T cell","41BB-Hi CD8+ T cell","MX1-Hi CD8+ T cell","NKT",
                                                 "FGFBP2- NK","FGFBP2+ NK"))
  links$target <- factor(links$target,levels = c("STAD: Proliferative T",
                                                 "STAD: Regulatory T","STAD: CD4 naive T cell","STAD: CD8 effector T cell","STAD: CD4 helper T",
                                                 "STAD: LNK2","STAD: LNK3","STAD: CD8 naive T",
                                                 "STAD: LNK1","STAD: LNK4"))
  # pdf(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_STAD_T_NK_cells_flow_20250506.pdf"),width=15,height = 6)
  p5 <- ggplot(data = links,
         aes(axis1 = source, axis2 = target, y = value)) +
    geom_alluvium(aes(fill = source),width = 1/10) +
    geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
    # geom_bar() +
    # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
    scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
    scale_fill_manual(values=my_colors) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
      aes(label = as.character(target)),
      stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
    ggrepel::geom_text_repel(
      # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
      aes(label = as.character(source)),
      stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
    theme(legend.position = "none")
  # dev.off()

pdf(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_allatlases_T_NK_cells_flow_20250506.pdf"),height=25,width=15)
print(CombinePlots(list(p1,p2,p3,p4,p5),ncol=1))
dev.off()

#### calculate the ARI
adjustedRandIndex(cell_summary$cell_types_from_Bi,cell_summary$to_TICA)
# [1] 0.3168938
adjustedRandIndex(cell_summary$cell_types_from_Bi,cell_summary$to_panblueprint)
# [1] 0.2950238
adjustedRandIndex(cell_summary$cell_types_from_Bi,cell_summary$to_ABTC)
# [1] 0.3048876c
adjustedRandIndex(cell_summary$cell_types_from_Bi,cell_summary$to_HCC)
# [1] 0.2375804
adjustedRandIndex(cell_summary$cell_types_from_Bi,cell_summary$to_STAD)
# [1] 0.2454156




#   ## previouse barplots
#   ##
#   ##
# index <- c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD")
# index2 <- c("TICA: ","pan- T: ","pan- blueprint: ","ABTC: ","HCC: ","STAD: ")
# # for T_cell:CD4+_effector_memory
# output0 <- NULL
# output1 <- NULL;output2 <- NULL;output3<-NULL
# pp <- NULL
# # for (i in 1:length(index)) {
#   i = 1
#   #### for p1
#   
#   cell_summary_subset1 <- subset(cell_summary,cell_types_from_Bi %in% c("T-Reg","Effector T-Helper","Memory T-Helper")) #
#   cell_summary_subset1$hpca_fine <- "CD4+ T"
#   # cell_summary_subset2 <- subset(cell_summary,cell_types_from_Bi %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
#   # cell_summary_subset3 <- subset(cell_summary,cell_types_from_Bi %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
#   cell_summary_subset <- cell_summary_subset1 #rbind(cell_summary_subset1,cell_summary_subset2,cell_summary_subset3)
#   output_pre1 <- as.data.frame(table(cell_summary_subset1[,c(index[i],"hpca_fine")]))
#   output_pre1 <- data.frame(group=index[i],output_pre1)
#   output_pre1[,2] <- as.character(output_pre1[,2])
#   output_pre1$Freq <- as.numeric(output_pre1$Freq)
#   # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
#   colnames(output_pre1)[2] <- "fill"
#   # output_pre1$fill <- as.character(new_cluster_names[output_pre1$fill])
#   output_pre1 <- data.frame(output_pre1,populations=output_pre1$fill)
#   output_pre1$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre1$fill),": ",2)))[,2])
#   output_pre1$Freq <- as.numeric(output_pre1$Freq)
#   output_pre1 <- output_pre1[order(output_pre1$Freq),]
#   output_pre1 <- output_pre1 %>%
#     group_by(group) %>%
#     mutate(label_y = cumsum(Freq)-0.5*Freq)
#   output_pre1 <- output_pre1[order(-output_pre1$Freq),]
#   output1 <- rbind(output1,output_pre1)
#   
#   output1 <- as.data.frame(output1)
#   output1 <- rbind(output1,c("RF","CD4+ T","CD4+ T",sum(output_pre1$Freq),"RF: CD4+ T",sum(output_pre1$Freq)/2))#;output1 <- output
#   # output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2));output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output1,output2,output3)
#   output1$fill <- factor(output1$fill,levels = c(output1$fill))
#   output1$populations <- factor(output1$populations,levels = c(output1$populations))
#   output1$group <- factor(output1$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output1$Freq <- as.numeric(output1$Freq)
#   output1$label_y <- as.numeric(output1$label_y)
#   output1[output1$Freq < 3,"fill"] <- "" #10
#   if (i ==1) {
#     cols <- hue_pal()(length(unique(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#     names(cols) <- paste0(index2[i],sort(unique(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   }
#   if (i ==2) {
#     cols <- hue_pal()(length(unique(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#     names(cols) <- paste0(index2[i],sort(unique(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   }
#   if (i ==3) {
#     cols <- hue_pal()(length(unique(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
#     # names(cols) <- sort(as.character(new_cluster_names[paste0("pan- blueprint: ",unique(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))]))
#     names(cols) <- sort(as.character(paste0("pan- blueprint: ",unique(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))))
#   }
#   if (i ==4) {
#     cols <- hue_pal()(length(unique(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
#     # names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
#     names(cols) <- sort(as.character(paste0("ABTC: ",unique(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells))))
#   }
#   if (i ==5) {
#     cols <- hue_pal()(length(unique(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
#     # names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
#     names(cols) <- as.character(paste0("HCC: ",sort(unique(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells))))
#   }
#   if (i ==6) {
#     cols <- hue_pal()(length(unique(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
#     # names(cols) <- sort(as.character(new_cluster_names[paste0("STAD: ",unique(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells))]))
#     names(cols) <- sort(as.character(paste0("STAD: ",unique(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells))))
#   }
#   
#   
#   
#   cols <- c("RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
#   p1 <- ggplot(output1,aes(x=group,y=Freq,fill=populations)) + 
#     geom_bar(stat="identity",colour="black") + 
#     scale_fill_manual(values=cols) +
#     geom_text(aes(y=label_y,label=fill),size =10) +
#     theme_classic()
#   
#   
#   ##### for p2
#   
#   # cell_summary_subset2 <- subset(cell_summary,cell_types_from_Bi %in% c("0_T memory cells","2_T memory cells")) #
#   # cell_summary_subset2$hpca_fine <- "CD4+ T"
#   cell_summary_subset2 <- subset(cell_summary,cell_types_from_Bi %in% c("41BB-Hi CD8+ T cell","41BB-Lo CD8+ T cell","Cycling CD8+ T cell","MX1-Hi CD8+ T cell"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
#   # cell_summary_subset3 <- subset(cell_summary,cell_types_from_Bi %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
#   # cell_summary_subset <- cell_summary_subset2 #rbind(cell_summary_subset2,cell_summary_subset2,cell_summary_subset3)
#   output_pre2 <- as.data.frame(table(cell_summary_subset2[,c(index[i],"hpca_fine")]))
#   output_pre2 <- data.frame(group=index[i],output_pre2)
#   output_pre2[,2] <- as.character(output_pre2[,2])
#   output_pre2$Freq <- as.numeric(output_pre2$Freq)
#   # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
#   colnames(output_pre2)[2] <- "fill"
#   # output_pre2$fill <- as.character(new_cluster_names[output_pre2$fill])
#   output_pre2 <- data.frame(output_pre2,populations=output_pre2$fill)
#   output_pre2$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre2$fill),": ",2)))[,2])
#   output_pre2$Freq <- as.numeric(output_pre2$Freq)
#   output_pre2 <- output_pre2[order(output_pre2$Freq),]
#   output_pre2 <- output_pre2 %>%
#     group_by(group) %>%
#     mutate(label_y = cumsum(Freq)-0.5*Freq)
#   output_pre2 <- output_pre2[order(-output_pre2$Freq),]
#   output2 <- rbind(output2,output_pre2)
#   
#   output2 <- as.data.frame(output2)
#   # output2 <- rbind(output2,c("RF","CD4+ T","CD4+ T",sum(output_pre2$Freq),"RF: CD4+ T",sum(output_pre2$Freq)/2))#;output2 <- output
#   output2 <- rbind(output2,c("RF","CD8+ T","CD8+ T",sum(output_pre2$Freq),"RF: CD8+ T",sum(output_pre2$Freq)/2))# ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2));output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output2,output2,output3)
#   output2$fill <- factor(output2$fill,levels = c(output2$fill))
#   output2$populations <- factor(output2$populations,levels = c(output2$populations))
#   output2$group <- factor(output2$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output2$Freq <- as.numeric(output2$Freq)
#   output2$label_y <- as.numeric(output2$label_y)
#   output2[output2$Freq < 3,"fill"] <- "" #20
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # names(cols) <- paste0(index2[i],sort(unique(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # names(cols) <- paste0(index2[i],sort(unique(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("pan- blueprint: ",sort(unique(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("STAD: ",sort(unique(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))])
#   
#   
#   
#   cols <- c("RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
#   p2 <- ggplot(output2,aes(x=group,y=Freq,fill=populations)) + 
#     geom_bar(stat="identity",colour="black") + 
#     scale_fill_manual(values=cols) +
#     geom_text(aes(y=label_y,label=fill),size =10) +
#     theme_classic()
#   
#   
#   #### for p3
#   
#   # cell_summary_subset3 <- subset(cell_summary,cell_types_from_Bi %in% c("0_T memory cells","2_T memory cells")) #
#   # cell_summary_subset3$hpca_fine <- "CD4+ T"
#   # cell_summary_subset2 <- subset(cell_summary,cell_types_from_Bi %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
#   cell_summary_subset3 <- subset(cell_summary,cell_types_from_Bi %in% c("NKT"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
#   # cell_summary_subset <- cell_summary_subset3 #rbind(cell_summary_subset3,cell_summary_subset2,cell_summary_subset3)
#   output_pre3 <- as.data.frame(table(cell_summary_subset3[,c(index[i],"hpca_fine")]))
#   output_pre3 <- data.frame(group=index[i],output_pre3)
#   output_pre3[,2] <- as.character(output_pre3[,2])
#   output_pre3$Freq <- as.numeric(output_pre3$Freq)
#   # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
#   colnames(output_pre3)[2] <- "fill"
#   # output_pre3$fill <- as.character(new_cluster_names[output_pre3$fill])
#   output_pre3 <- data.frame(output_pre3,populations=output_pre3$fill)
#   output_pre3$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre3$fill),": ",2)))[,2])
#   output_pre3$Freq <- as.numeric(output_pre3$Freq)
#   output_pre3 <- output_pre3[order(output_pre3$Freq),]
#   output_pre3 <- output_pre3 %>%
#     group_by(group) %>%
#     mutate(label_y = cumsum(Freq)-0.5*Freq)
#   output_pre3 <- output_pre3[order(-output_pre3$Freq),]
#   output3 <- rbind(output3,output_pre3)
#   
#   output3 <- as.data.frame(output3)
#   # output3 <- rbind(output3,c("RF","CD4+ T","CD4+ T",sum(output_pre3$Freq),"RF: CD4+ T",sum(output_pre3$Freq)/2))#;output3 <- output
#   # output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   output3 <- rbind(output3,c("RF","NK cell","NK",sum(output_pre3$Freq),"RF: NK",sum(output_pre3$Freq)/2))#;output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output3,output2,output3)
#   output3$fill <- factor(output3$fill,levels = c(output3$fill))
#   output3$populations <- factor(output3$populations,levels = c(output3$populations))
#   output3$group <- factor(output3$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output3$Freq <- as.numeric(output3$Freq)
#   output3$label_y <- as.numeric(output3$label_y)
#   output3[output3$Freq < 3,"fill"] <- "" #30
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # names(cols) <- paste0(index2[i],sort(unique(Bi_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # names(cols) <- paste0(index2[i],sort(unique(Bi_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("pan- blueprint: ",sort(unique(Bi_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(Bi_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(Bi_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("STAD: ",sort(unique(Bi_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))])
#   
#   
#   cols <- c("RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
#   p3 <- ggplot(output3,aes(x=group,y=Freq,fill=populations)) + 
#     geom_bar(stat="identity",colour="black") + 
#     scale_fill_manual(values=cols) +
#     geom_text(aes(y=label_y,label=fill),size =10) +
#     theme_classic()
#   
#   
#   # CombinePlots(plots=list(p1,p2,p3),nrow=1,legend="right")
#   
#   p1 <- p1 + theme(axis.text.x = element_text(size = 34,color = "black"),axis.text.y = element_text(size = 28,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="")
#   p2 <- p2 + theme(axis.text.x = element_text(size = 34,color = "black"),axis.text.y = element_text(size = 28,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="")
#   p3 <- p3 + theme(axis.text.x = element_text(size = 34,color = "black"),axis.text.y = element_text(size = 28,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="")
#   CombinePlots(plots=list(p1,p2,p3),nrow=1,legend="none")
#   # pdf(paste0(workdir,"test.pdf"))
#   # pp <- list(pp,CombinePlots(plots=list(p1,p2,p3),nrow=1,legend="right"))
#   # dev.off()
# }





####################################################
#### 2024.02.28 using Livnat Jerby-Arnon
####
#### 2025.05.05 further cluster 0,1,7,10 clusters
####################################################
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
"%ni%" <- Negate("%in%")



workdir <- "/Users/jingyang/Dropbox/work/Collabrator_Justin_Balko/Livnat Jerby-Arnon_melanoma/data/"
malignant <- readRDS(paste0(workdir,"ImmRes_Rfiles/scData/Mel.malignant.rds"))
non_malignant <- readRDS(paste0(workdir,"ImmRes_Rfiles/scData/Mel.all.data.QC.rds"))
SCLC <- CreateSeuratObject(counts = non_malignant$tpm, min.cells = 5, min.features = 10)
SCLC[["percent.mt"]] <- PercentageFeatureSet(SCLC, pattern = Mtpattern)
print(VlnPlot(SCLC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
plot1 <- FeatureScatter(SCLC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SCLC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(CombinePlots(plots = list(plot1, plot2)))
plot(sort(SCLC@meta.data$nCount_RNA,decreasing = T),ylab="UMI counts",lty=2,pch=16,log="xy")
abline(h=1000)
abline(h=500)
# SCLC <- NormalizeData(SCLC)
SCLC <- FindVariableFeatures(SCLC, selection.method = "vst", nfeatures = 2000)
SCLC <- ScaleData(SCLC)
SCLC <- RunPCA(SCLC)#, features = var.genes)
SCLC <- FindNeighbors(SCLC, dims = 1:30)
SCLC <- FindClusters(SCLC, resolution = 2)
SCLC <- RunUMAP(SCLC, dims = 1:30)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
SCLC <- CellCycleScoring(SCLC,s.genes,g2m.genes)
print(FeaturePlot(SCLC, features=c("percent.mt","nFeature_RNA","nCount_RNA","PC_1"),label=T))
plot1<- ggplot(data=data.frame(cluster=SCLC@active.ident,nUMI=SCLC@meta.data$nCount_RNA),aes(x=cluster,y=nUMI))+geom_boxplot()+scale_y_log10()   
plot2<- ggplot(data=data.frame(cluster=SCLC@active.ident,ngene=SCLC@meta.data$nFeature_RNA),aes(x=cluster,y=ngene))+geom_boxplot()+scale_y_log10()
plot3<-ggplot(data=data.frame(cluster=SCLC@active.ident,MtRNA=SCLC@meta.data$percent.mt),aes(x=cluster,y=MtRNA))+geom_boxplot()  
plot4<- ggplot(data=data.frame(cluster=SCLC@active.ident),aes(x=cluster))+geom_bar(stat="count")+theme(axis.text.x = element_text(angle = 45))+ylab("nCells")+ggtitle(paste0("total cells:",dim(SCLC)[2])) 
plot5<- DimPlot(SCLC, reduction = "pca",label=T,label.size=4) 
plot6<-DimPlot(SCLC, reduction = "umap",label=T,label.size=4)
print(CombinePlots(plots = list(plot1, plot2, plot3,plot4,plot5,plot6),rel_heights=c(1,1)))
markerfile <- "/Users/jingyang/Dropbox/work/single-cell RNA-seq 10x pipeline/scRNAseq_scripts/notations/PanglaoDB_markers_21_Oct_2019.tsv"
species <- "Hs"
add_Subtype=FALSE
celltype_predictmethod="cta"
marker<-data.frame(fread(markerfile))
species_ind<-regexpr(species,marker[,1])
marker_species<-marker[species_ind>0 & marker$ubiquitousness.index<0.05,]
##change the gene symbol only keep the first letter capitalize
if (species=="Mm") {
  marker_species$official.gene.symbol<-paste0(substr(marker_species$official.gene.symbol,1,1),substr(tolower(marker_species$official.gene.symbol),2,nchar(marker_species$official.gene.symbol)))
}
##
cellType<-tapply(marker_species$official.gene.symbol,marker_species$cell.type,list)
medianexp<-t(apply(GetAssayData(SCLC,slot="data"),1,function(x){tapply(x,SCLC@active.ident,mean)}))
if (nrow(medianexp)==1) medianexp<-t(medianexp)
source("markerCode_filter.R")
predict_celltype<-ORA_celltype(medianexp,cellType,method=celltype_predictmethod)
new.cluster.ids<-paste(levels(SCLC@active.ident),rownames(predict_celltype$predict_result),sep="_")

names(new.cluster.ids) <- levels(SCLC)
SCLC <- RenameIdents(SCLC, new.cluster.ids)
print(DimPlot(SCLC, reduction = "umap",label=T,label.size=6))

SCLC.markers <- FindAllMarkers(SCLC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SCLC@misc$markers<-SCLC.markers

top10 <- SCLC.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
if (nrow(top10)>200) {genesize=5} else {
  if (nrow(top10)>100) {genesize=6}  else {genesize=7} }

print(DoHeatmap(SCLC, features = top10$gene)+ theme(axis.text.y = element_text(size = genesize)))
print(DoHeatmap(SCLC, features = top10$gene,slot="data") + scale_fill_viridis_b() + theme(axis.text.y = element_text(size = genesize))) 

SCLC[["tSNE"]] <- CreateDimReducObject(embeddings=non_malignant$tsne,key="tSNE_",assay="RNA")
SCLC$cell_types_from_author <- non_malignant$cell.types
saveRDS(SCLC,file="/Users/jingyang/Dropbox/work/collabrator_Chi-Yan/reference_1_a_cancer_cell_program_promotes_T_cell_exclusion_and_resistance_to_checkpoint_blockade/ref1_data_20250506.rds")



# load("/Users/jingyang/Dropbox/work/cancer-immu for biological discoveries/code/scRNA-seq_from_Qi_object.Rdata")
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

LJA <- readRDS("/Users/jingyang/Dropbox/work/collabrator_Chi-Yan/reference_1_a_cancer_cell_program_promotes_T_cell_exclusion_and_resistance_to_checkpoint_blockade/ref1_data_20250506.rds")
# LJA <- subset(LJA, idents=c("0_T cells","1_T memory cells","7_T cells","10_NK cells"))
LJA <- subset(LJA, idents=c("0_T cells","1_T memory cells","3_T cells","4_T memory cells","5_T cells","7_T cells","10_T cells","14_NK cells"))
LJA <- RenameIdents(object = LJA, `0_T cells` = "Tgd_c3")
LJA <- RenameIdents(object = LJA, `1_T memory cells` = "CD4_c2_Tn")
LJA <- RenameIdents(object = LJA, `3_T cells` = "CD8_c1_Tex")
LJA <- RenameIdents(object = LJA, `4_T memory cells` = "CD4_c1_Treg")
LJA <- RenameIdents(object = LJA, `5_T cells` = "CD8_c2_Teff")
LJA <- RenameIdents(object = LJA, `7_T cells` = "NKT_c0")
LJA <- RenameIdents(object = LJA, `10_T cells` = "Proliferative")
LJA <- RenameIdents(object = LJA, `14_NK cells` = "NK")

LJA <- subset(LJA, cells=rownames(LJA@reductions$umap@cell.embeddings[LJA@reductions$umap@cell.embeddings[,"umap_1"] >0 & LJA@reductions$umap@cell.embeddings[,"umap_2"] < 5,])) #> -3,]))
LJA_T_NK <- LJA

### mapping to atlases using LJA dataset 
LJA_T_NK_cellsbububian <- LJA_T_NK

#### 1. TICA

atlas_TICA_T_NK_cells <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells.rds"))
atlas_TICA_T_NK_cells <- subset(atlas_TICA_T_NK_cells,downsample=1000)
LJA_T_NK_cells <- LJA_T_NK_cellsbububian
anchors_map_LJA_T_NK_cells_to_TICA_T_NK_cells <- FindTransferAnchors(
  reference = atlas_TICA_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = LJA_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
LJA_T_NK_cells <- MapQuery(
  anchorset = anchors_map_LJA_T_NK_cells_to_TICA_T_NK_cells,
  query = LJA_T_NK_cells,
  reference = atlas_TICA_T_NK_cells,
  refdata = list(
    to_TICA_T_NK_cells = "cell_types_from_author",
    to_TICA_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_TICA_T_NK_cells = "integrated"#"RNA"
  ),
  new.reduction.name = "to_TICA_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_LJA_T_NK_cells_to_TICA_T_NK_cells")
rm("atlas_TICA_T_NK_cells")
LJA_T_NK_cells[["to_TICA_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(LJA_T_NK_cells,reduction = "ref.umap"),
                                                                    key = "to_TICA_T_NK_cells_",
                                                                    assay = DefaultAssay(LJA_T_NK_cells))
saveRDS(LJA_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_TICA_T_NK_cells_integrated_20250506.rds"))
rm("LJA_T_NK_cells")

#### pancancer T cells CD8 CD4
LJA_T_NK_cells <- LJA_T_NK_cellsbububian
atlas_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_filled.rds"))
anchors_map_LJA_T_NK_cells_to_pancancer_T_cells_CD8_CD4 <- FindTransferAnchors(
  reference = atlas_pancancer_T_cells_CD8_CD4,#atlas_pancancer_T_cells_CD8,
  query = LJA_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
LJA_T_NK_cells <- MapQuery(
  anchorset = anchors_map_LJA_T_NK_cells_to_pancancer_T_cells_CD8_CD4,
  query = LJA_T_NK_cells,
  reference = atlas_pancancer_T_cells_CD8_CD4,
  refdata = list(
    to_pancancer_T_cells_CD8_CD4 = "cell_types_from_author", #"cell_types_from_author_simplized"
    to_pancancer_T_cells_CD8_CD4_our_calculation = "cellcluster_ann",
    predicted_to_pancancer_T_cells_CD8_CD4 = "RNA"
  ),
  new.reduction.name = "to_pancancer_T_cells_CD8_CD4_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_LJA_T_NK_cells_to_pancancer_T_cells_CD8_CD4")
rm("atlas_pancancer_T_cells_CD8_CD4")
LJA_T_NK_cells[["to_pancancer_T_cells_CD8_CD4_umap"]] <- CreateDimReducObject(embeddings = Embeddings(LJA_T_NK_cells,reduction = "ref.umap"),
                                                                              key = "to_pancancer_T_cells_CD8_CD4_",
                                                                              assay = DefaultAssay(LJA_T_NK_cells))
saveRDS(LJA_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_pancancer_T_cells_CD8_CD4_20250506.rds"))
rm("LJA_T_NK_cells")

# 3. map to pancancer_blueprint_T_NK_cells
atlas_pancancer_blueprint_T_NK_cells <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells.rds"))
DefaultAssay(atlas_pancancer_blueprint_T_NK_cells) <- "integrated"
LJA_T_NK_cells <- LJA_T_NK_cellsbububian
anchors_map_LJA_T_NK_cells_to_pancancer_blueprint_T_NK_cells <- FindTransferAnchors(
  reference = atlas_pancancer_blueprint_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = LJA_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
LJA_T_NK_cells <- MapQuery(
  anchorset = anchors_map_LJA_T_NK_cells_to_pancancer_blueprint_T_NK_cells,
  query = LJA_T_NK_cells,
  reference = atlas_pancancer_blueprint_T_NK_cells,
  refdata = list(
    to_pancancer_blueprint_T_NK_cells = "cell_types_from_author",
    to_pancancer_blueprint_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_pancancer_blueprint_T_NK_cells = "integrated"#"RNA"
  ),
  new.reduction.name = "to_pancancer_blueprint_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_LJA_T_NK_cells_to_pancancer_blueprint_T_NK_cells")
rm("atlas_pancancer_blueprint_T_NK_cells")
LJA_T_NK_cells[["to_pancancer_blueprint_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(LJA_T_NK_cells,reduction = "ref.umap"),
                                                                                   key = "to_pancancer_blueprint_T_NK_cells_",
                                                                                   assay = DefaultAssay(LJA_T_NK_cells))
saveRDS(LJA_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated_20250506.rds"))
rm("LJA_T_NK_cells")

# 4. map to ABTC_T_NK_cells
LJA_T_NK_cells <- LJA_T_NK_cellsbububian
atlas_ABTC_T_NK_cells <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells.rds"))
anchors_map_LJA_T_NK_cells_to_ABTC_T_NK_cells <- FindTransferAnchors(
  reference = atlas_ABTC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = LJA_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
LJA_T_NK_cells <- MapQuery(
  anchorset = anchors_map_LJA_T_NK_cells_to_ABTC_T_NK_cells,
  query = LJA_T_NK_cells,
  reference = atlas_ABTC_T_NK_cells,
  refdata = list(
    to_ABTC_T_NK_cells = "cell_types_from_author",
    to_ABTC_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_ABTC_T_NK_cells = "RNA"
  ),
  new.reduction.name = "to_ABTC_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_LJA_T_NK_cells_to_ABTC_T_NK_cells")
rm("atlas_ABTC_T_NK_cells")
LJA_T_NK_cells[["to_ABTC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(LJA_T_NK_cells,reduction = "ref.umap"),
                                                                    key = "to_ABTC_T_NK_cells_",
                                                                    assay = DefaultAssay(LJA_T_NK_cells))
saveRDS(LJA_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_ABTC_T_NK_cells_20250506.rds"))
rm("LJA_T_NK_cells")

# 5. map to HCC_T_NK_cells
LJA_T_NK_cells <- LJA_T_NK_cellsbububian
atlas_HCC_T_NK_cells <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells.rds"))
DefaultAssay(atlas_HCC_T_NK_cells) <- "integrated"
anchors_map_LJA_T_NK_cells_to_HCC_T_NK_cells <- FindTransferAnchors(
  reference = atlas_HCC_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = LJA_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
LJA_T_NK_cells <- MapQuery(
  anchorset = anchors_map_LJA_T_NK_cells_to_HCC_T_NK_cells,
  query = LJA_T_NK_cells,
  reference = atlas_HCC_T_NK_cells,
  refdata = list(
    to_HCC_T_NK_cells = "cell_types_from_author",
    to_HCC_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_HCC_T_NK_cells = "RNA"
  ),
  new.reduction.name = "to_HCC_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_LJA_T_NK_cells_to_HCC_T_NK_cells")
rm("atlas_HCC_T_NK_cells")
LJA_T_NK_cells[["to_HCC_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(LJA_T_NK_cells,reduction = "ref.umap"),
                                                                   key = "to_HCC_T_NK_cells_",
                                                                   assay = DefaultAssay(LJA_T_NK_cells))
saveRDS(LJA_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_HCC_T_NK_cells_20250506.rds"))
rm("LJA_T_NK_cells")

# 6. map to STAD_T_NK_cells
LJA_T_NK_cells <- LJA_T_NK_cellsbububian
atlas_STAD_T_NK_cells <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells.rds"))
anchors_map_LJA_T_NK_cells_to_STAD_T_NK_cells <- FindTransferAnchors(
  reference = atlas_STAD_T_NK_cells,#atlas_pancancer_T_cells_CD8,
  query = LJA_T_NK_cells,#atlas_ABTC, ## query's counts or data cann't contain NA
  # features = VariableFeatures(atlas_pancancer_T_cells_CD8),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),#rownames(atlas_pancancer_T_cells_CD8@assays$RNA@scale.data),
  # normalization.method = "LogNormalize",
  reference.reduction = "pca",
  # k.filter = NA,
  dims = 1:30 #reference dim
)
LJA_T_NK_cells <- MapQuery(
  anchorset = anchors_map_LJA_T_NK_cells_to_STAD_T_NK_cells,
  query = LJA_T_NK_cells,
  reference = atlas_STAD_T_NK_cells,
  refdata = list(
    to_STAD_T_NK_cells = "cell_types_from_author",
    to_STAD_T_NK_cells_our_calculation = "cellcluster_ann",
    predicted_to_STAD_T_NK_cells = "RNA"
  ),
  new.reduction.name = "to_STAD_T_NK_cells_pca",
  reference.reduction = "pca",
  reduction.model = "umap"
)
rm("anchors_map_LJA_T_NK_cells_to_STAD_T_NK_cells")
rm("atlas_STAD_T_NK_cells")
LJA_T_NK_cells[["to_STAD_T_NK_cells_umap"]] <- CreateDimReducObject(embeddings = Embeddings(LJA_T_NK_cells,reduction = "ref.umap"),
                                                                    key = "to_STAD_T_NK_cells_",
                                                                    assay = DefaultAssay(LJA_T_NK_cells))
saveRDS(LJA_T_NK_cells,file=paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_STAD_T_NK_cells_20250506.rds"))
rm("LJA_T_NK_cells")


##########################################
#### 2024.02.28 plot umpas
###########################################
rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
changenames <- c("CD4 T","CD4 T","CD8 T","NK")
names(changenames) <- c("0_T memory cells","2_T memory cells","4_T cells","6_Gamma delta T cells")
new_cluster_names <- c("Recently activated CD4 T cells", "Naive-memory CD4 T cells", "Transitional memory CD4 T cells","Cytotoxic CD8 T cells",
                       "Effector memory CD8 T cells", "Th17 cells", "NK", "Terminally exhausted CD8 T cells",
                       "Naive T cells", "Regulatory T cells", "T helper cells", "Proliferative T cells",
                       "Pre-exhausted CD8 T cells",
                       
                       
                       "CD8.c07.Temra.CX3CR1", "CD8.c09.Tk.KIR2DL4", "CD8.c05.Tem.CXCR5" ,
                       "CD8.c15.ISG.IFIT1", "CD8.c01.Tn.MAL", "CD8.c17.Tm.NME1", "CD8.c08.Tk.TYROBP" ,
                       "CD8.c12.Tex.CXCL13", "CD8.c02.Tm.IL7R", "CD8.c10.Trm.ZNF683", "CD8.c03.Tm.RPS12" ,
                       "CD8.c14.Tex.TCF7", "CD8.c06.Tem.GZMK", "CD8.c11.Tex.PDCD1", "CD8.c16.MAIT.SLC4A10" ,
                       "CD8.c04.Tm.CD52", "CD8.c13.Tex.myl12a", "CD4.c20.Treg.TNFRSF9", "CD4.c18.Treg.RTKN2" ,
                       "CD4.c06.Tm.ANXA1", "CD4.c22.ISG.IFIT1", "CD4.c01.Tn.TCF7", "CD4.c12.Tem.GZMK" ,
                       "CD4.c11.Tm.GZMA", "CD4.c13.Temra.CX3CR1", "CD4.c16.Tfh.CXCR5", "CD4.c07.Tm.ANXA2" ,
                       "CD4.c17.TfhTh1.CXCL13", "CD4.c19.Treg.S1PR1", "CD4.c10.Tm.CAPG", "CD4.c24.Mix.NME2" ,
                       "CD4.c21.Treg.OAS1", "CD4.c14.Th17.SLC4A10", "CD4.c09.Tm.CCL5", "CD4.c02.Tn.PASK" ,
                       "CD4.c03.Tn.ADSL", "CD4.c05.Tm.TNF", "CD4.c08.Tm.CREM", "CD4.c04.Tn.il7r" ,
                       "CD4.c15.Th17.IL23R", "CD4.c23.Mix.NME1",
                       
                       
                       "CD8 memory T", "CD8 exhausted cytotoxic T",
                       "CD8 pre-effector T",   "CD4 memory/effector T",
                       "NK_XCL1",   "CD4 regulatory T",
                       "CD4 naive T",   "CD8 effector T",
                       "CD4 exhausted effector T", "Cytotoxic NK",
                       
                       "Naive CD4 T cell",     "Exhausted CD8 T cell",
                       "Memory CD8 T cell",    "Treg" ,               
                       "DNT",                  "Effector CD8 T cell" ,
                       "NCAM1nFCGR3Ap NK",     "Naive CD8 T cell" ,   
                       "NKT and NK cell",      "NCAM1pFCGR3Ap NK"   ,"NCAM1pFCGR3An NK",
                       
                       "13_CD8_effector memory T", "2_NK", "2_CD8_MAIT",  "20_CD4_naive T",
                       "9_CD4_regulatory T",  "8_NK",  "35_CD4_central memory T", "18_NK",
                       "11_CD8_inter-states T", "14_CD8_TR_memory T", "32_CD4_other T", "1_CD8_cytotoxic T" ,
                       "28_NK", "51_CD8_other T",
                       
                       "LNK1", "LNK3", "LNK4", "CD8 naive T" ,
                       "CD4 naive T cell",  "LNK2", "CD4 helper T",  "Regulatory T" ,
                       "CD8 effector T cell",  "Proliferative T")

names(new_cluster_names) <- c("Recently activated CD4 T cells",   "Naive-memory CD4 T cells",         "Transitional memory CD4 T cells",
                              "Cytotoxic CD8 T cells",            "Effector memory CD8 T cells",      "Th17 cells",
                              "NK",                               "Terminally exhausted CD8 T cells", "Naive T cells",
                              "Regulatory T cells",               "T helper cells",                   "Proliferative T cells",
                              "Pre-exhausted CD8 T cells",
                              
                              "CD8.c07.Temra.CX3CR1",   "CD8.c09.Tk.KIR2DL4",
                              "CD8.c05.Tem.CXCR5",      "CD8.c15.ISG.IFIT1",      "CD8.c01.Tn.MAL",
                              "CD8.c17.Tm.NME1",        "CD8.c08.Tk.TYROBP",      "CD8.c12.Tex.CXCL13",
                              "CD8.c02.Tm.IL7R",        "CD8.c10.Trm.ZNF683",     "CD8.c03.Tm.RPS12",
                              "CD8.c14.Tex.TCF7",       "CD8.c06.Tem.GZMK",       "CD8.c11.Tex.PDCD1",
                              "CD8.c16.MAIT.SLC4A10",   "CD8.c04.Tm.CD52",        "CD8.c13.Tex.myl12a",
                              "CD4.c20.Treg.TNFRSF9",   "CD4.c18.Treg.RTKN2",     "CD4.c06.Tm.ANXA1",
                              "CD4.c22.ISG.IFIT1",      "CD4.c01.Tn.TCF7",        "CD4.c12.Tem.GZMK",
                              "CD4.c11.Tm.GZMA",        "CD4.c13.Temra.CX3CR1",   "CD4.c16.Tfh.CXCR5",
                              "CD4.c07.Tm.ANXA2",       "CD4.c17.TfhTh1.CXCL13",  "CD4.c19.Treg.S1PR1",
                              "CD4.c10.Tm.CAPG",        "CD4.c24.Mix.NME2",       "CD4.c21.Treg.OAS1",
                              "CD4.c14.Th17.SLC4A10",   "CD4.c09.Tm.CCL5",        "CD4.c02.Tn.PASK",
                              "CD4.c03.Tn.ADSL",        "CD4.c05.Tm.TNF",         "CD4.c08.Tm.CREM",
                              "CD4.c04.Tn.il7r",        "CD4.c15.Th17.IL23R",     "CD4.c23.Mix.NME1",
                              
                              
                              "C3_CD8_ZNF683",     "C1_CD8_HAVCR2",     "C2_CD8_GZMK",
                              "C6_CD4_GZMA",       "C10_NK_XCL1",       "C8_CD4_FOXP3",
                              "C5_CD4_CCR7",       "C4_CD8_CX3CR1",     "C7_CD4_CXCL13",
                              "C9_NK_FGFBP2",
                              
                              "Naive CD4 T cell",                 "Exhausted CD8 T cell",
                              "Memory CD8 T cell",                "Treg",                             "DNT",
                              "Effector CD8 T cell",              "NCAM1-FCGR3A+ NK",                 "Naive CD8 T cell",
                              "NKT and NK cell",                  "NCAM1+FCGR3A+ NK",                 "NCAM1+FCGR3A- NK",
                              
                              
                              "13_T/NK",                           "25_T/NK",                           "2_T/NK",
                              "20_T/NK",                           "9_T/NK",                            "8_T/NK",
                              "35_T/NK",                           "18_T/NK",                           "11_T/NK",
                              "14_T/NK",                           "32_T/NK",                           "1_T/NK",
                              "28_T/NK",                           "51_T/NK",                           "LNK1",
                              
                              
                              "LNK3",                             "LNK4",                             "LT2",
                              "LT3",                              "LNK2",                             "LT4",
                              "LT5",                              "LT1",                              "LPT")

LJA_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_TICA_T_NK_cells_integrated_20250506.rds"))
LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)])
DimPlot(LJA_T_NK_cells_TICA,pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(LJA_T_NK_cells_TICA,group.by = "predicted.to_TICA_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

LJA_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_pancancer_T_cells_CD8_CD4_20250506.rds"))
LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4 <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)])
# DimPlot(LJA_T_NK_cells_panT,group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(LJA_T_NK_cells_panT,group.by = "predicted.to_pancancer_T_cells_CD8_CD4",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

LJA_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated_20250506.rds"))
LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)])
# DimPlot(LJA_T_NK_cells_panblueprint,group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(LJA_T_NK_cells_panblueprint,group.by = "predicted.to_pancancer_blueprint_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

LJA_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_ABTC_T_NK_cells_20250506.rds"))
LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)])
# DimPlot(LJA_T_NK_cells_ABTC,group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(LJA_T_NK_cells_ABTC,group.by = "predicted.to_ABTC_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

LJA_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_HCC_T_NK_cells_20250506.rds"))
LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)])
LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells[LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells == "2_NK"] <- "25_NK"
# DimPlot(LJA_T_NK_cells_HCC,group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(LJA_T_NK_cells_HCC,group.by = "predicted.to_HCC_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )

LJA_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_STAD_T_NK_cells_20250506.rds"))
LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)])
# DimPlot(LJA_T_NK_cells_STAD,group.by = "cell_types_from_author",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )
DimPlot(LJA_T_NK_cells_STAD,group.by = "predicted.to_STAD_T_NK_cells",pt.size=2) + theme(legend.text=element_text(size=36)) + guides(color = guide_legend(override.aes = list(size=7), ncol=1) )


##############################
#### 2024.02.28 barplot for LJA dataset
#### 2025.05.06 sankey plots
##############################

rm(list=ls())
library(Seurat)
library(DESeq2)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
library(scales)
"%ni%" <- Negate("%in%")

workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
new_cluster_names <- c("Recently activated CD4 T cells", "Naive-memory CD4 T cells", "Transitional memory CD4 T cells","Cytotoxic CD8 T cells",
                       "Effector memory CD8 T cells", "Th17 cells", "NK", "Terminally exhausted CD8 T cells",
                       "Naive T cells", "Regulatory T cells", "T helper cells", "Proliferative T cells",
                       "Pre-exhausted CD8 T cells",
                       
                       
                       "CD8.c07.Temra.CX3CR1", "CD8.c09.Tk.KIR2DL4", "CD8.c05.Tem.CXCR5" ,
                       "CD8.c15.ISG.IFIT1", "CD8.c01.Tn.MAL", "CD8.c17.Tm.NME1", "CD8.c08.Tk.TYROBP" ,
                       "CD8.c12.Tex.CXCL13", "CD8.c02.Tm.IL7R", "CD8.c10.Trm.ZNF683", "CD8.c03.Tm.RPS12" ,
                       "CD8.c14.Tex.TCF7", "CD8.c06.Tem.GZMK", "CD8.c11.Tex.PDCD1", "CD8.c16.MAIT.SLC4A10" ,
                       "CD8.c04.Tm.CD52", "CD8.c13.Tex.myl12a", "CD4.c20.Treg.TNFRSF9", "CD4.c18.Treg.RTKN2" ,
                       "CD4.c06.Tm.ANXA1", "CD4.c22.ISG.IFIT1", "CD4.c01.Tn.TCF7", "CD4.c12.Tem.GZMK" ,
                       "CD4.c11.Tm.GZMA", "CD4.c13.Temra.CX3CR1", "CD4.c16.Tfh.CXCR5", "CD4.c07.Tm.ANXA2" ,
                       "CD4.c17.TfhTh1.CXCL13", "CD4.c19.Treg.S1PR1", "CD4.c10.Tm.CAPG", "CD4.c24.Mix.NME2" ,
                       "CD4.c21.Treg.OAS1", "CD4.c14.Th17.SLC4A10", "CD4.c09.Tm.CCL5", "CD4.c02.Tn.PASK" ,
                       "CD4.c03.Tn.ADSL", "CD4.c05.Tm.TNF", "CD4.c08.Tm.CREM", "CD4.c04.Tn.il7r" ,
                       "CD4.c15.Th17.IL23R", "CD4.c23.Mix.NME1",
                       
                       
                       "CD8 memory T", "CD8 exhausted cytotoxic T",
                       "CD8 pre-effector T",   "CD4 memory/effector T",
                       "NK_XCL1",   "CD4 regulatory T",
                       "CD4 naive T",   "CD8 effector T",
                       "CD4 exhausted effector T", "Cytotoxic NK",
                       
                       "Naive CD4 T cell",     "Exhausted CD8 T cell",
                       "Memory CD8 T cell",    "Treg" ,               
                       "DNT",                  "Effector CD8 T cell" ,
                       "NCAM1nFCGR3Ap NK",     "Naive CD8 T cell" ,   
                       "NKT and NK cell",      "NCAM1pFCGR3Ap NK"   ,"NCAM1pFCGR3An NK",
                       
                       "13_CD8_effector memory T", "2_NK", "2_CD8_MAIT",  "20_CD4_naive T",
                       "9_CD4_regulatory T",  "8_NK",  "35_CD4_central memory T", "18_NK",
                       "11_CD8_inter-states T", "14_CD8_TR_memory T", "32_CD4_other T", "1_CD8_cytotoxic T" ,
                       "28_NK", "51_CD8_other T",
                       
                       "LNK1", "LNK3", "LNK4", "CD8 naive T" ,
                       "CD4 naive T cell",  "LNK2", "CD4 helper T",  "Regulatory T" ,
                       "CD8 effector T cell",  "Proliferative T")

names(new_cluster_names) <- c("Recently activated CD4 T cells",   "Naive-memory CD4 T cells",         "Transitional memory CD4 T cells",
                              "Cytotoxic CD8 T cells",            "Effector memory CD8 T cells",      "Th17 cells",
                              "NK",                               "Terminally exhausted CD8 T cells", "Naive T cells",
                              "Regulatory T cells",               "T helper cells",                   "Proliferative T cells",
                              "Pre-exhausted CD8 T cells",
                              
                              "CD8.c07.Temra.CX3CR1",   "CD8.c09.Tk.KIR2DL4",
                              "CD8.c05.Tem.CXCR5",      "CD8.c15.ISG.IFIT1",      "CD8.c01.Tn.MAL",
                              "CD8.c17.Tm.NME1",        "CD8.c08.Tk.TYROBP",      "CD8.c12.Tex.CXCL13",
                              "CD8.c02.Tm.IL7R",        "CD8.c10.Trm.ZNF683",     "CD8.c03.Tm.RPS12",
                              "CD8.c14.Tex.TCF7",       "CD8.c06.Tem.GZMK",       "CD8.c11.Tex.PDCD1",
                              "CD8.c16.MAIT.SLC4A10",   "CD8.c04.Tm.CD52",        "CD8.c13.Tex.myl12a",
                              "CD4.c20.Treg.TNFRSF9",   "CD4.c18.Treg.RTKN2",     "CD4.c06.Tm.ANXA1",
                              "CD4.c22.ISG.IFIT1",      "CD4.c01.Tn.TCF7",        "CD4.c12.Tem.GZMK",
                              "CD4.c11.Tm.GZMA",        "CD4.c13.Temra.CX3CR1",   "CD4.c16.Tfh.CXCR5",
                              "CD4.c07.Tm.ANXA2",       "CD4.c17.TfhTh1.CXCL13",  "CD4.c19.Treg.S1PR1",
                              "CD4.c10.Tm.CAPG",        "CD4.c24.Mix.NME2",       "CD4.c21.Treg.OAS1",
                              "CD4.c14.Th17.SLC4A10",   "CD4.c09.Tm.CCL5",        "CD4.c02.Tn.PASK",
                              "CD4.c03.Tn.ADSL",        "CD4.c05.Tm.TNF",         "CD4.c08.Tm.CREM",
                              "CD4.c04.Tn.il7r",        "CD4.c15.Th17.IL23R",     "CD4.c23.Mix.NME1",
                              
                              
                              "C3_CD8_ZNF683",     "C1_CD8_HAVCR2",     "C2_CD8_GZMK",
                              "C6_CD4_GZMA",       "C10_NK_XCL1",       "C8_CD4_FOXP3",
                              "C5_CD4_CCR7",       "C4_CD8_CX3CR1",     "C7_CD4_CXCL13",
                              "C9_NK_FGFBP2",
                              
                              "Naive CD4 T cell",                 "Exhausted CD8 T cell",
                              "Memory CD8 T cell",                "Treg",                             "DNT",
                              "Effector CD8 T cell",              "NCAM1-FCGR3A+ NK",                 "Naive CD8 T cell",
                              "NKT and NK cell",                  "NCAM1+FCGR3A+ NK",                 "NCAM1+FCGR3A- NK",
                              
                              
                              "13_T/NK",                           "25_T/NK",                           "2_T/NK",
                              "20_T/NK",                           "9_T/NK",                            "8_T/NK",
                              "35_T/NK",                           "18_T/NK",                           "11_T/NK",
                              "14_T/NK",                           "32_T/NK",                           "1_T/NK",
                              "28_T/NK",                           "51_T/NK",                           "LNK1",
                              
                              
                              "LNK3",                             "LNK4",                             "LT2",
                              "LT3",                              "LNK2",                             "LT4",
                              "LT5",                              "LT1",                              "LPT")

LJA_T_NK_cells_TICA <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_TICA_T_NK_cells_integrated_20250506.rds"))
LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)])
cell_summary_TICA <- data.frame(Idents(LJA_T_NK_cells_TICA),paste0("TICA: ",LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells))
# rm("LJA_T_NK_cells_TICA")
LJA_T_NK_cells_panT <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_pancancer_T_cells_CD8_CD4_20250506.rds"))
LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4 <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)])
cell_summary_panT <- data.frame(Idents(LJA_T_NK_cells_TICA),paste0("pan- T: ",LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4))
# rm("LJA_T_NK_cells_panT")
LJA_T_NK_cells_panblueprint <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_pancancer_blueprint_T_NK_cells_integrated_20250506.rds"))
LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)])
cell_summary_panblueprint <- data.frame(Idents(LJA_T_NK_cells_TICA),paste0("pan- blueprint: ",LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))
# rm("LJA_T_NK_cells_panblueprint")
LJA_T_NK_cells_ABTC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_ABTC_T_NK_cells_20250506.rds"))
LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)])
cell_summary_ABTC <- data.frame(Idents(LJA_T_NK_cells_TICA),paste0("ABTC: ",LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells))
# rm("LJA_T_NK_cells_ABTC")
LJA_T_NK_cells_HCC <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_HCC_T_NK_cells_20250506.rds"))
LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)])
LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells[LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells == "2_NK"] <- "25_NK"
cell_summary_HCC <- data.frame(Idents(LJA_T_NK_cells_TICA),paste0("HCC: ",LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells))
# rm("LJA_T_NK_cells_HCC")
LJA_T_NK_cells_STAD <- readRDS(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_STAD_T_NK_cells_20250506.rds"))
LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells <- as.character(new_cluster_names[as.character(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)])
cell_summary_STAD <- data.frame(Idents(LJA_T_NK_cells_TICA),paste0("STAD: ",LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells))
# rm("LJA_T_NK_cells_STAD")

cell_summary <- cbind(cell_summary_TICA,to_panT=cell_summary_panT[,2],
                      to_panblueprint=cell_summary_panblueprint[,2],to_ABTC=cell_summary_ABTC[,2],to_HCC=cell_summary_HCC[,2],
                      to_STAD=cell_summary_STAD[,2])
colnames(cell_summary)[1:2] <-c("cell_types_from_LJA","to_TICA")
# cell_summary[1:3,]
# cell_summary[,"cell_types_from_LJA"] <- as.character(cell_summary[,"cell_types_from_LJA"])
# cell_summary[,"hpca_fine"] <- as.character(cell_summary[,"hpca_fine"])

unique(cell_summary$cell_types_from_LJA)
# [1] CD8_c1_Tex    CD4_c1_Treg   Tgd_c3        CD8_c2_Teff   CD4_c2_Tn     NKT_c0        Proliferative NK       




library(networkD3)
library(ggalluvial)
colnames_atlas <- colnames(cell_summary)
my_colors <- c(
  "Proliferative"         = "#E41A1C",
  "CD4_c1_Treg"     = "#377EB8",
  "CD4_c2_Tn"   = "#4DAF4A",
  "CD8_c1_Tex"      = "#FF7F00",
  "CD8_c2_Teff"      = "#984EA3",
  "Tgd_c3"        = "#00CED1",
  "NKT_c0" = "#FFC107",
  "NK"   = "#6A3D9A"
)
### For TICA
links <- dplyr::count(cell_summary,cell_types_from_LJA,to_TICA)
sankeyCols <- c("source", "target", "value")
colnames(links) <- sankeyCols
links <- links[links$value >10,]
links$source <- factor(links$source,levels = c("Proliferative",
                                               "CD4_c1_Treg","CD4_c2_Tn",
                                               "CD8_c1_Tex","CD8_c2_Teff","Tgd_c3","NKT_c0",
                                               "NK"))
links$target <- factor(links$target,levels = c("TICA: Proliferative T cells",
                                               "TICA: Regulatory T cells","TICA: T helper cells","TICA: Naive-memory CD4 T cells","TICA: Th17 cells","TICA: Transitional memory CD4 T cells","TICA: Naive T cells",
                                               "TICA: Cytotoxic CD8 T cells","TICA: Recently activated CD4 T cells","TICA: Effector memory CD8 T cells","TICA: Pre-exhausted CD8 T cells","TICA: Terminally exhausted CD8 T cells",
                                               "TICA: NK" ))

# pdf(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.pdf"),width=10)
# png(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.png"),width=1800)
p6 <- ggplot(data = links,
       aes(axis1 = source, axis2 = target, y = value)) +
  geom_alluvium(aes(fill = source),width = 1/10) +
  geom_stratum(alpha = .1,width=1/10,linewidth=1) +#control bar width
  # geom_bar() + # remove empty lines
  # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
  scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
  scale_fill_manual(values=my_colors) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
    aes(label = as.character(target)),
    stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
    aes(label = as.character(source)),
    stat = "stratum", size =6, direction = "y", nudge_x = -0.3) +
  theme(legend.position = "none")
# dev.off()

### For panT
links <- dplyr::count(cell_summary,cell_types_from_LJA,to_panT)
sankeyCols <- c("source", "target", "value")
colnames(links) <- sankeyCols
# links <- links[links$value >10,]
links$source <- factor(links$source,levels = c("Proliferative",
                                               "CD4_c1_Treg","CD4_c2_Tn",
                                               "CD8_c1_Tex","CD8_c2_Teff","Tgd_c3","NKT_c0",
                                               "NK"))
links$target <- factor(links$target,levels = c("pan- T: CD4.c03.Tn.ADSL","pan- T: CD4.c10.Tm.CAPG","pan- T: CD4.c24.Mix.NME2","pan- T: CD8.c07.Temra.CX3CR1","pan- T: CD8.c10.Trm.ZNF683",  
                                               "pan- T: CD4.c20.Treg.TNFRSF9","pan- T: CD8.c15.ISG.IFIT1","pan- T: CD4.c22.ISG.IFIT1"   ))

# pdf(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.pdf"),width=10)
# png(paste0(workdir,"2.5_mapping_known_dataset_Bi_T_NK_cells_to_TICA_T_NK_cells_flow_20250506.png"),width=1800)
p0_nouse <- ggplot(data = links,
             aes(axis1 = source, axis2 = target, y = value)) +
  geom_alluvium(aes(fill = source),width = 1/10) +
  geom_stratum(alpha = .1,width=1/10,linewidth=1) +#control bar width
  # geom_bar() + # remove empty lines
  # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
  scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
  scale_fill_manual(values=my_colors) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
    aes(label = as.character(target)),
    stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
    aes(label = as.character(source)),
    stat = "stratum", size =6, direction = "y", nudge_x = -0.3) +
  theme(legend.position = "none")
# dev.off()


## For blueprint
links <- dplyr::count(cell_summary,cell_types_from_LJA,to_panblueprint)
sankeyCols <- c("source", "target", "value")
colnames(links) <- sankeyCols
links <- links[links$value >10,]
links$source <- factor(links$source,levels = c("Proliferative",
                                               "CD4_c1_Treg","CD4_c2_Tn",
                                               "CD8_c1_Tex","CD8_c2_Teff","Tgd_c3","NKT_c0",
                                               "NK"))
links$target <- factor(links$target,levels = c("pan- blueprint: CD4 regulatory T","pan- blueprint: CD4 memory/effector T","pan- blueprint: CD4 exhausted effector T","pan- blueprint: CD4 naive T",
                                               "pan- blueprint: CD8 pre-effector T","pan- blueprint: CD8 exhausted cytotoxic T","pan- blueprint: CD8 memory T","pan- blueprint: CD8 effector T",
                                               "pan- blueprint: NK_XCL1","pan- blueprint: Cytotoxic NK"))
# pdf(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_pancancer_blueprint_T_NK_cells_flow_20250506.pdf"),width=10)
p7 <- ggplot(data = links,
       aes(axis1 = source, axis2 = target, y = value)) +
  geom_alluvium(aes(fill = source),width = 1/10) +
  geom_stratum(alpha = .1,width=1/10,linewidth=1.5) + #control bar width
  # geom_bar() +
  # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
  scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
  scale_fill_manual(values=my_colors) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
    aes(label = as.character(target)),
    stat = "stratum", size = 6, direction = "y", nudge_x = 0.7) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
    aes(label = as.character(source)),
    stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
  theme(legend.position = "none")
# dev.off()


## For ABTC
links <- dplyr::count(cell_summary,cell_types_from_LJA,to_ABTC)
sankeyCols <- c("source", "target", "value")
colnames(links) <- sankeyCols
links <- links[links$value >10,]
links$source <- factor(links$source,levels = c("Proliferative",
                                               "CD4_c1_Treg","CD4_c2_Tn",
                                               "CD8_c1_Tex","CD8_c2_Teff","Tgd_c3","NKT_c0",
                                               "NK"))
links$target <- factor(links$target,levels = c("ABTC: Treg","ABTC: Naive CD4 T cell","ABTC: DNT",
                                               "ABTC: Exhausted CD8 T cell","ABTC: Effector CD8 T cell","ABTC: Memory CD8 T cell","ABTC: Naive CD8 T cell",
                                               "ABTC: NCAM1pFCGR3An NK","ABTC: NCAM1pFCGR3Ap NK","ABTC: NCAM1nFCGR3Ap NK" ))
# pdf(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_ABTC_T_NK_cells_flow_20250506.pdf"),width=10)
p8 <- ggplot(data = links,
       aes(axis1 = source, axis2 = target, y = value)) +
  geom_alluvium(aes(fill = source),width = 1/10) +
  geom_stratum(alpha = .1,width=1/10,linewidth=1.5) + #control bar width
  # geom_bar() +
  # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
  scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
  scale_fill_manual(values=my_colors) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
    aes(label = as.character(target)),
    stat = "stratum", size = 6, direction = "y", nudge_x = 0.6) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
    aes(label = as.character(source)),
    stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
  theme(legend.position = "none")
# dev.off()

## For HCC
links <- dplyr::count(cell_summary,cell_types_from_LJA,to_HCC)
sankeyCols <- c("source", "target", "value")
colnames(links) <- sankeyCols
links <- links[links$value >10,]
links$source <- factor(links$source,levels = c("Proliferative",
                                               "CD4_c1_Treg","CD4_c2_Tn",
                                               "CD8_c1_Tex","CD8_c2_Teff","Tgd_c3","NKT_c0",
                                               "NK"))
links$target <- factor(links$target,levels = c("HCC: 32_CD4_other T","HCC: 51_CD8_other T",
                                               "HCC: 9_CD4_regulatory T","HCC: 20_CD4_naive T","HCC: 35_CD4_central memory T",
                                               "HCC: 28_NK","HCC: 13_CD8_effector memory T","HCC: 11_CD8_inter-states T","HCC: 1_CD8_cytotoxic T","HCC: 2_CD8_MAIT","HCC: 14_CD8_TR_memory T",
                                               "HCC: 8_NK","HCC: 25_NK", "HCC: 18_NK"))
# pdf(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_HCC_T_NK_cells_flow_20250506.pdf"),width=10,height=4)
p9 <- ggplot(data = links,
       aes(axis1 = source, axis2 = target, y = value)) +
  geom_alluvium(aes(fill = source),width = 1/10) +
  geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
  # geom_bar() +
  # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
  scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
  scale_fill_manual(values=my_colors) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
    aes(label = as.character(target)),
    stat = "stratum", size = 6, direction = "y", nudge_x = 0.6) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
    aes(label = as.character(source)),
    stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
  theme(legend.position = "none")
# dev.off()

## For STAD
links <- dplyr::count(cell_summary,cell_types_from_LJA,to_STAD)
sankeyCols <- c("source", "target", "value")
colnames(links) <- sankeyCols
links <- links[links$value >10,]
links$source <- factor(links$source,levels = c("Proliferative",
                                               "CD4_c1_Treg","CD4_c2_Tn",
                                               "CD8_c1_Tex","CD8_c2_Teff","Tgd_c3","NKT_c0",
                                               "NK"))
links$target <- factor(links$target,levels = c("STAD: Proliferative T",
                                               "STAD: Regulatory T","STAD: CD4 naive T cell","STAD: CD8 effector T cell","STAD: CD4 helper T",
                                               "STAD: LNK2","STAD: LNK3","STAD: CD8 naive T",
                                               "STAD: LNK1","STAD: LNK4"))
# pdf(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_STAD_T_NK_cells_flow_20250506.pdf"),width=15,height = 6)
p10 <- ggplot(data = links,
       aes(axis1 = source, axis2 = target, y = value)) +
  geom_alluvium(aes(fill = source),width = 1/10) +
  geom_stratum(alpha = .1,width=1/10,linewidth=1) + #control bar width
  # geom_bar() +
  # geom_text(stat = "stratum",size =6, aes(label = after_stat(stratum)),min.y=10) +
  scale_x_discrete(limits = c("source", "target"),expand = c(0.75, 0.75)) + theme_void() +
  scale_fill_manual(values=my_colors) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 3, as.character(target), NA)),
    aes(label = as.character(target)),
    stat = "stratum", size = 6, direction = "y", nudge_x = 0.5) +
  ggrepel::geom_text_repel(
    # aes(label = ifelse(as.numeric(source) == 1, after_stat(stratum), after_stat(stratum))),
    aes(label = as.character(source)),
    stat = "stratum", size = 6, direction = "y", nudge_x = -0.3) +
  theme(legend.position = "none")
# dev.off()

pdf(paste0(workdir,"2.5_mapping_known_dataset_LJA_T_NK_cells_to_allatlases_T_NK_cells_flow_20250506.pdf"),height=25,width=15)
print(CombinePlots(list(p6,p7,p8,p9,p10),ncol=1))
dev.off()

#### calculate the ARI
adjustedRandIndex(cell_summary$cell_types_from_LJA,cell_summary$to_TICA)
# [1] 0.3168938
adjustedRandIndex(cell_summary$cell_types_from_LJA,cell_summary$to_panblueprint)
# [1] 0.2950238
adjustedRandIndex(cell_summary$cell_types_from_LJA,cell_summary$to_ABTC)
# [1] 0.3048876c
adjustedRandIndex(cell_summary$cell_types_from_LJA,cell_summary$to_HCC)
# [1] 0.2375804
adjustedRandIndex(cell_summary$cell_types_from_LJA,cell_summary$to_STAD)
# [1] 0.2454156


# index <- c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD")
# index2 <- c("TICA: ","pan- T: ","pan- blueprint: ","ABTC: ","HCC: ","STAD: ")
# # for T_cell:CD4+_effector_memory
# output0 <- NULL; output1 <- NULL;output2 <- NULL;output3<-NULL
# pp <- NULL
# # for (i in 1:length(index)) {
#   
#     i=1
#   
#   #### for p1
#   
#   cell_summary_subset1 <- subset(cell_summary,cell_types_from_LJA %in% c("T.CD4")) #
#   cell_summary_subset1$hpca_fine <- "CD4+ T"
#   # cell_summary_subset2 <- subset(cell_summary,cell_types_from_LJA %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
#   # cell_summary_subset3 <- subset(cell_summary,cell_types_from_LJA %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
#   cell_summary_subset <- cell_summary_subset1 #rbind(cell_summary_subset1,cell_summary_subset2,cell_summary_subset3)
#   output_pre1 <- as.data.frame(table(cell_summary_subset1[,c(index[i],"hpca_fine")]))
#   output_pre1 <- data.frame(group=index[i],output_pre1)
#   output_pre1[,2] <- as.character(output_pre1[,2])
#   output_pre1$Freq <- as.numeric(output_pre1$Freq)
#   # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
#   colnames(output_pre1)[2] <- "fill"
#   # output_pre1$fill <- as.character(new_cluster_names[output_pre1$fill])
#   output_pre1 <- data.frame(output_pre1,populations=output_pre1$fill)
#   output_pre1$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre1$fill),": ",2)))[,2])
#   output_pre1$Freq <- as.numeric(output_pre1$Freq)
#   output_pre1 <- output_pre1[order(output_pre1$Freq),]
#   output_pre1 <- output_pre1 %>%
#     group_by(group) %>%
#     mutate(label_y = cumsum(Freq)-0.5*Freq)
#   output_pre1 <- output_pre1[order(-output_pre1$Freq),]
#   output1 <- rbind(output1,output_pre1)
#   
#   output1 <- as.data.frame(output1)
#     # output1 <- rbind(output1,c("RF","CD4+ T","CD4+ T",sum(output_pre1$Freq),"RF: CD4+ T",sum(output_pre1$Freq)/2))#;output1 <- output
#   # output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2));output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output1,output2,output3)
#   output1$fill <- factor(output1$fill,levels = c(output1$fill))
#   output1$populations <- factor(output1$populations,levels = c(output1$populations))
#   # output1$group <- factor(output1$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output1$group <- factor(output1$group,levels = c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output1$Freq <- as.numeric(output1$Freq)
#   output1$label_y <- as.numeric(output1$label_y)
#   output1[output1$Freq < 3,"fill"] <- "" #10
#   if (i ==1) {
#     cols <- hue_pal()(length(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#     names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   }
#   if (i ==2) {
#     cols <- hue_pal()(length(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#     names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   }
#   if (i ==3) {
#     cols <- hue_pal()(length(unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
#     # names(cols) <- sort(as.character(new_cluster_names[paste0("pan- blueprint: ",unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))]))
#     names(cols) <- sort(as.character(paste0("pan- blueprint: ",unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells))))
#   }
#   if (i ==4) {
#     cols <- hue_pal()(length(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
#     # names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
#     names(cols) <- sort(as.character(paste0("ABTC: ",unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells))))
#   }
#   if (i ==5) {
#     cols <- hue_pal()(length(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
#     # names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
#     names(cols) <- as.character(paste0("HCC: ",sort(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells))))
#   }
#   if (i ==6) {
#     cols <- hue_pal()(length(unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
#     # names(cols) <- sort(as.character(new_cluster_names[paste0("STAD: ",unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells))]))
#     names(cols) <- sort(as.character(paste0("STAD: ",unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells))))
#   }
#   
#   
#   
#   # cols <- c("RF: Proliferative cells" = "purple", "RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
#   p1 <- ggplot(output1,aes(x=group,y=Freq,fill=populations)) + 
#     geom_bar(stat="identity",colour="black") + 
#     scale_fill_manual(values=cols) +
#     # geom_text(aes(y=label_y,label=fill),size =10) +
#     theme_classic()
#   
#   #### for p0
#   
#   # cell_summary_subset3 <- subset(cell_summary,cell_types_from_LJA %in% c("0_T memory cells","2_T memory cells")) #
#   # cell_summary_subset3$hpca_fine <- "CD4+ T"
#   # cell_summary_subset2 <- subset(cell_summary,cell_types_from_LJA %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
#   cell_summary_subset0 <- subset(cell_summary,cell_types_from_LJA %in% c("T.cell"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   cell_summary_subset0$hpca_fine <- "Proliferative cells"#"CD4+ T"
#   # cell_summary_subset <- cell_summary_subset3 #rbind(cell_summary_subset3,cell_summary_subset2,cell_summary_subset3)
#   output_pre0 <- as.data.frame(table(cell_summary_subset0[,c(index[i],"hpca_fine")]))
#   output_pre0 <- data.frame(group=index[i],output_pre0)
#   output_pre0[,2] <- as.character(output_pre0[,2])
#   output_pre0$Freq <- as.numeric(output_pre0$Freq)
#   # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
#   colnames(output_pre0)[2] <- "fill"
#   # output_pre3$fill <- as.character(new_cluster_names[output_pre3$fill])
#   output_pre0 <- data.frame(output_pre0,populations=output_pre0$fill)
#   output_pre0$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre0$fill),": ",2)))[,2])
#   output_pre0$Freq <- as.numeric(output_pre0$Freq)
#   output_pre0 <- output_pre0[order(output_pre0$Freq),]
#   output_pre0 <- output_pre0 %>%
#     group_by(group) %>%
#     mutate(label_y = cumsum(Freq)-0.5*Freq)
#   output_pre0 <- output_pre0[order(-output_pre0$Freq),]
#   output0 <- rbind(output0,output_pre0)
#   
#   output0 <- as.data.frame(output0)
#   # output3 <- rbind(output3,c("RF","CD4+ T","CD4+ T",sum(output_pre3$Freq),"RF: CD4+ T",sum(output_pre3$Freq)/2))#;output3 <- output
#   # output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output0 <- rbind(output0,c("RF","Proliferative cell","Proliferative cells",sum(output_pre0$Freq),"RF: Proliferative cells",sum(output_pre0$Freq)/2))#;output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output3,output2,output3)
#   output0$fill <- factor(output0$fill,levels = c(output0$fill))
#   output0$populations <- factor(output0$populations,levels = c(output0$populations))
#   # output0$group <- factor(output0$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output0$group <- factor(output0$group,levels = c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output0$Freq <- as.numeric(output0$Freq)
#   output0$label_y <- as.numeric(output0$label_y)
#   output0[output0$Freq < 3,"fill"] <- "" #30
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("pan- blueprint: ",sort(unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("STAD: ",sort(unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))])
#   
#   
#   # cols <- c("RF: Proliferative cells" = "purple", "RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
#   p0 <- ggplot(output0,aes(x=group,y=Freq,fill=populations)) + 
#     geom_bar(stat="identity",colour="black") + 
#     scale_fill_manual(values=cols) +
#     # geom_text(aes(y=label_y,label=fill),size =10) +
#     theme_classic()
#   
#   ##### for p2
#   
#   # cell_summary_subset2 <- subset(cell_summary,cell_types_from_LJA %in% c("0_T memory cells","2_T memory cells")) #
#   # cell_summary_subset2$hpca_fine <- "CD4+ T"
#   cell_summary_subset2 <- subset(cell_summary,cell_types_from_LJA %in% c("T.CD8"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
#   # cell_summary_subset3 <- subset(cell_summary,cell_types_from_LJA %in% c("6_Gamma delta T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
#   # cell_summary_subset <- cell_summary_subset2 #rbind(cell_summary_subset2,cell_summary_subset2,cell_summary_subset3)
#   output_pre2 <- as.data.frame(table(cell_summary_subset2[,c(index[i],"hpca_fine")]))
#   output_pre2 <- data.frame(group=index[i],output_pre2)
#   output_pre2[,2] <- as.character(output_pre2[,2])
#   output_pre2$Freq <- as.numeric(output_pre2$Freq)
#   # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
#   colnames(output_pre2)[2] <- "fill"
#   # output_pre2$fill <- as.character(new_cluster_names[output_pre2$fill])
#   output_pre2 <- data.frame(output_pre2,populations=output_pre2$fill)
#   output_pre2$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre2$fill),": ",2)))[,2])
#   output_pre2$Freq <- as.numeric(output_pre2$Freq)
#   output_pre2 <- output_pre2[order(output_pre2$Freq),]
#   output_pre2 <- output_pre2 %>%
#     group_by(group) %>%
#     mutate(label_y = cumsum(Freq)-0.5*Freq)
#   output_pre2 <- output_pre2[order(-output_pre2$Freq),]
#   output2 <- rbind(output2,output_pre2)
#   
#   output2 <- as.data.frame(output2)
#   # output2 <- rbind(output2,c("RF","CD4+ T","CD4+ T",sum(output_pre2$Freq),"RF: CD4+ T",sum(output_pre2$Freq)/2))#;output2 <- output
#   # output2 <- rbind(output2,c("RF","CD8+ T","CD8+ T",sum(output_pre2$Freq),"RF: CD8+ T",sum(output_pre2$Freq)/2))
#   # ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output,c("RF","NK cell","NK",sum(output_pre$Freq),"RF: NK",sum(output_pre$Freq)/2));output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output2,output2,output3)
#   output2$fill <- factor(output2$fill,levels = c(output2$fill))
#   output2$populations <- factor(output2$populations,levels = c(output2$populations))
#   # output2$group <- factor(output2$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output2$group <- factor(output2$group,levels = c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output2$Freq <- as.numeric(output2$Freq)
#   output2$label_y <- as.numeric(output2$label_y)
#   output2[output2$Freq < 3,"fill"] <- "" #20
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("pan- blueprint: ",sort(unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("STAD: ",sort(unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))])
#   
#   
#   
#   # cols <- c("RF: Proliferative cells" = "purple", "RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
#   p2 <- ggplot(output2,aes(x=group,y=Freq,fill=populations)) + 
#     geom_bar(stat="identity",colour="black") + 
#     scale_fill_manual(values=cols) +
#     # geom_text(aes(y=label_y,label=fill),size =10) +
#     theme_classic()
#   
#   
#   #### for p3
#   
#   # cell_summary_subset3 <- subset(cell_summary,cell_types_from_LJA %in% c("0_T memory cells","2_T memory cells")) #
#   # cell_summary_subset3$hpca_fine <- "CD4+ T"
#   # cell_summary_subset2 <- subset(cell_summary,cell_types_from_LJA %in% c("4_T cells"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   # cell_summary_subset2$hpca_fine <- "CD8+ T"#"CD4+ T"
#   cell_summary_subset3 <- subset(cell_summary,cell_types_from_LJA %in% c("NK"))#c("T_cell:CD4+_effector_memory","T_cell:CD4+","T_cell:CD4+_central_memory","T_cell:CD4+_Naive"))
#   cell_summary_subset3$hpca_fine <- "NK"#"CD4+ T"
#   # cell_summary_subset <- cell_summary_subset3 #rbind(cell_summary_subset3,cell_summary_subset2,cell_summary_subset3)
#   output_pre3 <- as.data.frame(table(cell_summary_subset3[,c(index[i],"hpca_fine")]))
#   output_pre3 <- data.frame(group=index[i],output_pre3)
#   output_pre3[,2] <- as.character(output_pre3[,2])
#   output_pre3$Freq <- as.numeric(output_pre3$Freq)
#   # output_pre <- rbind(output_pre,c("T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory","T_cell:CD4+_effector_memory",sum(output_pre$Freq)))
#   colnames(output_pre3)[2] <- "fill"
#   # output_pre3$fill <- as.character(new_cluster_names[output_pre3$fill])
#   output_pre3 <- data.frame(output_pre3,populations=output_pre3$fill)
#   output_pre3$fill <- as.character(t(as.data.frame(str_split(as.character(output_pre3$fill),": ",2)))[,2])
#   output_pre3$Freq <- as.numeric(output_pre3$Freq)
#   output_pre3 <- output_pre3[order(output_pre3$Freq),]
#   output_pre3 <- output_pre3 %>%
#     group_by(group) %>%
#     mutate(label_y = cumsum(Freq)-0.5*Freq)
#   output_pre3 <- output_pre3[order(-output_pre3$Freq),]
#   output3 <- rbind(output3,output_pre3)
#   
#   output3 <- as.data.frame(output3)
#   # output3 <- rbind(output3,c("RF","CD4+ T","CD4+ T",sum(output_pre3$Freq),"RF: CD4+ T",sum(output_pre3$Freq)/2))#;output3 <- output
#   # output <- rbind(output,c("RF","CD8+ T","CD8+ T",sum(output_pre$Freq),"RF: CD8+ T",sum(output_pre$Freq)/2)) ; output2 <- output#c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output3 <- rbind(output3,c("RF","NK cell","NK",sum(output_pre3$Freq),"RF: NK",sum(output_pre3$Freq)/2))#;output3<- output #c("original","CD4+ T","CD4+ T",sum(output_pre$Freq),sum(output_pre$Freq)/2))
#   # output <- rbind(output3,output2,output3)
#   output3$fill <- factor(output3$fill,levels = c(output3$fill))
#   output3$populations <- factor(output3$populations,levels = c(output3$populations))
#   # output3$group <- factor(output3$group,levels = c("RF","to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   # output3$group <- paste0("NK cells mapped to ",strsplit(as.character(unique(output3$group)),"_")[[1]][2])
#   # output3$group <- factor(output3$group,levels = c("NK cells mapped to TICA","NK cells mapped to panT","NK cells mapped to panblueprint","NK cells mapped to ABTC","NK cells mapped to HCC","NK cells mapped to STAD"))
#   output3$group <- factor(output3$group,levels = c("to_TICA","to_panT","to_panblueprint","to_ABTC","to_HCC","to_STAD"))
#   output3$Freq <- as.numeric(output3$Freq)
#   output3$label_y <- as.numeric(output3$label_y)
#   output3[output3$Freq < 3,"fill"] <- "" #30
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_TICA$predicted.to_TICA_T_NK_cells)))
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # names(cols) <- paste0(index2[i],sort(unique(LJA_T_NK_cells_panT$predicted.to_pancancer_T_cells_CD8_CD4)))
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("pan- blueprint: ",sort(unique(LJA_T_NK_cells_panblueprint$predicted.to_pancancer_blueprint_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("ABTC: ",sort(unique(LJA_T_NK_cells_ABTC$predicted.to_ABTC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("HCC: ",sort(unique(LJA_T_NK_cells_HCC$predicted.to_HCC_T_NK_cells)))])
#   # cols <- hue_pal()(length(unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))
#   # names(cols) <- as.character(new_cluster_names[paste0("STAD: ",sort(unique(LJA_T_NK_cells_STAD$predicted.to_STAD_T_NK_cells)))])
#   
#   
#   # cols <- c("RF: Proliferative cells" = "purple", "RF: CD4+ T" = "green","RF: CD8+ T" = "blue","RF: NK" = "red",cols)
#   p3 <- ggplot(output3,aes(x=group,y=Freq,fill=populations)) + 
#     geom_bar(stat="identity",colour="black") + 
#     scale_fill_manual(values=cols) +
#     # geom_text(aes(y=label_y,label=fill),size =10) +
#     theme_classic()
#   
#   
#   # ComLJAnePlots(plots=list(p1,p2,p3),nrow=1,legend="right")
#   p0 <- p0 + ylab ("T.cell (proliferative cells)") + theme(axis.title.y = element_text (size =30),axis.text.x = element_text(size = 30,color = "black"),axis.text.y = element_text(size = 18,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )
#   p1 <- p1 + ylab ("T.CD4 cells") + theme(axis.title.y = element_text (size =30),axis.text.x = element_text(size = 30,color = "black"),axis.text.y = element_text(size = 18,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )
#   p2 <- p2 + ylab ("T.CD8 cells") + theme(axis.title.y = element_text (size =30),axis.text.x = element_text(size = 30,color = "black"),axis.text.y = element_text(size = 18,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )
#   p3 <- p3 + ylab ("NK cells") + theme(axis.title.y = element_text (size =30), axis.text.x = element_text(size = 30,color = "black"),axis.text.y = element_text(size = 18,color = "black"),legend.text=element_text(size=20)) + guides(color = guide_legend(override.aes = list(size=14), ncol=1) )
#   CombinePlots(plots=list(p0,p1,p2,p3),nrow=1,legend="none")
#   # pdf(paste0(workdir,"test.pdf"))
#   # pp <- list(pp,ComLJAnePlots(plots=list(p1,p2,p3),nrow=1,legend="right"))
#   # dev.off()
# # }


  #### calculate the ARI

  library(Seurat)
  library(SeuratObject)
  library(lisi)
  library(mclust)
  library(cluster)
  to_TICA <- c("CD4","CD8","CD4","Proliferative","CD8","CD8","CD4")
  names(to_TICA) <- c("TICA: Regulatory T cells","TICA: Terminally exhausted CD8 T cells",
                      "TICA: Naive T cells","TICA: Proliferative T cells",           
                      "TICA: Pre-exhausted CD8 T cells","TICA: Cytotoxic CD8 T cells",           
                      "TICA: T helper cells")
  
  to_panT <- c("CD8","CD4","CD4","CD4","CD8","CD4","CD8","CD8")
  names(to_panT) <- c("pan- T: CD8.c10.Trm.ZNF683","pan- T: CD4.c10.Tm.CAPG","pan- T: CD4.c24.Mix.NME2",    
                      "pan- T: CD4.c03.Tn.ADSL","pan- T: CD8.c15.ISG.IFIT1","pan- T: CD4.c20.Treg.TNFRSF9",
                      "pan- T: CD8.c07.Temra.CX3CR1","pan- T: CD8.c06.Tem.GZMK" )
  
  to_panblueprint <- c("CD8","CD4","CD8","CD4","CD4","CD4","CD8","CD8","NK","NK")
  names(to_panblueprint) <- c("pan- blueprint: CD8 exhausted cytotoxic T","pan- blueprint: CD4 regulatory T",         
                              "pan- blueprint: CD8 pre-effector T","pan- blueprint: CD4 memory/effector T" ,   
                              "pan- blueprint: CD4 exhausted effector T","pan- blueprint: CD4 naive T",              
                              "pan- blueprint: CD8 effector T","pan- blueprint: CD8 memory T",             
                              "pan- blueprint: NK_XCL1","pan- blueprint: Cytotoxic NK" )
  
  to_ABTC <- c("CD8","CD4","CD8","CD4","CD4","NK","CD8","CD8")
  names(to_ABTC) <- c("ABTC: Exhausted CD8 T cell","ABTC: Naive CD4 T cell","ABTC: Memory CD8 T cell","ABTC: DNT",                 
                      "ABTC: Treg","ABTC: NCAM1nFCGR3Ap NK","ABTC: Naive CD8 T cell","ABTC: Effector CD8 T cell" )
  
  to_HCC <- c("CD4","CD8","NK","CD8","CD8","CD8","proliferative","CD8","NK","NK")
  names(to_HCC) <- c("HCC: 9_CD4_regulatory T","HCC: 14_CD8_TR_memory T","HCC: 25_NK","HCC: 2_CD8_MAIT",           
                     "HCC: 51_CD8_other T","HCC: 28_NK","HCC: 32_CD4_other T","HCC: 11_CD8_inter-states T",
                     "HCC: 18_NK","HCC: 8_NK")
  
  to_STAD <- c("CD8","CD4","CD8","proliferative","NK","CD4","CD8","CD4","NK")
  names(to_STAD) <- c("STAD: LNK3","STAD: Regulatory T","STAD: LNK2","STAD: Proliferative T", 
                      "STAD: LNK1","STAD: CD4 naive T cell","STAD: CD8 naive T","STAD: CD4 helper T",    
                     "STAD: LNK4" )
  
  # cell_summary_subset1 <- subset(cell_summary,cell_types_from_LJA %in% c("T.CD4")) #
  # cell_summary_subset1$hpca_fine <- "CD4+ T"
  # adjustedRandIndex(c(cell_summary_subset1$hpca_fine,"CD8+ T"),c(as.character(to_TICA[cell_summary_subset1$to_TICA]),"CD8+ T"))
  # adjustedRandIndex(c(cell_summary_subset1$hpca_fine,"CD8+ T"),c(as.character(to_panT[cell_summary_subset1$to_panT]),"CD8+ T"))
  # adjustedRandIndex(c(cell_summary_subset1$hpca_fine,"CD8+ T"),c(as.character(to_panblueprint[cell_summary_subset1$to_panblueprint]),"CD8+ T"))
  # adjustedRandIndex(c(cell_summary_subset1$hpca_fine,"CD8+ T"),c(as.character(to_ABTC[cell_summary_subset1$to_ABTC]),"CD8+ T"))
  # adjustedRandIndex(c(cell_summary_subset1$hpca_fine,"CD8+ T"),c(as.character(to_HCC[cell_summary_subset1$to_HCC]),"CD8+ T"))
  # adjustedRandIndex(c(cell_summary_subset1$hpca_fine,"CD8+ T"),c(as.character(to_STAD[cell_summary_subset1$to_STAD]),"CD8+ T"))
  
  
  
  adjustedRandIndex(cell_summary$cell_types_from_LJA,as.character(to_TICA[cell_summary$to_TICA]))
  adjustedRandIndex(cell_summary$cell_types_from_LJA,as.character(to_panT[cell_summary$to_panT]))
  adjustedRandIndex(cell_summary$cell_types_from_LJA,as.character(to_panblueprint[cell_summary$to_panblueprint]))
  adjustedRandIndex(cell_summary$cell_types_from_LJA,as.character(to_ABTC[cell_summary$to_ABTC]))
  adjustedRandIndex(cell_summary$cell_types_from_LJA,as.character(to_HCC[cell_summary$to_HCC]))
  adjustedRandIndex(cell_summary$cell_types_from_LJA,as.character(to_STAD[cell_summary$to_STAD]))

