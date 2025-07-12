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

####################################################
#### 2024.02.28 using Livnat Jerby-Arnon
####
#### 2025.05.05 Livnat Jerby-Arnon with granular cell annotations
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




