##############################################
#### 2023.07.06 plot marker gene redetection rates by groups
###############################################
library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
library(ggplot2)
library(ggpubr)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
output <- read.table(paste0(workdir,"3_mapping_query_datasets_criteria_redetections_and_correlations_within_populations_20240319.txt"),sep="\t",header=T,check.names = F)

CD4CD8_proliferative<- c("TICA: Proliferative T cells",
                         "HCC: 32_CD4_other T",
                         "HCC: 51_CD8_other T",
                         "STAD: Proliferative T")
CD4_regulatory <- c("TICA: Regulatory T cells",
                    "pan- blueprint: CD4 regulatory T",
                    "ABTC: Treg",
                    "HCC: 9_CD4_regulatory T",
                    "STAD: Regulatory T")
CD4_naive <- c("TICA: Naive T cells",
               "pan- blueprint: CD4 naive T",
               "STAD: CD4 naive T cell",
               "ABTC: Naive CD4 T cell",
               "HCC: 20_CD4_naive T",
               "HCC: 35_CD4_central memory T")
CD4_memory <- c("pan- blueprint: CD4 memory/effector T",
                "TICA: Naive-memory CD4 T cells",
                "STAD: CD8 effector T cell",
                "TICA: Transitional memory CD4 T cells",
                "ABTC: DNT")
CD4_exhausted <- c("TICA: T helper cells",
                   "pan- blueprint: CD4 exhausted effector T",
                   "STAD: CD4 helper T")
CD8_exhausted <- c("ABTC: Exhausted CD8 T cell",
                   "pan- blueprint: CD8 exhausted cytotoxic T",
                   "TICA: Terminally exhausted CD8 T cells",
                   "STAD: LNK2"
)
CD8_pre_exhausted <- "TICA: Pre-exhausted CD8 T cells"
CD8_effector <- c("HCC: 28_NK",
                  "ABTC: Effector CD8 T cell",
                  "pan- blueprint: CD8 effector T",
                  "HCC: 13_CD8_effector memory T",
                  "ABTC: NKT and NK cell",
                  "STAD: LNK3")
CD8_pre_effector <- c("TICA: Cytotoxic CD8 T cells",
                      "HCC: 11_CD8_inter-states T",
                      "pan- blueprint: CD8 pre-effector T",
                      "ABTC: Memory CD8 T cell",
                      "HCC: 1_CD8_cytotoxic T")
CD8_memory <- c("pan- blueprint: CD8 memory T",
                "TICA: Effector memory CD8 T cells",
                "TICA: Recently activated CD4 T cells"
)
CD8_naive <- c("ABTC: Naive CD8 T cell",
               "HCC: 2_CD8_MAIT",
               "HCC: 14_CD8_TR_memory T",
               "TICA: Th17 cells",
               "STAD: CD8 naive T"
)
NK_pn <- c("ABTC: NCAM1pFCGR3An NK",
           "HCC: 8_NK",
           "pan- blueprint: NK_XCL1"
)
NK_pp <- c("ABTC: NCAM1pFCGR3Ap NK",
           "HCC: 25_NK",
           "STAD: LNK1",
           "TICA: NK",
           "pan- blueprint: Cytotoxic NK"
)
NK_np <- "ABTC: NCAM1nFCGR3Ap NK"
NK_nn <- c("HCC: 18_NK",
           "STAD: LNK4"
)

plotdata <- data.frame(
groups = c(rep("CD4/CD8 Prolifereative T cells",length(CD4CD8_proliferative)),
rep("CD4 Regulatory T cells",length(CD4_regulatory)),
rep("CD4 naive T cells",length(CD4_naive)),
rep("CD4 memory T cells",length(CD4_memory)),
rep("CD4 exhausted T cells",length(CD4_exhausted)),
rep("CD8 exhausted T cells",length(CD8_exhausted)),
rep("CD8 pre-exhausted T cells",length(CD8_pre_exhausted)),
rep("CD8 effector T cells",length(CD8_effector)),
rep("CD8 pre-effector T cells",length(CD8_pre_effector)),
rep("CD8 memory T cells",length(CD8_memory)),
rep("CD8 naive T cells",length(CD8_naive)),
rep("NCAM1+FCGR3A- NK",length(NK_pn)),
rep("NCAM1+FCGR3A+ NK",length(NK_pp)),
rep("NCAM1-FCGR3A+ NK",length(NK_np)),
rep("NCAM1-FCGR3A- NK",length(NK_nn))),
values = apply(output[1:4,c(CD4CD8_proliferative,CD4_regulatory,CD4_naive,CD4_memory,CD4_exhausted,CD8_exhausted,CD8_pre_exhausted,CD8_effector,CD8_pre_effector,CD8_memory,CD8_naive,NK_pn,NK_pp,NK_np,NK_nn)],
               2,function(x) mean(x,na.rm=T)))
plotdata$groups <- factor(plotdata$groups,levels = unique(plotdata$groups))

# plotdata <- data.frame(
#   groups = c(rep("group 1",3),rep("group 2", 7),rep("group 3",9),rep("group 4", 9),rep("group 5",4),rep("group 6", 8),
#              rep("group 7",3),rep("group 8",3),rep("group 9",5)),
#   values = c(0.417371511,
#              0.556679968,
#              0.442217768,
#              0.23559765,
#              0.256162381,
#              0.427069892,
#              0.381433561,
#              0.37916154,
#              0.353081162,
#              0.48493833,
#              0.53506165,
#              0.418253491,
#              0.33170624,
#              0.439842052,
#              0.432844282,
#              0.319523958,
#              0.195086477,
#              0.404843438,
#              0.331427563,
#              0.429578446,
#              0.266907483,
#              0.394484373,
#              0.222873957,
#              0.19915843,
#              0.20484672,
#              0.260776357,
#              0.391743625,
#              0.2696559,
#              0.244274809,
#              0.193619699,
#              0.266869773,
#              0.307706961,
#              0.36590922,
#              0.156414208,
#              0.25994974,
#              0.261319168,
#              0.411883223,
#              0.152619555,
#              0.387519105,
#              0.241267847,
#              0.249231211,
#              0.210417761,
#              0.347399556,
#              0.291300873,
#              0.162922796,
#              0.302198153,
#              0.324036695,
#              0.457787836,
#              0.5212776,
#              0.339045054,
#              0.327849293)
# )
# plotdata$groups <- as.factor(plotdata$groups)
my_comparisons = list( #c("group 1", "group 4"),c("group 1", "group 5"), c("group 1", "group 6"), c("group 1", "group 7"),c("group 1", "group 8"),c("group 1", "group 9"),
                       c("group 5","group 9"),c("group 6","group 9") )
ggplot(plotdata, aes(x = groups, y = values)) + 
  geom_boxplot(outlier.color = "red",fill="#FFC000") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 28,color = "black",angle = 270),
        axis.text.y = element_text(size = 28,color = "black"),
        axis.title.y = element_text(size = 28),
        legend.text=element_text(size=20)) + 
  guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="Marker gene redetection rates",size =28) +
  stat_compare_means(method="t.test",comparisons = my_comparisons, label.y = c(0.65, 0.75),bracket.size = 1, size=10) +#, 0.85,0.95,1.05,1.15,1.25,1.35))+
  stat_compare_means(label.y = 1)


plot_ove <- read.table(file=paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_ove_plot_ordered_20230606.txt"),header=T,sep="\t")
plot_ove[lower.tri(plot_ove, diag = TRUE)] <- NA
plot_ove_data <- data.frame(
  groups = c(rep("group 1",(4*4-4)/2),rep("group 2", (7*7-7)/2),rep("group 3",(10*10-10)/2),rep("group 4", (9*9-9)/2),rep("group 5",(4*4-4)/2),rep("group 6", (8*8-8)/2),
             rep("group 7",(3*3-3)/2),rep("group 8",(3*3-3)/2),rep("group 9",(5*5-5)/2)),
  values = c(as.numeric(na.omit(as.numeric(unlist(plot_ove[1:4,1:4])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[5:11,5:11])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[12:21,12:21])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[22:30,22:30])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[31:34,31:34])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[36:43,35:43])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[44:46,44:46])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[47:49,47:49])))),
             as.numeric(na.omit(as.numeric(unlist(plot_ove[52:56,52:56])))))
)
my_comparisons = list( c("group 1", "group 4"),c("group 1", "group 5"), c("group 1", "group 6"), c("group 1", "group 7"),c("group 1", "group 8"),c("group 1", "group 9"),c("group 5","group 9") )
ggplot(plot_ove_data, aes(x = groups, y = values)) + 
  geom_boxplot(outlier.color = "red",fill="#FFC000") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 28,color = "black",angle = 270),
        axis.text.y = element_text(size = 28,color = "black"),
        legend.text=element_text(size=20)) + 
  guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ labs(y="") +
  stat_compare_means(method="t.test",comparisons = my_comparisons, label.y = c(0.65, 0.75, 0.85,0.95,1.05,1.15))+
  stat_compare_means(label.y = 1)
##############################################
#### 2023.07.03 check mapping score
###############################################
library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"


querydata <- c("Kathryn","Baolin","Moshe","Caushi")
atlases <- c("TICA","pancancer_blueprint","ABTC","HCC","STAD")
output <- matrix(NA,length(querydata),ncol=length(atlases))
rownames(output) <- querydata
colnames(output) <- atlases
for (i in 1:length(querydata)) {
  print(i)
  for (j in 1:length(atlases)) {
    print(j)
    cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_",querydata[i],"_T_NK_cells_to_",atlases[j],"_T_NK_cells_mapscore.rds"))
    output[i,j] <- mean(unlist(cells1@meta.data[paste0("mapscore_to_",atlases[j],"_T_NK_cells")]))
    rm(cells1)
  }
}






##############################################
#### 2023.01.30 check different resolution
###############################################
library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.1, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.2, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.3, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.4, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.5, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.6, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.7, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.8, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 0.9, n.start = 10)
Kathryn_T_NK_cells1 <- FindClusters(Kathryn_T_NK_cells1, resolution = 1, n.start = 10)
saveRDS(Kathryn_T_NK_cells1, file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
Kathryn_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Kathryn_T_NK_cells2$RNA_snn_res.0.1 <- Kathryn_T_NK_cells1$RNA_snn_res.0.1
Kathryn_T_NK_cells2$RNA_snn_res.0.2 <- Kathryn_T_NK_cells1$RNA_snn_res.0.2
Kathryn_T_NK_cells2$RNA_snn_res.0.3 <- Kathryn_T_NK_cells1$RNA_snn_res.0.3
Kathryn_T_NK_cells2$RNA_snn_res.0.4 <- Kathryn_T_NK_cells1$RNA_snn_res.0.4
Kathryn_T_NK_cells2$RNA_snn_res.0.5 <- Kathryn_T_NK_cells1$RNA_snn_res.0.5
Kathryn_T_NK_cells2$RNA_snn_res.0.6 <- Kathryn_T_NK_cells1$RNA_snn_res.0.6
Kathryn_T_NK_cells2$RNA_snn_res.0.7 <- Kathryn_T_NK_cells1$RNA_snn_res.0.7
Kathryn_T_NK_cells2$RNA_snn_res.0.8 <- Kathryn_T_NK_cells1$RNA_snn_res.0.8
Kathryn_T_NK_cells2$RNA_snn_res.0.9 <- Kathryn_T_NK_cells1$RNA_snn_res.0.9
Kathryn_T_NK_cells2$RNA_snn_res.1 <- Kathryn_T_NK_cells1$RNA_snn_res.1
saveRDS(Kathryn_T_NK_cells2, file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
rm(Kathryn_T_NK_cells2)
Kathryn_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Kathryn_T_NK_cells3$RNA_snn_res.0.1 <- Kathryn_T_NK_cells1$RNA_snn_res.0.1
Kathryn_T_NK_cells3$RNA_snn_res.0.2 <- Kathryn_T_NK_cells1$RNA_snn_res.0.2
Kathryn_T_NK_cells3$RNA_snn_res.0.3 <- Kathryn_T_NK_cells1$RNA_snn_res.0.3
Kathryn_T_NK_cells3$RNA_snn_res.0.4 <- Kathryn_T_NK_cells1$RNA_snn_res.0.4
Kathryn_T_NK_cells3$RNA_snn_res.0.5 <- Kathryn_T_NK_cells1$RNA_snn_res.0.5
Kathryn_T_NK_cells3$RNA_snn_res.0.6 <- Kathryn_T_NK_cells1$RNA_snn_res.0.6
Kathryn_T_NK_cells3$RNA_snn_res.0.7 <- Kathryn_T_NK_cells1$RNA_snn_res.0.7
Kathryn_T_NK_cells3$RNA_snn_res.0.8 <- Kathryn_T_NK_cells1$RNA_snn_res.0.8
Kathryn_T_NK_cells3$RNA_snn_res.0.9 <- Kathryn_T_NK_cells1$RNA_snn_res.0.9
Kathryn_T_NK_cells3$RNA_snn_res.1 <- Kathryn_T_NK_cells1$RNA_snn_res.1
saveRDS(Kathryn_T_NK_cells3, file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
rm(Kathryn_T_NK_cells3)
Kathryn_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Kathryn_T_NK_cells4$RNA_snn_res.0.1 <- Kathryn_T_NK_cells1$RNA_snn_res.0.1
Kathryn_T_NK_cells4$RNA_snn_res.0.2 <- Kathryn_T_NK_cells1$RNA_snn_res.0.2
Kathryn_T_NK_cells4$RNA_snn_res.0.3 <- Kathryn_T_NK_cells1$RNA_snn_res.0.3
Kathryn_T_NK_cells4$RNA_snn_res.0.4 <- Kathryn_T_NK_cells1$RNA_snn_res.0.4
Kathryn_T_NK_cells4$RNA_snn_res.0.5 <- Kathryn_T_NK_cells1$RNA_snn_res.0.5
Kathryn_T_NK_cells4$RNA_snn_res.0.6 <- Kathryn_T_NK_cells1$RNA_snn_res.0.6
Kathryn_T_NK_cells4$RNA_snn_res.0.7 <- Kathryn_T_NK_cells1$RNA_snn_res.0.7
Kathryn_T_NK_cells4$RNA_snn_res.0.8 <- Kathryn_T_NK_cells1$RNA_snn_res.0.8
Kathryn_T_NK_cells4$RNA_snn_res.0.9 <- Kathryn_T_NK_cells1$RNA_snn_res.0.9
Kathryn_T_NK_cells4$RNA_snn_res.1 <- Kathryn_T_NK_cells1$RNA_snn_res.1
saveRDS(Kathryn_T_NK_cells4, file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
rm(Kathryn_T_NK_cells4)
Kathryn_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells.rds"))
Kathryn_T_NK_cells5$RNA_snn_res.0.1 <- Kathryn_T_NK_cells1$RNA_snn_res.0.1
Kathryn_T_NK_cells5$RNA_snn_res.0.2 <- Kathryn_T_NK_cells1$RNA_snn_res.0.2
Kathryn_T_NK_cells5$RNA_snn_res.0.3 <- Kathryn_T_NK_cells1$RNA_snn_res.0.3
Kathryn_T_NK_cells5$RNA_snn_res.0.4 <- Kathryn_T_NK_cells1$RNA_snn_res.0.4
Kathryn_T_NK_cells5$RNA_snn_res.0.5 <- Kathryn_T_NK_cells1$RNA_snn_res.0.5
Kathryn_T_NK_cells5$RNA_snn_res.0.6 <- Kathryn_T_NK_cells1$RNA_snn_res.0.6
Kathryn_T_NK_cells5$RNA_snn_res.0.7 <- Kathryn_T_NK_cells1$RNA_snn_res.0.7
Kathryn_T_NK_cells5$RNA_snn_res.0.8 <- Kathryn_T_NK_cells1$RNA_snn_res.0.8
Kathryn_T_NK_cells5$RNA_snn_res.0.9 <- Kathryn_T_NK_cells1$RNA_snn_res.0.9
Kathryn_T_NK_cells5$RNA_snn_res.1 <- Kathryn_T_NK_cells1$RNA_snn_res.1
saveRDS(Kathryn_T_NK_cells5, file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells.rds"))
rm(Kathryn_T_NK_cells5)
Kathryn_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells.rds"))
Kathryn_T_NK_cells6$RNA_snn_res.0.1 <- Kathryn_T_NK_cells1$RNA_snn_res.0.1
Kathryn_T_NK_cells6$RNA_snn_res.0.2 <- Kathryn_T_NK_cells1$RNA_snn_res.0.2
Kathryn_T_NK_cells6$RNA_snn_res.0.3 <- Kathryn_T_NK_cells1$RNA_snn_res.0.3
Kathryn_T_NK_cells6$RNA_snn_res.0.4 <- Kathryn_T_NK_cells1$RNA_snn_res.0.4
Kathryn_T_NK_cells6$RNA_snn_res.0.5 <- Kathryn_T_NK_cells1$RNA_snn_res.0.5
Kathryn_T_NK_cells6$RNA_snn_res.0.6 <- Kathryn_T_NK_cells1$RNA_snn_res.0.6
Kathryn_T_NK_cells6$RNA_snn_res.0.7 <- Kathryn_T_NK_cells1$RNA_snn_res.0.7
Kathryn_T_NK_cells6$RNA_snn_res.0.8 <- Kathryn_T_NK_cells1$RNA_snn_res.0.8
Kathryn_T_NK_cells6$RNA_snn_res.0.9 <- Kathryn_T_NK_cells1$RNA_snn_res.0.9
Kathryn_T_NK_cells6$RNA_snn_res.1 <- Kathryn_T_NK_cells1$RNA_snn_res.1
saveRDS(Kathryn_T_NK_cells6, file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells.rds"))
rm(Kathryn_T_NK_cells6)


Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.1, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.2, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.3, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.4, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.5, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.6, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.7, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.8, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 0.9, n.start = 10)
Baolin_T_NK_cells1 <- FindClusters(Baolin_T_NK_cells1, resolution = 1, n.start = 10)
saveRDS(Baolin_T_NK_cells1, file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
Baolin_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Baolin_T_NK_cells2$RNA_snn_res.0.1 <- Baolin_T_NK_cells1$RNA_snn_res.0.1
Baolin_T_NK_cells2$RNA_snn_res.0.2 <- Baolin_T_NK_cells1$RNA_snn_res.0.2
Baolin_T_NK_cells2$RNA_snn_res.0.3 <- Baolin_T_NK_cells1$RNA_snn_res.0.3
Baolin_T_NK_cells2$RNA_snn_res.0.4 <- Baolin_T_NK_cells1$RNA_snn_res.0.4
Baolin_T_NK_cells2$RNA_snn_res.0.5 <- Baolin_T_NK_cells1$RNA_snn_res.0.5
Baolin_T_NK_cells2$RNA_snn_res.0.6 <- Baolin_T_NK_cells1$RNA_snn_res.0.6
Baolin_T_NK_cells2$RNA_snn_res.0.7 <- Baolin_T_NK_cells1$RNA_snn_res.0.7
Baolin_T_NK_cells2$RNA_snn_res.0.8 <- Baolin_T_NK_cells1$RNA_snn_res.0.8
Baolin_T_NK_cells2$RNA_snn_res.0.9 <- Baolin_T_NK_cells1$RNA_snn_res.0.9
Baolin_T_NK_cells2$RNA_snn_res.1 <- Baolin_T_NK_cells1$RNA_snn_res.1
saveRDS(Baolin_T_NK_cells2, file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
rm(Baolin_T_NK_cells2)
Baolin_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Baolin_T_NK_cells3$RNA_snn_res.0.1 <- Baolin_T_NK_cells1$RNA_snn_res.0.1
Baolin_T_NK_cells3$RNA_snn_res.0.2 <- Baolin_T_NK_cells1$RNA_snn_res.0.2
Baolin_T_NK_cells3$RNA_snn_res.0.3 <- Baolin_T_NK_cells1$RNA_snn_res.0.3
Baolin_T_NK_cells3$RNA_snn_res.0.4 <- Baolin_T_NK_cells1$RNA_snn_res.0.4
Baolin_T_NK_cells3$RNA_snn_res.0.5 <- Baolin_T_NK_cells1$RNA_snn_res.0.5
Baolin_T_NK_cells3$RNA_snn_res.0.6 <- Baolin_T_NK_cells1$RNA_snn_res.0.6
Baolin_T_NK_cells3$RNA_snn_res.0.7 <- Baolin_T_NK_cells1$RNA_snn_res.0.7
Baolin_T_NK_cells3$RNA_snn_res.0.8 <- Baolin_T_NK_cells1$RNA_snn_res.0.8
Baolin_T_NK_cells3$RNA_snn_res.0.9 <- Baolin_T_NK_cells1$RNA_snn_res.0.9
Baolin_T_NK_cells3$RNA_snn_res.1 <- Baolin_T_NK_cells1$RNA_snn_res.1
saveRDS(Baolin_T_NK_cells3, file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
rm(Baolin_T_NK_cells3)
Baolin_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Baolin_T_NK_cells4$RNA_snn_res.0.1 <- Baolin_T_NK_cells1$RNA_snn_res.0.1
Baolin_T_NK_cells4$RNA_snn_res.0.2 <- Baolin_T_NK_cells1$RNA_snn_res.0.2
Baolin_T_NK_cells4$RNA_snn_res.0.3 <- Baolin_T_NK_cells1$RNA_snn_res.0.3
Baolin_T_NK_cells4$RNA_snn_res.0.4 <- Baolin_T_NK_cells1$RNA_snn_res.0.4
Baolin_T_NK_cells4$RNA_snn_res.0.5 <- Baolin_T_NK_cells1$RNA_snn_res.0.5
Baolin_T_NK_cells4$RNA_snn_res.0.6 <- Baolin_T_NK_cells1$RNA_snn_res.0.6
Baolin_T_NK_cells4$RNA_snn_res.0.7 <- Baolin_T_NK_cells1$RNA_snn_res.0.7
Baolin_T_NK_cells4$RNA_snn_res.0.8 <- Baolin_T_NK_cells1$RNA_snn_res.0.8
Baolin_T_NK_cells4$RNA_snn_res.0.9 <- Baolin_T_NK_cells1$RNA_snn_res.0.9
Baolin_T_NK_cells4$RNA_snn_res.1 <- Baolin_T_NK_cells1$RNA_snn_res.1
saveRDS(Baolin_T_NK_cells4, file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
rm(Baolin_T_NK_cells4)
Baolin_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"))
Baolin_T_NK_cells5$RNA_snn_res.0.1 <- Baolin_T_NK_cells1$RNA_snn_res.0.1
Baolin_T_NK_cells5$RNA_snn_res.0.2 <- Baolin_T_NK_cells1$RNA_snn_res.0.2
Baolin_T_NK_cells5$RNA_snn_res.0.3 <- Baolin_T_NK_cells1$RNA_snn_res.0.3
Baolin_T_NK_cells5$RNA_snn_res.0.4 <- Baolin_T_NK_cells1$RNA_snn_res.0.4
Baolin_T_NK_cells5$RNA_snn_res.0.5 <- Baolin_T_NK_cells1$RNA_snn_res.0.5
Baolin_T_NK_cells5$RNA_snn_res.0.6 <- Baolin_T_NK_cells1$RNA_snn_res.0.6
Baolin_T_NK_cells5$RNA_snn_res.0.7 <- Baolin_T_NK_cells1$RNA_snn_res.0.7
Baolin_T_NK_cells5$RNA_snn_res.0.8 <- Baolin_T_NK_cells1$RNA_snn_res.0.8
Baolin_T_NK_cells5$RNA_snn_res.0.9 <- Baolin_T_NK_cells1$RNA_snn_res.0.9
Baolin_T_NK_cells5$RNA_snn_res.1 <- Baolin_T_NK_cells1$RNA_snn_res.1
saveRDS(Baolin_T_NK_cells5, file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"))
rm(Baolin_T_NK_cells5)
Baolin_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"))
Baolin_T_NK_cells6$RNA_snn_res.0.1 <- Baolin_T_NK_cells1$RNA_snn_res.0.1
Baolin_T_NK_cells6$RNA_snn_res.0.2 <- Baolin_T_NK_cells1$RNA_snn_res.0.2
Baolin_T_NK_cells6$RNA_snn_res.0.3 <- Baolin_T_NK_cells1$RNA_snn_res.0.3
Baolin_T_NK_cells6$RNA_snn_res.0.4 <- Baolin_T_NK_cells1$RNA_snn_res.0.4
Baolin_T_NK_cells6$RNA_snn_res.0.5 <- Baolin_T_NK_cells1$RNA_snn_res.0.5
Baolin_T_NK_cells6$RNA_snn_res.0.6 <- Baolin_T_NK_cells1$RNA_snn_res.0.6
Baolin_T_NK_cells6$RNA_snn_res.0.7 <- Baolin_T_NK_cells1$RNA_snn_res.0.7
Baolin_T_NK_cells6$RNA_snn_res.0.8 <- Baolin_T_NK_cells1$RNA_snn_res.0.8
Baolin_T_NK_cells6$RNA_snn_res.0.9 <- Baolin_T_NK_cells1$RNA_snn_res.0.9
Baolin_T_NK_cells6$RNA_snn_res.1 <- Baolin_T_NK_cells1$RNA_snn_res.1
saveRDS(Baolin_T_NK_cells6, file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"))
rm(Baolin_T_NK_cells6)



Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.1, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.2, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.3, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.4, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.5, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.6, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.7, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.8, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 0.9, n.start = 10)
Moshe_T_NK_cells1 <- FindClusters(Moshe_T_NK_cells1, resolution = 1, n.start = 10)
saveRDS(Moshe_T_NK_cells1, file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Moshe_T_NK_cells2$RNA_snn_res.0.1 <- Moshe_T_NK_cells1$RNA_snn_res.0.1
Moshe_T_NK_cells2$RNA_snn_res.0.2 <- Moshe_T_NK_cells1$RNA_snn_res.0.2
Moshe_T_NK_cells2$RNA_snn_res.0.3 <- Moshe_T_NK_cells1$RNA_snn_res.0.3
Moshe_T_NK_cells2$RNA_snn_res.0.4 <- Moshe_T_NK_cells1$RNA_snn_res.0.4
Moshe_T_NK_cells2$RNA_snn_res.0.5 <- Moshe_T_NK_cells1$RNA_snn_res.0.5
Moshe_T_NK_cells2$RNA_snn_res.0.6 <- Moshe_T_NK_cells1$RNA_snn_res.0.6
Moshe_T_NK_cells2$RNA_snn_res.0.7 <- Moshe_T_NK_cells1$RNA_snn_res.0.7
Moshe_T_NK_cells2$RNA_snn_res.0.8 <- Moshe_T_NK_cells1$RNA_snn_res.0.8
Moshe_T_NK_cells2$RNA_snn_res.0.9 <- Moshe_T_NK_cells1$RNA_snn_res.0.9
Moshe_T_NK_cells2$RNA_snn_res.1 <- Moshe_T_NK_cells1$RNA_snn_res.1
saveRDS(Moshe_T_NK_cells2, file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
rm(Moshe_T_NK_cells2)
Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Moshe_T_NK_cells3$RNA_snn_res.0.1 <- Moshe_T_NK_cells1$RNA_snn_res.0.1
Moshe_T_NK_cells3$RNA_snn_res.0.2 <- Moshe_T_NK_cells1$RNA_snn_res.0.2
Moshe_T_NK_cells3$RNA_snn_res.0.3 <- Moshe_T_NK_cells1$RNA_snn_res.0.3
Moshe_T_NK_cells3$RNA_snn_res.0.4 <- Moshe_T_NK_cells1$RNA_snn_res.0.4
Moshe_T_NK_cells3$RNA_snn_res.0.5 <- Moshe_T_NK_cells1$RNA_snn_res.0.5
Moshe_T_NK_cells3$RNA_snn_res.0.6 <- Moshe_T_NK_cells1$RNA_snn_res.0.6
Moshe_T_NK_cells3$RNA_snn_res.0.7 <- Moshe_T_NK_cells1$RNA_snn_res.0.7
Moshe_T_NK_cells3$RNA_snn_res.0.8 <- Moshe_T_NK_cells1$RNA_snn_res.0.8
Moshe_T_NK_cells3$RNA_snn_res.0.9 <- Moshe_T_NK_cells1$RNA_snn_res.0.9
Moshe_T_NK_cells3$RNA_snn_res.1 <- Moshe_T_NK_cells1$RNA_snn_res.1
saveRDS(Moshe_T_NK_cells3, file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
rm(Moshe_T_NK_cells3)
Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Moshe_T_NK_cells4$RNA_snn_res.0.1 <- Moshe_T_NK_cells1$RNA_snn_res.0.1
Moshe_T_NK_cells4$RNA_snn_res.0.2 <- Moshe_T_NK_cells1$RNA_snn_res.0.2
Moshe_T_NK_cells4$RNA_snn_res.0.3 <- Moshe_T_NK_cells1$RNA_snn_res.0.3
Moshe_T_NK_cells4$RNA_snn_res.0.4 <- Moshe_T_NK_cells1$RNA_snn_res.0.4
Moshe_T_NK_cells4$RNA_snn_res.0.5 <- Moshe_T_NK_cells1$RNA_snn_res.0.5
Moshe_T_NK_cells4$RNA_snn_res.0.6 <- Moshe_T_NK_cells1$RNA_snn_res.0.6
Moshe_T_NK_cells4$RNA_snn_res.0.7 <- Moshe_T_NK_cells1$RNA_snn_res.0.7
Moshe_T_NK_cells4$RNA_snn_res.0.8 <- Moshe_T_NK_cells1$RNA_snn_res.0.8
Moshe_T_NK_cells4$RNA_snn_res.0.9 <- Moshe_T_NK_cells1$RNA_snn_res.0.9
Moshe_T_NK_cells4$RNA_snn_res.1 <- Moshe_T_NK_cells1$RNA_snn_res.1
saveRDS(Moshe_T_NK_cells4, file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
rm(Moshe_T_NK_cells4)
Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
Moshe_T_NK_cells5$RNA_snn_res.0.1 <- Moshe_T_NK_cells1$RNA_snn_res.0.1
Moshe_T_NK_cells5$RNA_snn_res.0.2 <- Moshe_T_NK_cells1$RNA_snn_res.0.2
Moshe_T_NK_cells5$RNA_snn_res.0.3 <- Moshe_T_NK_cells1$RNA_snn_res.0.3
Moshe_T_NK_cells5$RNA_snn_res.0.4 <- Moshe_T_NK_cells1$RNA_snn_res.0.4
Moshe_T_NK_cells5$RNA_snn_res.0.5 <- Moshe_T_NK_cells1$RNA_snn_res.0.5
Moshe_T_NK_cells5$RNA_snn_res.0.6 <- Moshe_T_NK_cells1$RNA_snn_res.0.6
Moshe_T_NK_cells5$RNA_snn_res.0.7 <- Moshe_T_NK_cells1$RNA_snn_res.0.7
Moshe_T_NK_cells5$RNA_snn_res.0.8 <- Moshe_T_NK_cells1$RNA_snn_res.0.8
Moshe_T_NK_cells5$RNA_snn_res.0.9 <- Moshe_T_NK_cells1$RNA_snn_res.0.9
Moshe_T_NK_cells5$RNA_snn_res.1 <- Moshe_T_NK_cells1$RNA_snn_res.1
saveRDS(Moshe_T_NK_cells5, file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
rm(Moshe_T_NK_cells5)
Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
Moshe_T_NK_cells6$RNA_snn_res.0.1 <- Moshe_T_NK_cells1$RNA_snn_res.0.1
Moshe_T_NK_cells6$RNA_snn_res.0.2 <- Moshe_T_NK_cells1$RNA_snn_res.0.2
Moshe_T_NK_cells6$RNA_snn_res.0.3 <- Moshe_T_NK_cells1$RNA_snn_res.0.3
Moshe_T_NK_cells6$RNA_snn_res.0.4 <- Moshe_T_NK_cells1$RNA_snn_res.0.4
Moshe_T_NK_cells6$RNA_snn_res.0.5 <- Moshe_T_NK_cells1$RNA_snn_res.0.5
Moshe_T_NK_cells6$RNA_snn_res.0.6 <- Moshe_T_NK_cells1$RNA_snn_res.0.6
Moshe_T_NK_cells6$RNA_snn_res.0.7 <- Moshe_T_NK_cells1$RNA_snn_res.0.7
Moshe_T_NK_cells6$RNA_snn_res.0.8 <- Moshe_T_NK_cells1$RNA_snn_res.0.8
Moshe_T_NK_cells6$RNA_snn_res.0.9 <- Moshe_T_NK_cells1$RNA_snn_res.0.9
Moshe_T_NK_cells6$RNA_snn_res.1 <- Moshe_T_NK_cells1$RNA_snn_res.1
saveRDS(Moshe_T_NK_cells6, file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
rm(Moshe_T_NK_cells6)




Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.1, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.2, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.3, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.4, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.5, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.6, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.7, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.8, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 0.9, n.start = 10)
Caushi_T_NK_cells1 <- FindClusters(Caushi_T_NK_cells1, resolution = 1, n.start = 10)
saveRDS(Caushi_T_NK_cells1, file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
Caushi_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Caushi_T_NK_cells2$RNA_snn_res.0.1 <- Caushi_T_NK_cells1$RNA_snn_res.0.1
Caushi_T_NK_cells2$RNA_snn_res.0.2 <- Caushi_T_NK_cells1$RNA_snn_res.0.2
Caushi_T_NK_cells2$RNA_snn_res.0.3 <- Caushi_T_NK_cells1$RNA_snn_res.0.3
Caushi_T_NK_cells2$RNA_snn_res.0.4 <- Caushi_T_NK_cells1$RNA_snn_res.0.4
Caushi_T_NK_cells2$RNA_snn_res.0.5 <- Caushi_T_NK_cells1$RNA_snn_res.0.5
Caushi_T_NK_cells2$RNA_snn_res.0.6 <- Caushi_T_NK_cells1$RNA_snn_res.0.6
Caushi_T_NK_cells2$RNA_snn_res.0.7 <- Caushi_T_NK_cells1$RNA_snn_res.0.7
Caushi_T_NK_cells2$RNA_snn_res.0.8 <- Caushi_T_NK_cells1$RNA_snn_res.0.8
Caushi_T_NK_cells2$RNA_snn_res.0.9 <- Caushi_T_NK_cells1$RNA_snn_res.0.9
Caushi_T_NK_cells2$RNA_snn_res.1 <- Caushi_T_NK_cells1$RNA_snn_res.1
saveRDS(Caushi_T_NK_cells2, file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
rm(Caushi_T_NK_cells2)
Caushi_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Caushi_T_NK_cells3$RNA_snn_res.0.1 <- Caushi_T_NK_cells1$RNA_snn_res.0.1
Caushi_T_NK_cells3$RNA_snn_res.0.2 <- Caushi_T_NK_cells1$RNA_snn_res.0.2
Caushi_T_NK_cells3$RNA_snn_res.0.3 <- Caushi_T_NK_cells1$RNA_snn_res.0.3
Caushi_T_NK_cells3$RNA_snn_res.0.4 <- Caushi_T_NK_cells1$RNA_snn_res.0.4
Caushi_T_NK_cells3$RNA_snn_res.0.5 <- Caushi_T_NK_cells1$RNA_snn_res.0.5
Caushi_T_NK_cells3$RNA_snn_res.0.6 <- Caushi_T_NK_cells1$RNA_snn_res.0.6
Caushi_T_NK_cells3$RNA_snn_res.0.7 <- Caushi_T_NK_cells1$RNA_snn_res.0.7
Caushi_T_NK_cells3$RNA_snn_res.0.8 <- Caushi_T_NK_cells1$RNA_snn_res.0.8
Caushi_T_NK_cells3$RNA_snn_res.0.9 <- Caushi_T_NK_cells1$RNA_snn_res.0.9
Caushi_T_NK_cells3$RNA_snn_res.1 <- Caushi_T_NK_cells1$RNA_snn_res.1
saveRDS(Caushi_T_NK_cells3, file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
rm(Caushi_T_NK_cells3)
Caushi_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Caushi_T_NK_cells4$RNA_snn_res.0.1 <- Caushi_T_NK_cells1$RNA_snn_res.0.1
Caushi_T_NK_cells4$RNA_snn_res.0.2 <- Caushi_T_NK_cells1$RNA_snn_res.0.2
Caushi_T_NK_cells4$RNA_snn_res.0.3 <- Caushi_T_NK_cells1$RNA_snn_res.0.3
Caushi_T_NK_cells4$RNA_snn_res.0.4 <- Caushi_T_NK_cells1$RNA_snn_res.0.4
Caushi_T_NK_cells4$RNA_snn_res.0.5 <- Caushi_T_NK_cells1$RNA_snn_res.0.5
Caushi_T_NK_cells4$RNA_snn_res.0.6 <- Caushi_T_NK_cells1$RNA_snn_res.0.6
Caushi_T_NK_cells4$RNA_snn_res.0.7 <- Caushi_T_NK_cells1$RNA_snn_res.0.7
Caushi_T_NK_cells4$RNA_snn_res.0.8 <- Caushi_T_NK_cells1$RNA_snn_res.0.8
Caushi_T_NK_cells4$RNA_snn_res.0.9 <- Caushi_T_NK_cells1$RNA_snn_res.0.9
Caushi_T_NK_cells4$RNA_snn_res.1 <- Caushi_T_NK_cells1$RNA_snn_res.1
saveRDS(Caushi_T_NK_cells4, file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"))
rm(Caushi_T_NK_cells4)
Caushi_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells.rds"))
Caushi_T_NK_cells5$RNA_snn_res.0.1 <- Caushi_T_NK_cells1$RNA_snn_res.0.1
Caushi_T_NK_cells5$RNA_snn_res.0.2 <- Caushi_T_NK_cells1$RNA_snn_res.0.2
Caushi_T_NK_cells5$RNA_snn_res.0.3 <- Caushi_T_NK_cells1$RNA_snn_res.0.3
Caushi_T_NK_cells5$RNA_snn_res.0.4 <- Caushi_T_NK_cells1$RNA_snn_res.0.4
Caushi_T_NK_cells5$RNA_snn_res.0.5 <- Caushi_T_NK_cells1$RNA_snn_res.0.5
Caushi_T_NK_cells5$RNA_snn_res.0.6 <- Caushi_T_NK_cells1$RNA_snn_res.0.6
Caushi_T_NK_cells5$RNA_snn_res.0.7 <- Caushi_T_NK_cells1$RNA_snn_res.0.7
Caushi_T_NK_cells5$RNA_snn_res.0.8 <- Caushi_T_NK_cells1$RNA_snn_res.0.8
Caushi_T_NK_cells5$RNA_snn_res.0.9 <- Caushi_T_NK_cells1$RNA_snn_res.0.9
Caushi_T_NK_cells5$RNA_snn_res.1 <- Caushi_T_NK_cells1$RNA_snn_res.1
saveRDS(Caushi_T_NK_cells5, file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells.rds"))
rm(Caushi_T_NK_cells5)
Caushi_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells.rds"))
Caushi_T_NK_cells6$RNA_snn_res.0.1 <- Caushi_T_NK_cells1$RNA_snn_res.0.1
Caushi_T_NK_cells6$RNA_snn_res.0.2 <- Caushi_T_NK_cells1$RNA_snn_res.0.2
Caushi_T_NK_cells6$RNA_snn_res.0.3 <- Caushi_T_NK_cells1$RNA_snn_res.0.3
Caushi_T_NK_cells6$RNA_snn_res.0.4 <- Caushi_T_NK_cells1$RNA_snn_res.0.4
Caushi_T_NK_cells6$RNA_snn_res.0.5 <- Caushi_T_NK_cells1$RNA_snn_res.0.5
Caushi_T_NK_cells6$RNA_snn_res.0.6 <- Caushi_T_NK_cells1$RNA_snn_res.0.6
Caushi_T_NK_cells6$RNA_snn_res.0.7 <- Caushi_T_NK_cells1$RNA_snn_res.0.7
Caushi_T_NK_cells6$RNA_snn_res.0.8 <- Caushi_T_NK_cells1$RNA_snn_res.0.8
Caushi_T_NK_cells6$RNA_snn_res.0.9 <- Caushi_T_NK_cells1$RNA_snn_res.0.9
Caushi_T_NK_cells6$RNA_snn_res.1 <- Caushi_T_NK_cells1$RNA_snn_res.1
saveRDS(Caushi_T_NK_cells6, file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells.rds"))
rm(Caushi_T_NK_cells6)


#######################################################
#### 2023.07.03 recalculate lisi
#### 2023。01.30 criteria on different resolution
########################################################

library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells_mapscore.rds"))
res_num <- c(0.1,0.15,0.16,0.2,0.217,0.22,00.25,0.28,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
TICA_ari <- NULL
TICA_lisi <- NULL
TICA_more50 <- NULL
for (i in 1:length(res_num)) {
  # TICA_ari <- c(TICA_ari,adjustedRandIndex(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Kathryn_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),data.frame(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # TICA_lisi <- c(TICA_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # TICA_more50 <- c(TICA_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
                                     # 2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  TICA_lisi <- c(TICA_lisi,mean(tmp[,1]))
  }

# Kathryn_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# panT_ari <- NULL
# panT_lisi <- NULL
# panT_more50 <- NULL
# for (i in 1:length(res_num)) {
#   # panT_ari <- c(panT_ari,adjustedRandIndex(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
#   # tmp <- compute_lisi(Embeddings(Kathryn_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),data.frame(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Kathryn_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4",paste0("RNA_snn_res.",res_num[i])))
#   # max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
#   # panT_lisi <- c(panT_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
#   panT_more50 <- c(panT_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
#                                          2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
#   
# }
# rm(Kathryn_T_NK_cells2)
Kathryn_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells_mapscore.rds"))
panblueprint_ari <- NULL
panblueprint_lisi <- NULL
panblueprint_more50 <- NULL
for (i in 1:length(res_num)) {
  # panblueprint_ari <- c(panblueprint_ari,adjustedRandIndex(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Kathryn_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),data.frame(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Kathryn_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # panblueprint_lisi <- c(panblueprint_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # panblueprint_more50 <- c(panblueprint_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  panblueprint_lisi <- c(panblueprint_lisi,mean(tmp[,1]))
}
rm(Kathryn_T_NK_cells3)
Kathryn_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells_mapscore.rds"))
ABTC_ari <- NULL
ABTC_lisi <- NULL
ABTC_more50 <- NULL
for (i in 1:length(res_num)) {
  # ABTC_ari <- c(ABTC_ari,adjustedRandIndex(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Kathryn_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),data.frame(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Kathryn_T_NK_cells4.predicted.to_ABTC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # ABTC_lisi <- c(ABTC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # ABTC_more50 <- c(ABTC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  ABTC_lisi <- c(ABTC_lisi,mean(tmp[,1]))
}
rm(Kathryn_T_NK_cells4)
Kathryn_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells_mapscore.rds"))
HCC_ari <- NULL
HCC_lisi <- NULL
HCC_more50 <- NULL
for (i in 1:length(res_num)) {
  # HCC_ari <- c(HCC_ari,adjustedRandIndex(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Kathryn_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),data.frame(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Kathryn_T_NK_cells5.predicted.to_HCC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # HCC_lisi <- c(HCC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # HCC_more50 <- c(HCC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  HCC_lisi <- c(HCC_lisi,mean(tmp[,1]))
}
rm(Kathryn_T_NK_cells5)
Kathryn_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells_mapscore.rds"))
STAD_ari <- NULL
STAD_lisi <- NULL
STAD_more50 <- NULL
for (i in 1:length(res_num)) {
  # STAD_ari <- c(STAD_ari,adjustedRandIndex(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Kathryn_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),data.frame(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Kathryn_T_NK_cells6.predicted.to_STAD_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # STAD_lisi <- c(STAD_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # STAD_more50 <- c(STAD_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  STAD_lisi <- c(STAD_lisi,mean(tmp[,1]))
}
rm(Kathryn_T_NK_cells6)
# data.frame(TICA_ari,panT_ari,panblueprint_ari,ABTC_ari,HCC_ari,STAD_ari)
# data.frame(TICA_lisi,panT_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
data.frame(TICA_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
# data.frame(TICA_more50,panT_more50,panblueprint_more50,ABTC_more50,HCC_more50,STAD_more50)


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_mapscore.rds"))
res_num <- c(0.1,0.15,0.2,0.21,0.22,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
TICA_ari <- NULL
TICA_lisi <- NULL
TICA_more50 <- NULL
for (i in 1:length(res_num)) {
  # TICA_ari <- c(TICA_ari,adjustedRandIndex(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Baolin_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),data.frame(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # TICA_lisi <- c(TICA_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # TICA_more50 <- c(TICA_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  # 2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  TICA_lisi <- c(TICA_lisi,mean(tmp[,1]))
}

# Baolin_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# panT_ari <- NULL
# panT_lisi <- NULL
# panT_more50 <- NULL
# for (i in 1:length(res_num)) {
#   # panT_ari <- c(panT_ari,adjustedRandIndex(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
#   # tmp <- compute_lisi(Embeddings(Baolin_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),data.frame(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Baolin_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4",paste0("RNA_snn_res.",res_num[i])))
#   # max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
#   # panT_lisi <- c(panT_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
#   panT_more50 <- c(panT_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
#                                          2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
#   
# }
# rm(Baolin_T_NK_cells2)
Baolin_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_mapscore.rds"))
panblueprint_ari <- NULL
panblueprint_lisi <- NULL
panblueprint_more50 <- NULL
for (i in 1:length(res_num)) {
  # panblueprint_ari <- c(panblueprint_ari,adjustedRandIndex(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Baolin_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),data.frame(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Baolin_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # panblueprint_lisi <- c(panblueprint_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # panblueprint_more50 <- c(panblueprint_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  panblueprint_lisi <- c(panblueprint_lisi,mean(tmp[,1]))
}
rm(Baolin_T_NK_cells3)
Baolin_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_mapscore.rds"))
ABTC_ari <- NULL
ABTC_lisi <- NULL
ABTC_more50 <- NULL
for (i in 1:length(res_num)) {
  # ABTC_ari <- c(ABTC_ari,adjustedRandIndex(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Baolin_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),data.frame(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Baolin_T_NK_cells4.predicted.to_ABTC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # ABTC_lisi <- c(ABTC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # ABTC_more50 <- c(ABTC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  ABTC_lisi <- c(ABTC_lisi,mean(tmp[,1]))
}
rm(Baolin_T_NK_cells4)
Baolin_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells_mapscore.rds"))
HCC_ari <- NULL
HCC_lisi <- NULL
HCC_more50 <- NULL
for (i in 1:length(res_num)) {
  # HCC_ari <- c(HCC_ari,adjustedRandIndex(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Baolin_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),data.frame(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Baolin_T_NK_cells5.predicted.to_HCC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # HCC_lisi <- c(HCC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # HCC_more50 <- c(HCC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  HCC_lisi <- c(HCC_lisi,mean(tmp[,1]))
}
rm(Baolin_T_NK_cells5)
Baolin_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells_mapscore.rds"))
STAD_ari <- NULL
STAD_lisi <- NULL
STAD_more50 <- NULL
for (i in 1:length(res_num)) {
  # STAD_ari <- c(STAD_ari,adjustedRandIndex(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Baolin_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),data.frame(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Baolin_T_NK_cells6.predicted.to_STAD_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # STAD_lisi <- c(STAD_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # STAD_more50 <- c(STAD_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  STAD_lisi <- c(STAD_lisi,mean(tmp[,1]))
}
rm(Baolin_T_NK_cells6)
# data.frame(TICA_ari,panT_ari,panblueprint_ari,ABTC_ari,HCC_ari,STAD_ari)
# data.frame(TICA_lisi,panT_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
data.frame(TICA_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
# data.frame(TICA_more50,panT_more50,panblueprint_more50,ABTC_more50,HCC_more50,STAD_more50)


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_mapscore.rds"))
res_num <- c(0.1,0.2,0.3,0.35,0.4,0.5,0.52,0.6,0.7,0.8,0.82,0.8245,0.9,1)
TICA_ari <- NULL
TICA_lisi <- NULL
TICA_more50 <- NULL
for (i in 1:length(res_num)) {
  # TICA_ari <- c(TICA_ari,adjustedRandIndex(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Moshe_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),data.frame(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # TICA_lisi <- c(TICA_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # TICA_more50 <- c(TICA_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  # 2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  TICA_lisi <- c(TICA_lisi,mean(tmp[,1]))
}

# Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# panT_ari <- NULL
# panT_lisi <- NULL
# panT_more50 <- NULL
# for (i in 1:length(res_num)) {
#   # panT_ari <- c(panT_ari,adjustedRandIndex(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
#   # tmp <- compute_lisi(Embeddings(Moshe_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),data.frame(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Moshe_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4",paste0("RNA_snn_res.",res_num[i])))
#   # max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
#   # panT_lisi <- c(panT_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
#   panT_more50 <- c(panT_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
#                                          2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
#   
# }
# rm(Moshe_T_NK_cells2)
Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_mapscore.rds"))
panblueprint_ari <- NULL
panblueprint_lisi <- NULL
panblueprint_more50 <- NULL
for (i in 1:length(res_num)) {
  # panblueprint_ari <- c(panblueprint_ari,adjustedRandIndex(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Moshe_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),data.frame(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Moshe_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # panblueprint_lisi <- c(panblueprint_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # panblueprint_more50 <- c(panblueprint_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  panblueprint_lisi <- c(panblueprint_lisi,mean(tmp[,1]))
}
rm(Moshe_T_NK_cells3)
Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_mapscore.rds"))
ABTC_ari <- NULL
ABTC_lisi <- NULL
ABTC_more50 <- NULL
for (i in 1:length(res_num)) {
  # ABTC_ari <- c(ABTC_ari,adjustedRandIndex(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Moshe_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),data.frame(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Moshe_T_NK_cells4.predicted.to_ABTC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # ABTC_lisi <- c(ABTC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # ABTC_more50 <- c(ABTC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  ABTC_lisi <- c(ABTC_lisi,mean(tmp[,1]))
}
rm(Moshe_T_NK_cells4)
Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells_mapscore.rds"))
HCC_ari <- NULL
HCC_lisi <- NULL
HCC_more50 <- NULL
for (i in 1:length(res_num)) {
  # HCC_ari <- c(HCC_ari,adjustedRandIndex(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Moshe_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),data.frame(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Moshe_T_NK_cells5.predicted.to_HCC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # HCC_lisi <- c(HCC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # HCC_more50 <- c(HCC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  HCC_lisi <- c(HCC_lisi,mean(tmp[,1]))
}
rm(Moshe_T_NK_cells5)
Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells_mapscore.rds"))
STAD_ari <- NULL
STAD_lisi <- NULL
STAD_more50 <- NULL
for (i in 1:length(res_num)) {
  # STAD_ari <- c(STAD_ari,adjustedRandIndex(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Moshe_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),data.frame(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Moshe_T_NK_cells6.predicted.to_STAD_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # STAD_lisi <- c(STAD_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # STAD_more50 <- c(STAD_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  STAD_lisi <- c(STAD_lisi,mean(tmp[,1]))
}
rm(Moshe_T_NK_cells6)
# data.frame(TICA_ari,panT_ari,panblueprint_ari,ABTC_ari,HCC_ari,STAD_ari)
# data.frame(TICA_lisi,panT_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
data.frame(TICA_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
# data.frame(TICA_more50,panT_more50,panblueprint_more50,ABTC_more50,HCC_more50,STAD_more50)



library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_mapscore.rds"))
res_num <- c(0.1,0.2,0.208,0.21,0.3,0.31,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1)
TICA_ari <- NULL
TICA_lisi <- NULL
TICA_more50 <- NULL
for (i in 1:length(res_num)) {
  # TICA_ari <- c(TICA_ari,adjustedRandIndex(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Caushi_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),data.frame(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # TICA_lisi <- c(TICA_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # TICA_more50 <- c(TICA_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  # 2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  TICA_lisi <- c(TICA_lisi,mean(tmp[,1]))
}

# Caushi_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# panT_ari <- NULL
# panT_lisi <- NULL
# panT_more50 <- NULL
# for (i in 1:length(res_num)) {
#   # panT_ari <- c(panT_ari,adjustedRandIndex(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
#   # tmp <- compute_lisi(Embeddings(Caushi_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),data.frame(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Caushi_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4",paste0("RNA_snn_res.",res_num[i])))
#   # max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
#   # panT_lisi <- c(panT_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
#   panT_more50 <- c(panT_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
#                                          2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
#   
# }
# rm(Caushi_T_NK_cells2)
Caushi_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_mapscore.rds"))
panblueprint_ari <- NULL
panblueprint_lisi <- NULL
panblueprint_more50 <- NULL
for (i in 1:length(res_num)) {
  # panblueprint_ari <- c(panblueprint_ari,adjustedRandIndex(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Caushi_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),data.frame(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Caushi_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # panblueprint_lisi <- c(panblueprint_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # panblueprint_more50 <- c(panblueprint_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  panblueprint_lisi <- c(panblueprint_lisi,mean(tmp[,1]))
}
rm(Caushi_T_NK_cells3)
Caushi_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_mapscore.rds"))
ABTC_ari <- NULL
ABTC_lisi <- NULL
ABTC_more50 <- NULL
for (i in 1:length(res_num)) {
  # ABTC_ari <- c(ABTC_ari,adjustedRandIndex(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Caushi_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),data.frame(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Caushi_T_NK_cells4.predicted.to_ABTC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # ABTC_lisi <- c(ABTC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # ABTC_more50 <- c(ABTC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  ABTC_lisi <- c(ABTC_lisi,mean(tmp[,1]))
}
rm(Caushi_T_NK_cells4)
Caushi_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells_mapscore.rds"))
HCC_ari <- NULL
HCC_lisi <- NULL
HCC_more50 <- NULL
for (i in 1:length(res_num)) {
  # HCC_ari <- c(HCC_ari,adjustedRandIndex(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Caushi_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),data.frame(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Caushi_T_NK_cells5.predicted.to_HCC_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # HCC_lisi <- c(HCC_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # HCC_more50 <- c(HCC_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  HCC_lisi <- c(HCC_lisi,mean(tmp[,1]))
}
rm(Caushi_T_NK_cells5)
Caushi_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells_mapscore.rds"))
STAD_ari <- NULL
STAD_lisi <- NULL
STAD_more50 <- NULL
for (i in 1:length(res_num)) {
  # STAD_ari <- c(STAD_ari,adjustedRandIndex(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]][,1]))
  tmp <- compute_lisi(Embeddings(Caushi_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),data.frame(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]),c("Caushi_T_NK_cells6.predicted.to_STAD_T_NK_cells",paste0("RNA_snn_res.",res_num[i])))
  max_tmp <- max(abs(tmp[,1]-tmp[,2]));min_tmp <- min(abs(tmp[,1]-tmp[,2]))
  # STAD_lisi <- c(STAD_lisi,mean((abs(tmp[,1]-tmp[,2])-min_tmp)/(max_tmp-min_tmp)))
  # STAD_more50 <- c(STAD_more50,sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]])))),2,function(x) x/sum(x)),
  #                                        2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells1[[paste0("RNA_snn_res.",res_num[i])]]))))))
  
  # max_tmp <- max(tmp[,1]); min_tmp <- min(tmp[,2])
  tmp[,1] <- (max_tmp-abs(tmp[,1]-tmp[,2]))/(max_tmp-min_tmp)
  STAD_lisi <- c(STAD_lisi,mean(tmp[,1]))
}
rm(Caushi_T_NK_cells6)
# data.frame(TICA_ari,panT_ari,panblueprint_ari,ABTC_ari,HCC_ari,STAD_ari)
# data.frame(TICA_lisi,panT_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
data.frame(TICA_lisi,panblueprint_lisi,ABTC_lisi,HCC_lisi,STAD_lisi)
# data.frame(TICA_more50,panT_more50,panblueprint_more50,ABTC_more50,HCC_more50,STAD_more50)


#######################################################
#### 2023.02.02 re-caculate ASW
#####################################################


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# workdir <- "/scratch/beauchamp_lab/JingYang/tumor_immune_cell_atlas_estimation/data/"
# TICA_T_NK_cells <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells.rds"))
# TICA_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = TICA_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(TICA_T_NK_cells[["pca"]])))[,3])
# print(TICA_ASW)
# rm(TICA_T_NK_cells)
# pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_simplized.rds"))
# panT_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = pancancer_T_cells_CD8_CD4$cell_types_from_author_simplized)), dist = dist(Embeddings(pancancer_T_cells_CD8_CD4[["pca"]])))[,3])
# print(panT_ASW)
# rm(pancancer_T_cells_CD8_CD4)
# pancancer_blueprint_T_NK_cells <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells.rds"))
# pancancer_blueprint_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = pancancer_blueprint_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(pancancer_blueprint_T_NK_cells[["pca"]])))[,3])
# print(pancancer_blueprint_ASW)
# rm(pancancer_blueprint_T_NK_cells)
# ABTC_T_NK_cells <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells.rds"))
# ABTC_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = ABTC_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(ABTC_T_NK_cells[["pca"]])))[,3])
# print(ABTC_ASW)
# rm(ABTC_T_NK_cells)
# HCC_T_NK_cells <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells.rds"))
# HCC_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = HCC_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(HCC_T_NK_cells[["pca"]])))[,3])
# print(HCC_ASW)
# rm(HCC_T_NK_cells)
# STAD_T_NK_cells <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells.rds"))
# STAD_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = STAD_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(STAD_T_NK_cells[["pca"]])))[,3])
# print(STAD_ASW)
# rm(STAD_T_NK_cells)


Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
Baolin_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = Baolin_T_NK_cells1$RNA_snn_res.0.21)), dist = dist(Embeddings(Baolin_T_NK_cells1[["pca"]])))[,3])
print(Baolin_ASW)
rm(Baolin_T_NK_cells1)

Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
Moshe_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = Moshe_T_NK_cells1$RNA_snn_res.0.35)), dist = dist(Embeddings(Moshe_T_NK_cells1[["pca"]])))[,3])
print(Moshe_ASW)
rm(Moshe_T_NK_cells1)

Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
Caushi_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = Caushi_T_NK_cells1$RNA_snn_res.0.3)), dist = dist(Embeddings(Caushi_T_NK_cells1[["pca"]])))[,3])
print(Caushi_ASW)
rm(Caushi_T_NK_cells1)

TICA <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells.rds"))
TICA <- subset(TICA, downsample = 3000)
TICA_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = TICA$cell_types_from_author)), dist = dist(Embeddings(TICA[["pca"]])))[,3])
print(TICA_ASW)
rm(TICA)

pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_simplized.rds"))
pancancer_T_cells_CD8_CD4_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = pancancer_T_cells_CD8_CD4$cell_types_from_author)), dist = dist(Embeddings(pancancer_T_cells_CD8_CD4[["pca"]])))[,3])
print(pancancer_T_cells_CD8_CD4_ASW)
rm(pancancer_T_cells_CD8_CD4)

pancancer_blueprint_T_NK_cells <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells.rds"))
pancancer_blueprint_T_NK_cells_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = pancancer_blueprint_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(pancancer_blueprint_T_NK_cells[["pca"]])))[,3])
print(pancancer_blueprint_T_NK_cells_ASW)
rm(pancancer_blueprint_T_NK_cells)

ABTC_T_NK_cells <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells.rds"))
ABTC_T_NK_cells_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = ABTC_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(ABTC_T_NK_cells[["pca"]])))[,3])
print(ABTC_T_NK_cells_ASW)
rm(ABTC_T_NK_cells)

HCC_T_NK_cells <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells.rds"))
HCC_T_NK_cells_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = HCC_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(HCC_T_NK_cells[["pca"]])))[,3])
print(HCC_T_NK_cells_ASW)
rm(HCC_T_NK_cells)

STAD_T_NK_cells <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells.rds"))
STAD_T_NK_cells_ASW <- mean(silhouette(x = as.numeric(x = as.factor(x = STAD_T_NK_cells$cell_types_from_author)), dist = dist(Embeddings(STAD_T_NK_cells[["pca"]])))[,3])
print(STAD_T_NK_cells_ASW)
rm(STAD_T_NK_cells)
############################################################
#### 2022.11.09 criteria for query datasets
############################################################

library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# workdir <- "/scratch/beauchamp_lab/JingYang/tumor_immune_cell_atlas_estimation/data/"

Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
Kathryn_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Kathryn_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Kathryn_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Kathryn_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells.rds"))
Kathryn_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells.rds"))
#ARI
ARI_output <- c(adjustedRandIndex(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1$cellcluster_ann),
                adjustedRandIndex(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells2$cellcluster_ann),
                adjustedRandIndex(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells3$cellcluster_ann),
                adjustedRandIndex(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells4$cellcluster_ann),
                adjustedRandIndex(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells5$cellcluster_ann),
                adjustedRandIndex(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells6$cellcluster_ann))
print("Kathryn ARI: ")
print(ARI_output)
# lisi
lisi_output <- c(mean(compute_lisi(
  Embeddings(Kathryn_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),
  data.frame(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1$cellcluster_ann),
  c("Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Kathryn_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),
    data.frame(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells2$cellcluster_ann),
    c("Kathryn_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4"))[,1]),
  mean(compute_lisi(
    Embeddings(Kathryn_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),
    data.frame(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells3$cellcluster_ann),
    c("Kathryn_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Kathryn_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),
    data.frame(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells4$cellcluster_ann),
    c("Kathryn_T_NK_cells4.predicted.to_ABTC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Kathryn_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),
    data.frame(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells5$cellcluster_ann),
    c("Kathryn_T_NK_cells5.predicted.to_HCC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Kathryn_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),
    data.frame(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells6$cellcluster_ann),
    c("Kathryn_T_NK_cells6.predicted.to_STAD_T_NK_cells"))[,1])
)
print("Kathryn lisi: "); print(lisi_output)
# ASW
ASW_output <- c(mean(silhouette(x = as.numeric(x = as.factor(x = Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells)), dist = dist(Embeddings(Kathryn_T_NK_cells1[["to_TICA_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4)), dist = dist(Embeddings(Kathryn_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells)), dist = dist(Embeddings(Kathryn_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells)), dist = dist(Embeddings(Kathryn_T_NK_cells4[["to_ABTC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells)), dist = dist(Embeddings(Kathryn_T_NK_cells5[["to_HCC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells)), dist = dist(Embeddings(Kathryn_T_NK_cells6[["to_STAD_T_NK_cells_pca"]])))[,3])
)
print("Kathryn ASW: "); print(ASW_output)
# the probability of each cluster whose >50% cells can be annotated in one cell type of reference
annotatable_output <- c(sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,Kathryn_T_NK_cells1$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells2$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Kathryn_T_NK_cells2$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells3$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Kathryn_T_NK_cells3$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells4$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Kathryn_T_NK_cells4$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells5$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells,Kathryn_T_NK_cells5$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells6$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells,Kathryn_T_NK_cells6$cellcluster_ann)))))
)
print("Kathryn annotatable: "); print(annotatable_output)
rm(list=ls())


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
workdir <- "/scratch/beauchamp_lab/JingYang/tumor_immune_cell_atlas_estimation/data/"
Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
Baolin_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Baolin_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Baolin_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"));Baolin_T_NK_cells4 <- subset(Baolin_T_NK_cells4,downsample=3000)
Baolin_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"));Baolin_T_NK_cells5 <- subset(Baolin_T_NK_cells5,downsample=3000)
Baolin_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"));Baolin_T_NK_cells6 <- subset(Baolin_T_NK_cells6,downsample=3000)
#ARI
ARI_output <- c(adjustedRandIndex(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1$cellcluster_ann),
                adjustedRandIndex(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells2$cellcluster_ann),
                adjustedRandIndex(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells3$cellcluster_ann),
                adjustedRandIndex(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells4$cellcluster_ann),
                adjustedRandIndex(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells5$cellcluster_ann),
                adjustedRandIndex(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells6$cellcluster_ann))
print("Baolin ARI: ");print(ARI_output)
# lisi
lisi_output <- c(mean(compute_lisi(
  Embeddings(Baolin_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),
  data.frame(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1$cellcluster_ann),
  c("Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Baolin_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),
    data.frame(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells2$cellcluster_ann),
    c("Baolin_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4"))[,1]),
  mean(compute_lisi(
    Embeddings(Baolin_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),
    data.frame(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells3$cellcluster_ann),
    c("Baolin_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Baolin_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),
    data.frame(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells4$cellcluster_ann),
    c("Baolin_T_NK_cells4.predicted.to_ABTC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Baolin_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),
    data.frame(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells5$cellcluster_ann),
    c("Baolin_T_NK_cells5.predicted.to_HCC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Baolin_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),
    data.frame(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells6$cellcluster_ann),
    c("Baolin_T_NK_cells6.predicted.to_STAD_T_NK_cells"))[,1])
)
print("Baolin lisi: "); print(lisi_output)
ASW
ASW_output <- c(mean(silhouette(x = as.numeric(x = as.factor(x = Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells)), dist = dist(Embeddings(Baolin_T_NK_cells1[["to_TICA_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4)), dist = dist(Embeddings(Baolin_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells)), dist = dist(Embeddings(Baolin_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells)), dist = dist(Embeddings(Baolin_T_NK_cells4[["to_ABTC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells)), dist = dist(Embeddings(Baolin_T_NK_cells5[["to_HCC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells)), dist = dist(Embeddings(Baolin_T_NK_cells6[["to_STAD_T_NK_cells_pca"]])))[,3])
)
# print("Baolin ASW: ");print(ASW_output)
# the probability of each cluster whose >50% cells can be annotated in one cell type of reference
annotatable_output <- c(sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,Baolin_T_NK_cells1$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells2$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Baolin_T_NK_cells2$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells3$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Baolin_T_NK_cells3$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells4$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Baolin_T_NK_cells4$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells5$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells,Baolin_T_NK_cells5$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells6$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells,Baolin_T_NK_cells6$cellcluster_ann)))))
)
print("Baolin annotatable: ");print(annotatable_output)
rm(list=ls())

# 
library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
workdir <- "/scratch/beauchamp_lab/JingYang/tumor_immune_cell_atlas_estimation/data/"
Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
# #ARI
ARI_output <- c(adjustedRandIndex(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1$cellcluster_ann),
                adjustedRandIndex(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells2$cellcluster_ann),
                adjustedRandIndex(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells3$cellcluster_ann),
                adjustedRandIndex(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells4$cellcluster_ann),
                adjustedRandIndex(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells5$cellcluster_ann),
                adjustedRandIndex(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells6$cellcluster_ann))
print("Moshe ARI: "); print(ARI_output)
# lisi
lisi_output <- c(mean(compute_lisi(
  Embeddings(Moshe_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),
  data.frame(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1$cellcluster_ann),
  c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Moshe_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),
    data.frame(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells2$cellcluster_ann),
    c("Moshe_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4"))[,1]),
  mean(compute_lisi(
    Embeddings(Moshe_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),
    data.frame(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells3$cellcluster_ann),
    c("Moshe_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Moshe_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),
    data.frame(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells4$cellcluster_ann),
    c("Moshe_T_NK_cells4.predicted.to_ABTC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Moshe_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),
    data.frame(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells5$cellcluster_ann),
    c("Moshe_T_NK_cells5.predicted.to_HCC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Moshe_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),
    data.frame(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells6$cellcluster_ann),
    c("Moshe_T_NK_cells6.predicted.to_STAD_T_NK_cells"))[,1])
)
print("Moshe lisi: "); print(lisi_output)
# ASW
ASW_output <- c(mean(silhouette(x = as.numeric(x = as.factor(x = Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells)), dist = dist(Embeddings(Moshe_T_NK_cells1[["to_TICA_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4)), dist = dist(Embeddings(Moshe_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells)), dist = dist(Embeddings(Moshe_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells)), dist = dist(Embeddings(Moshe_T_NK_cells4[["to_ABTC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells)), dist = dist(Embeddings(Moshe_T_NK_cells5[["to_HCC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells)), dist = dist(Embeddings(Moshe_T_NK_cells6[["to_STAD_T_NK_cells_pca"]])))[,3])
)
print("Moshe ASW: "); print(ASW_output)
# the probability of each cluster whose >50% cells can be annotated in one cell type of reference
annotatable_output <- c(sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,Moshe_T_NK_cells1$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells2$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Moshe_T_NK_cells2$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells3$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Moshe_T_NK_cells3$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells4$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Moshe_T_NK_cells4$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells5$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells,Moshe_T_NK_cells5$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells6$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells,Moshe_T_NK_cells6$cellcluster_ann)))))
)
print("Moshe annotatable: "); print(annotatable_output)
rm(list=ls())


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
workdir <- "/scratch/beauchamp_lab/JingYang/tumor_immune_cell_atlas_estimation/data/"
Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"));Caushi_T_NK_cells1 <- subset(Caushi_T_NK_cells1,downsample=3000)
Caushi_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"));Caushi_T_NK_cells2 <- subset(Caushi_T_NK_cells2,downsample=3000)
Caushi_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"));Caushi_T_NK_cells3 <- subset(Caushi_T_NK_cells3,downsample=3000)
Caushi_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"));Caushi_T_NK_cells4 <- subset(Caushi_T_NK_cells4,downsample=3000)
Caushi_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells.rds"));Caushi_T_NK_cells5 <- subset(Caushi_T_NK_cells5,downsample=3000)
Caushi_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells.rds"));Caushi_T_NK_cells6 <- subset(Caushi_T_NK_cells6,downsample=3000)
#ARI
ARI_output <- c(adjustedRandIndex(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1$cellcluster_ann),
                adjustedRandIndex(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells2$cellcluster_ann),
                adjustedRandIndex(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells3$cellcluster_ann),
                adjustedRandIndex(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells4$cellcluster_ann),
                adjustedRandIndex(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells5$cellcluster_ann),
                adjustedRandIndex(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells6$cellcluster_ann))
print("Caush ARI: "); print(ARI_output)
# lisi
lisi_output <- c(mean(compute_lisi(
  Embeddings(Caushi_T_NK_cells1[["to_TICA_T_NK_cells_umap"]]),
  data.frame(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1$cellcluster_ann),
  c("Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Caushi_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_umap"]]),
    data.frame(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells2$cellcluster_ann),
    c("Caushi_T_NK_cells2.predicted.to_pancancer_T_cells_CD8_CD4"))[,1]),
  mean(compute_lisi(
    Embeddings(Caushi_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_umap"]]),
    data.frame(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells3$cellcluster_ann),
    c("Caushi_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Caushi_T_NK_cells4[["to_ABTC_T_NK_cells_umap"]]),
    data.frame(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells4$cellcluster_ann),
    c("Caushi_T_NK_cells4.predicted.to_ABTC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Caushi_T_NK_cells5[["to_HCC_T_NK_cells_umap"]]),
    data.frame(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells5$cellcluster_ann),
    c("Caushi_T_NK_cells5.predicted.to_HCC_T_NK_cells"))[,1]),
  mean(compute_lisi(
    Embeddings(Caushi_T_NK_cells6[["to_STAD_T_NK_cells_umap"]]),
    data.frame(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells6$cellcluster_ann),
    c("Caushi_T_NK_cells6.predicted.to_STAD_T_NK_cells"))[,1])
)
print("Caushi lisi: "); print(lisi_output)
# ASW
ASW_output <- c(mean(silhouette(x = as.numeric(x = as.factor(x = Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells)), dist = dist(Embeddings(Caushi_T_NK_cells1[["to_TICA_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4)), dist = dist(Embeddings(Caushi_T_NK_cells2[["to_pancancer_T_cells_CD8_CD4_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells)), dist = dist(Embeddings(Caushi_T_NK_cells3[["to_pancancer_blueprint_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells)), dist = dist(Embeddings(Caushi_T_NK_cells4[["to_ABTC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells)), dist = dist(Embeddings(Caushi_T_NK_cells5[["to_HCC_T_NK_cells_pca"]])))[,3]),
                mean(silhouette(x = as.numeric(x = as.factor(x = Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells)), dist = dist(Embeddings(Caushi_T_NK_cells6[["to_STAD_T_NK_cells_pca"]])))[,3])
)
print("Caush ASW: "); print(ASW)
# the probability of each cluster whose >50% cells can be annotated in one cell type of reference
annotatable_output <- c(sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,Caushi_T_NK_cells1$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells2$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4,Caushi_T_NK_cells2$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells3$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells,Caushi_T_NK_cells3$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells4$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells,Caushi_T_NK_cells4$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells5$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells,Caushi_T_NK_cells5$cellcluster_ann))))),
                        sum(apply(apply(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells6$cellcluster_ann)))),2,function(x) x/sum(x)),
                                  2,function(x) sum(na.omit(x)>0.5)>0))/ncol(as.matrix(unlist(table(data.frame(Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells,Caushi_T_NK_cells6$cellcluster_ann)))))
)
print("Caush annotatable: "); print(annotatable_output)

############################################
#### 2022.11.22 ICB response
############################################
library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# workdir <- "/scratch/beauchamp_lab/JingYang/tumor_immune_cell_atlas_estimation/data/"
SampleInfo_Kathryn <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Clonal replacement of tumor-specific T cells following PD-1 blockade/clinical_meta.txt",header=T,sep="\t")

Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
Kathryn_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Kathryn_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Kathryn_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Kathryn_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells.rds"))
Kathryn_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells.rds"))
Kathryn_T_NK_cells <- Kathryn_T_NK_cells6

count.table<-table(data.frame(Kathryn_T_NK_cells$SampleId,Kathryn_T_NK_cells$predicted.to_STAD_T_NK_cells))
# count.table<-table(data.frame(Kathryn_T_NK_cells$SampleId,Kathryn_T_NK_cells$cellcluster_ann))
# count.table<-table(data.frame(Kathryn_T_NK_cells$SampleId,Kathryn_T_NK_cells$RNA_snn_res.0.16))
prop_celltype<-count.table/apply(count.table,1,sum)
prop_celltype<-as.matrix(unclass(prop_celltype))
prop_celltype <- merge(SampleInfo_Kathryn[,c("SampleId","condition")],prop_celltype,by.x="SampleId",by.y="row.names",all.y=T)
prop_celltype <- rbind(prop_celltype,c("p.value_wilcox_test","NA",apply(prop_celltype[,3:ncol(prop_celltype)],2,function(x) wilcox.test(x[prop_celltype[,2]=="R" & !is.na(prop_celltype[,2])],x[prop_celltype[,2]=="NR" & !is.na(prop_celltype[,2])])$p.value)))



SampleInfo_Baolin <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Temporal single-cell tracing reveals clonal revival and expansion of precursor exhausted T cells during anti-PD-1 therapy in lung cancer/GSE179994_Tcell.metadata_with_annotation.txt",header=T,sep="\t");dim(SampleInfo_Baolin)
Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
Baolin_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Baolin_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Baolin_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Baolin_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"))
Baolin_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"))
Baolin_T_NK_cells <- Baolin_T_NK_cells6

count.table<-table(data.frame(Baolin_T_NK_cells$SampleId,Baolin_T_NK_cells$predicted.to_STAD_T_NK_cells))
# count.table<-table(data.frame(Baolin_T_NK_cells$SampleId,Baolin_T_NK_cells$cellcluster_ann))
prop_celltype<-count.table/apply(count.table,1,sum)
prop_celltype<-as.matrix(unclass(prop_celltype))
prop_celltype <- merge(SampleInfo_Baolin[,c("SampleId","condition")],prop_celltype,by.x="SampleId",by.y="row.names",all.y=T)
prop_celltype <- rbind(prop_celltype,c("p.value_wilcox_test","NA",apply(prop_celltype[,3:ncol(prop_celltype)],2,function(x) wilcox.test(x[prop_celltype[,2]=="R" & !is.na(prop_celltype[,2])],x[prop_celltype[,2]=="NR" & !is.na(prop_celltype[,2])])$p.value)))


SampleInfo_Moshe <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From Board/Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma/cells_all_scp.txt",header=T,sep="\t")
SampleInfo_Moshe <- SampleInfo_Moshe[-1,]
SampleInfo_Moshe <- unique(SampleInfo_Moshe[,c(2:3,8,10:13)])
SampleInfo_Moshe <- SampleInfo_Moshe[SampleInfo_Moshe$response %in% c("R","NR"),]
SampleInfo_Moshe <- data.frame(SampleId = SampleInfo_Moshe$patient,condition=SampleInfo_Moshe$response,treatment=SampleInfo_Moshe$therapy,prepost=SampleInfo_Moshe$prepost,gender=SampleInfo_Moshe$gender,age=SampleInfo_Moshe$age,os=SampleInfo_Moshe$survival_days,stringsAsFactors=F)
SampleInfo_Moshe$prepost <- tolower(SampleInfo_Moshe$prepost)
SampleInfo_Moshe <- SampleInfo_Moshe[c(1:10,12:48),]
Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
Moshe_T_NK_cells <- Moshe_T_NK_cells6

count.table<-table(data.frame(Moshe_T_NK_cells$SampleId,Moshe_T_NK_cells$predicted.to_STAD_T_NK_cells))
# count.table<-table(data.frame(Moshe_T_NK_cells$SampleId,Moshe_T_NK_cells$cellcluster_ann))
prop_celltype<-count.table/apply(count.table,1,sum)
prop_celltype<-as.matrix(unclass(prop_celltype))
prop_celltype <- merge(SampleInfo_Moshe[,c("SampleId","condition")],prop_celltype,by.x="SampleId",by.y="row.names",all.y=T)
prop_celltype <- rbind(prop_celltype,c("p.value_wilcox_test","NA",apply(prop_celltype[,3:ncol(prop_celltype)],2,function(x) wilcox.test(x[prop_celltype[,2]=="R" & !is.na(prop_celltype[,2])],x[prop_celltype[,2]=="NR" & !is.na(prop_celltype[,2])])$p.value)))



SampleInfo_Caushi <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Transcriptional programs of neoantigen-specific TIL in anti-PD-1-treated lung cancers/clinical_meta4_simplified.txt",header=T,sep="\t")
Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
Caushi_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Caushi_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Caushi_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Caushi_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells.rds"))
Caushi_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells.rds"))
Caushi_T_NK_cells <- Caushi_T_NK_cells6

count.table<-table(data.frame(Caushi_T_NK_cells$SampleId,Caushi_T_NK_cells$predicted.to_STAD_T_NK_cells))
# count.table<-table(data.frame(Caushi_T_NK_cells$SampleId,Caushi_T_NK_cells$cellcluster_ann))
prop_celltype<-count.table/apply(count.table,1,sum)
prop_celltype<-as.matrix(unclass(prop_celltype))
prop_celltype <- merge(SampleInfo_Caushi[,c("SampleId","condition")],prop_celltype,by.x="SampleId",by.y="row.names",all.y=T)
prop_celltype <- rbind(prop_celltype,c("p.value_wilcox_test","NA",apply(prop_celltype[,3:ncol(prop_celltype)],2,function(x) wilcox.test(x[prop_celltype[,2]=="MPR" & !is.na(prop_celltype[,2])],x[prop_celltype[,2]=="non-MPR" & !is.na(prop_celltype[,2])])$p.value)))


########################################
## 2023.01.11 re-detection rate of marker genes
#########################################
library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
library(plyr)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells
# Kathryn_T_NK_cells1.markers <- FindAllMarkers(Kathryn_T_NK_cells1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Kathryn_T_NK_cells1.markers,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
Kathryn_T_NK_cells1.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
TICA_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells_markers.rds"))
output_ove1 <- NULL
types_got <- as.character(unique(Kathryn_T_NK_cells1.markers$cluster))
  for (i in 1:length(types_got)) {
    X <- unique(Kathryn_T_NK_cells1.markers[Kathryn_T_NK_cells1.markers$cluster==types_got[i],"gene"])
    Y <- unique(TICA_T_NK_cells_markers[TICA_T_NK_cells_markers$cluster==types_got[i],"gene"])
    X_Y <- intersect(X,Y)
    output_ove1 <- c(output_ove1,2*length(X_Y)/(length(X)+length(Y)))
  }
names(output_ove1) <- paste0("TICA: ",types_got)
rm("Kathryn_T_NK_cells1")

# Kathryn_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Kathryn_T_NK_cells2) <- Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
# Kathryn_T_NK_cells2.markers <- FindAllMarkers(Kathryn_T_NK_cells2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Kathryn_T_NK_cells2.markers,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
Kathryn_T_NK_cells2.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
pancancer_T_cells_CD8_CD4_markers <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_markers_simplized.rds"))
output_ove2 <- NULL
types_got <- as.character(unique(Kathryn_T_NK_cells2.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Kathryn_T_NK_cells2.markers[Kathryn_T_NK_cells2.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_T_cells_CD8_CD4_markers[pancancer_T_cells_CD8_CD4_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove2 <- c(output_ove2,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove2) <- types_got
rm("Kathryn_T_NK_cells2")

# Kathryn_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells3) <- Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
# Kathryn_T_NK_cells3.markers <- FindAllMarkers(Kathryn_T_NK_cells3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Kathryn_T_NK_cells3.markers,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
Kathryn_T_NK_cells3.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
pancancer_blueprint_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells_markers.rds"))
output_ove3 <- NULL
types_got <- as.character(unique(Kathryn_T_NK_cells3.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Kathryn_T_NK_cells3.markers[Kathryn_T_NK_cells3.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_blueprint_T_NK_cells_markers[pancancer_blueprint_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove3 <- c(output_ove3,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove3) <- paste0("pan- blueprint: ",types_got)
rm("Kathryn_T_NK_cells3")


# Kathryn_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells4) <- Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells
# Kathryn_T_NK_cells4.markers <- FindAllMarkers(Kathryn_T_NK_cells4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Kathryn_T_NK_cells4.markers,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
Kathryn_T_NK_cells4.markers <- paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells_markers.rds")
ABTC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells_markers.rds"))
output_ove4 <- NULL
types_got <- as.character(unique(Kathryn_T_NK_cells4.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Kathryn_T_NK_cells4.markers[Kathryn_T_NK_cells4.markers$cluster==types_got[i],"gene"])
  Y <- unique(ABTC_T_NK_cells_markers[ABTC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove4 <- c(output_ove4,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove4) <- paste0("ABTC: ",types_got)
rm("Kathryn_T_NK_cells4")


# Kathryn_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells5) <- Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells
# Kathryn_T_NK_cells5.markers <- FindAllMarkers(Kathryn_T_NK_cells5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Kathryn_T_NK_cells5.markers,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
Kathryn_T_NK_cells5.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
HCC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells_markers.rds"))
output_ove5 <- NULL
types_got <- as.character(unique(Kathryn_T_NK_cells5.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Kathryn_T_NK_cells5.markers[Kathryn_T_NK_cells5.markers$cluster==types_got[i],"gene"])
  Y <- unique(HCC_T_NK_cells_markers[HCC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove5 <- c(output_ove5,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove5) <- paste0("HCC: ",types_got)
rm("Kathryn_T_NK_cells5")


# Kathryn_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells6) <- Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells
# Kathryn_T_NK_cells6.markers <- FindAllMarkers(Kathryn_T_NK_cells6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Kathryn_T_NK_cells6.markers,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
Kathryn_T_NK_cells6.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
STAD_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells_markers.rds"))
output_ove6 <- NULL
types_got <- as.character(unique(Kathryn_T_NK_cells6.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Kathryn_T_NK_cells6.markers[Kathryn_T_NK_cells6.markers$cluster==types_got[i],"gene"])
  Y <- unique(STAD_T_NK_cells_markers[STAD_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove6 <- c(output_ove6,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove6) <- paste0("STAD: ",types_got)
rm("Kathryn_T_NK_cells6")
# c(mean(output_ove1),mean(output_ove2),mean(output_ove3),
#   mean(output_ove4),mean(output_ove5),mean(output_ove6))
Kathryn_redetection <- c(output_ove1,output_ove2,output_ove3,output_ove4,output_ove5,output_ove6)

library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells
# Baolin_T_NK_cells1.markers <- FindAllMarkers(Baolin_T_NK_cells1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Baolin_T_NK_cells1.markers,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
Baolin_T_NK_cells1.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
TICA_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells_markers.rds"))
output_ove1 <- NULL
types_got <- as.character(unique(Baolin_T_NK_cells1.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Baolin_T_NK_cells1.markers[Baolin_T_NK_cells1.markers$cluster==types_got[i],"gene"])
  Y <- unique(TICA_T_NK_cells_markers[TICA_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove1 <- c(output_ove1,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove1) <- paste0("TICA: ",types_got)
rm("Baolin_T_NK_cells1")

# Baolin_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Baolin_T_NK_cells2) <- Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
# Baolin_T_NK_cells2.markers <- FindAllMarkers(Baolin_T_NK_cells2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Baolin_T_NK_cells2.markers,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
Baolin_T_NK_cells2.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
pancancer_T_cells_CD8_CD4_markers <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_markers_simplized.rds"))
output_ove2 <- NULL
types_got <- as.character(unique(Baolin_T_NK_cells2.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Baolin_T_NK_cells2.markers[Baolin_T_NK_cells2.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_T_cells_CD8_CD4_markers[pancancer_T_cells_CD8_CD4_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove2 <- c(output_ove2,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove2) <- types_got
rm("Baolin_T_NK_cells2")

# Baolin_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells3) <- Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
# Baolin_T_NK_cells3.markers <- FindAllMarkers(Baolin_T_NK_cells3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Baolin_T_NK_cells3.markers,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
Baolin_T_NK_cells3.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
pancancer_blueprint_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells_markers.rds"))
output_ove3 <- NULL
types_got <- as.character(unique(Baolin_T_NK_cells3.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Baolin_T_NK_cells3.markers[Baolin_T_NK_cells3.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_blueprint_T_NK_cells_markers[pancancer_blueprint_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove3 <- c(output_ove3,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove3) <- paste0("pan- blueprint: ",types_got)
rm("Baolin_T_NK_cells3")

# Baolin_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells4) <- Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells
# Baolin_T_NK_cells4.markers <- FindAllMarkers(Baolin_T_NK_cells4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Baolin_T_NK_cells4.markers,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
Baolin_T_NK_cells4.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
ABTC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells_markers.rds"))
output_ove4 <- NULL
types_got <- as.character(unique(Baolin_T_NK_cells4.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Baolin_T_NK_cells4.markers[Baolin_T_NK_cells4.markers$cluster==types_got[i],"gene"])
  Y <- unique(ABTC_T_NK_cells_markers[ABTC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove4 <- c(output_ove4,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove4) <- paste0("ABTC: ",types_got)
rm("Baolin_T_NK_cells4")

# Baolin_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells5) <- Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells
# Baolin_T_NK_cells5.markers <- FindAllMarkers(Baolin_T_NK_cells5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Baolin_T_NK_cells5.markers,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
Baolin_T_NK_cells5.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
HCC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells_markers.rds"))
output_ove5 <- NULL
types_got <- as.character(unique(Baolin_T_NK_cells5.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Baolin_T_NK_cells5.markers[Baolin_T_NK_cells5.markers$cluster==types_got[i],"gene"])
  Y <- unique(HCC_T_NK_cells_markers[HCC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove5 <- c(output_ove5,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove5) <- paste0("HCC: ",types_got)
rm("Baolin_T_NK_cells5")


# Baolin_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells6) <- Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells
# Baolin_T_NK_cells6.markers <- FindAllMarkers(Baolin_T_NK_cells6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Baolin_T_NK_cells6.markers,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
Baolin_T_NK_cells6.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
STAD_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells_markers.rds"))
output_ove6 <- NULL
types_got <- as.character(unique(Baolin_T_NK_cells6.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Baolin_T_NK_cells6.markers[Baolin_T_NK_cells6.markers$cluster==types_got[i],"gene"])
  Y <- unique(STAD_T_NK_cells_markers[STAD_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove6 <- c(output_ove6,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove6) <- paste0("STAD: ",types_got)
rm("Baolin_T_NK_cells6")
# c(mean(output_ove1),mean(output_ove2),mean(output_ove3),
#   mean(output_ove4),mean(output_ove5),mean(output_ove6))
Baolin_redetection <- c(output_ove1,output_ove2,output_ove3,output_ove4,output_ove5,output_ove6)

library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells
# Moshe_T_NK_cells1.markers <- FindAllMarkers(Moshe_T_NK_cells1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Moshe_T_NK_cells1.markers,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
Moshe_T_NK_cells1.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
TICA_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells_markers.rds"))
output_ove1 <- NULL
types_got <- as.character(unique(Moshe_T_NK_cells1.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Moshe_T_NK_cells1.markers[Moshe_T_NK_cells1.markers$cluster==types_got[i],"gene"])
  Y <- unique(TICA_T_NK_cells_markers[TICA_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove1 <- c(output_ove1,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove1) <- paste0("TICA: ",types_got)
rm("Moshe_T_NK_cells1")

# Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Moshe_T_NK_cells2) <- Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
# Moshe_T_NK_cells2.markers <- FindAllMarkers(Moshe_T_NK_cells2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Moshe_T_NK_cells2.markers,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
Moshe_T_NK_cells2.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
pancancer_T_cells_CD8_CD4_markers <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_markers_simplized.rds"))
output_ove2 <- NULL
types_got <- as.character(unique(Moshe_T_NK_cells2.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Moshe_T_NK_cells2.markers[Moshe_T_NK_cells2.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_T_cells_CD8_CD4_markers[pancancer_T_cells_CD8_CD4_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove2 <- c(output_ove2,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove2) <- types_got
rm("Moshe_T_NK_cells2")

# Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells3) <- Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
# Moshe_T_NK_cells3.markers <- FindAllMarkers(Moshe_T_NK_cells3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Moshe_T_NK_cells3.markers,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
Moshe_T_NK_cells3.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
pancancer_blueprint_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells_markers.rds"))
output_ove3 <- NULL
types_got <- as.character(unique(Moshe_T_NK_cells3.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Moshe_T_NK_cells3.markers[Moshe_T_NK_cells3.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_blueprint_T_NK_cells_markers[pancancer_blueprint_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove3 <- c(output_ove3,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove3) <- paste0("pan- blueprint: ",types_got)
rm("Moshe_T_NK_cells3")

# Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells4) <- Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells
# Moshe_T_NK_cells4.markers <- FindAllMarkers(Moshe_T_NK_cells4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Moshe_T_NK_cells4.markers,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
Moshe_T_NK_cells4.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
ABTC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells_markers.rds"))
output_ove4 <- NULL
types_got <- as.character(unique(Moshe_T_NK_cells4.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Moshe_T_NK_cells4.markers[Moshe_T_NK_cells4.markers$cluster==types_got[i],"gene"])
  Y <- unique(ABTC_T_NK_cells_markers[ABTC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove4 <- c(output_ove4,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove4) <- paste0("ABTC: ",types_got)
rm("Moshe_T_NK_cells4")

# Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells5) <- Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells
# Moshe_T_NK_cells5.markers <- FindAllMarkers(Moshe_T_NK_cells5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Moshe_T_NK_cells5.markers,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
Moshe_T_NK_cells5.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
HCC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells_markers.rds"))
output_ove5 <- NULL
types_got <- as.character(unique(Moshe_T_NK_cells5.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Moshe_T_NK_cells5.markers[Moshe_T_NK_cells5.markers$cluster==types_got[i],"gene"])
  Y <- unique(HCC_T_NK_cells_markers[HCC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove5 <- c(output_ove5,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove5) <- paste0("HCC: ",types_got)
rm("Moshe_T_NK_cells5")


# Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells6) <- Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells
# Moshe_T_NK_cells6.markers <- FindAllMarkers(Moshe_T_NK_cells6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Moshe_T_NK_cells6.markers,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
Moshe_T_NK_cells6.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
STAD_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells_markers.rds"))
output_ove6 <- NULL
types_got <- as.character(unique(Moshe_T_NK_cells6.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Moshe_T_NK_cells6.markers[Moshe_T_NK_cells6.markers$cluster==types_got[i],"gene"])
  Y <- unique(STAD_T_NK_cells_markers[STAD_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove6 <- c(output_ove6,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove6) <- paste0("STAD: ",types_got)
rm("Moshe_T_NK_cells6")
# c(mean(output_ove1),mean(output_ove2),mean(output_ove3),
#   mean(output_ove4),mean(output_ove5),mean(output_ove6))
Moshe_redetection <- c(output_ove1,output_ove2,output_ove3,output_ove4,output_ove5,output_ove6)


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells
# Caushi_T_NK_cells1.markers <- FindAllMarkers(Caushi_T_NK_cells1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Caushi_T_NK_cells1.markers,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
Caushi_T_NK_cells1.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
TICA_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells_markers.rds"))
output_ove1 <- NULL
types_got <- as.character(unique(Caushi_T_NK_cells1.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Caushi_T_NK_cells1.markers[Caushi_T_NK_cells1.markers$cluster==types_got[i],"gene"])
  Y <- unique(TICA_T_NK_cells_markers[TICA_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove1 <- c(output_ove1,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove1) <- paste0("TICA: ",types_got)
rm("Caushi_T_NK_cells1")

# Caushi_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Caushi_T_NK_cells2) <- Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
# Caushi_T_NK_cells2.markers <- FindAllMarkers(Caushi_T_NK_cells2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Caushi_T_NK_cells2.markers,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
Caushi_T_NK_cells2.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
pancancer_T_cells_CD8_CD4_markers <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_markers_simplized.rds"))
output_ove2 <- NULL
types_got <- as.character(unique(Caushi_T_NK_cells2.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Caushi_T_NK_cells2.markers[Caushi_T_NK_cells2.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_T_cells_CD8_CD4_markers[pancancer_T_cells_CD8_CD4_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove2 <- c(output_ove2,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove2) <- types_got
rm("Caushi_T_NK_cells2")

# Caushi_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells3) <- Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
# Caushi_T_NK_cells3.markers <- FindAllMarkers(Caushi_T_NK_cells3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Caushi_T_NK_cells3.markers,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
Caushi_T_NK_cells3.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
pancancer_blueprint_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells_markers.rds"))
output_ove3 <- NULL
types_got <- as.character(unique(Caushi_T_NK_cells3.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Caushi_T_NK_cells3.markers[Caushi_T_NK_cells3.markers$cluster==types_got[i],"gene"])
  Y <- unique(pancancer_blueprint_T_NK_cells_markers[pancancer_blueprint_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove3 <- c(output_ove3,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove3) <- paste0("pan- blueprint: ",types_got)
rm("Caushi_T_NK_cells3")

# Caushi_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells4) <- Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells
# Caushi_T_NK_cells4.markers <- FindAllMarkers(Caushi_T_NK_cells4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Caushi_T_NK_cells4.markers,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
Caushi_T_NK_cells4.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
ABTC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells_markers.rds"))
output_ove4 <- NULL
types_got <- as.character(unique(Caushi_T_NK_cells4.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Caushi_T_NK_cells4.markers[Caushi_T_NK_cells4.markers$cluster==types_got[i],"gene"])
  Y <- unique(ABTC_T_NK_cells_markers[ABTC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove4 <- c(output_ove4,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove4) <- paste0("ABTC: ",types_got)
rm("Caushi_T_NK_cells4")

# Caushi_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells5) <- Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells
# Caushi_T_NK_cells5.markers <- FindAllMarkers(Caushi_T_NK_cells5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Caushi_T_NK_cells5.markers,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
Caushi_T_NK_cells5.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
HCC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells_markers.rds"))
output_ove5 <- NULL
types_got <- as.character(unique(Caushi_T_NK_cells5.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Caushi_T_NK_cells5.markers[Caushi_T_NK_cells5.markers$cluster==types_got[i],"gene"])
  Y <- unique(HCC_T_NK_cells_markers[HCC_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove5 <- c(output_ove5,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove5) <- paste0("HCC: ",types_got)
rm("Caushi_T_NK_cells5")


# Caushi_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells6) <- Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells
# Caushi_T_NK_cells6.markers <- FindAllMarkers(Caushi_T_NK_cells6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Caushi_T_NK_cells6.markers,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
Caushi_T_NK_cells6.markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
STAD_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells_markers.rds"))
output_ove6 <- NULL
types_got <- as.character(unique(Caushi_T_NK_cells6.markers$cluster))
for (i in 1:length(types_got)) {
  X <- unique(Caushi_T_NK_cells6.markers[Caushi_T_NK_cells6.markers$cluster==types_got[i],"gene"])
  Y <- unique(STAD_T_NK_cells_markers[STAD_T_NK_cells_markers$cluster==types_got[i],"gene"])
  X_Y <- intersect(X,Y)
  output_ove6 <- c(output_ove6,2*length(X_Y)/(length(X)+length(Y)))
}
names(output_ove6) <- paste0("STAD: ",types_got)
rm("Caushi_T_NK_cells6")
# c(mean(output_ove1),mean(output_ove2),mean(output_ove3),
#   mean(output_ove4),mean(output_ove5),mean(output_ove6))
Caushi_redetection <- c(output_ove1,output_ove2,output_ove3,output_ove4,output_ove5,output_ove6)

# output <-rbind.fill(output_redetection,as.data.frame(t(as.matrix(output_correlation_within_population))))
all_types <- c(paste0("TICA: ",unique(TICA_T_NK_cells_markers$cluster)),#paste0("pancancer_T_CD8_CD4: ",unique(pancancer_T_cells_CD8_CD4_markers$cluster)),
               paste0("pan- blueprint: ",unique(pancancer_blueprint_T_NK_cells_markers$cluster)),
               paste0("ABTC: ",unique(ABTC_T_NK_cells_markers$cluster)),
               paste0("HCC: ",unique(HCC_T_NK_cells_markers$cluster)),
               paste0("STAD: ",unique(STAD_T_NK_cells_markers$cluster)))
names(all_types) <- all_types
output_redetection <- rbind.fill(as.data.frame(t(as.matrix(all_types))),as.data.frame(t(as.matrix(Kathryn_redetection))),as.data.frame(t(as.matrix(Baolin_redetection))),
           as.data.frame(t(as.matrix(Moshe_redetection))),as.data.frame(t(as.matrix(Caushi_redetection))))


#####################################################
#### 2023.01.12 correlations within mapped populations
#### change "TICA_T_NK_cells" to others
#####################################################

# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#     Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$predicted_to_TICA_T_NK_cells@data
#     Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#     Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#     rm(Kathryn_sub_exprs)
#     Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$predicted_to_TICA_T_NK_cells@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$predicted_to_TICA_T_NK_cells@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$predicted_to_TICA_T_NK_cells@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$predicted_to_pancancer_blueprint_T_NK_cells@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$predicted_to_pancancer_blueprint_T_NK_cells@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$predicted_to_pancancer_blueprint_T_NK_cells@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$predicted_to_pancancer_blueprint_T_NK_cells@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$predicted_to_ABTC_T_NK_cells@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$predicted_to_ABTC_T_NK_cells@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$predicted_to_ABTC_T_NK_cells@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$predicted_to_ABTC_T_NK_cells@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# ##
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# Baolin_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# Moshe_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# Caushi_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_cluster_mean.rds"))
# # types <- unique(c(names(Kathryn_TICA),names(Baolin_TICA),names(Moshe_TICA),names(Caushi_TICA)))
# types_TICA <- unique(intersect(intersect(names(Kathryn_TICA),names(Baolin_TICA)),intersect(names(Moshe_TICA),names(Caushi_TICA))))
# output_TICA <- NULL
# for (i in 1:length(types_TICA)) {
#   tmp1 <- merge(as.data.frame(Kathryn_TICA[[types_TICA[i]]]),as.data.frame(Baolin_TICA[[types_TICA[i]]]),by="row.names")
#   tmp2 <- merge(tmp1,as.data.frame(Moshe_TICA[[types_TICA[i]]]),by.x="Row.names",by.y="row.names")
#   tmp3 <- merge(tmp2,as.data.frame(Caushi_TICA[[types_TICA[i]]]),by.x="Row.names",by.y="row.names")
#   tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
#   tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
#   output_TICA <- c(output_TICA,tmp_cor_mean)
# }
# mean(output_TICA)
#   
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# Baolin_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# Moshe_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# Caushi_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_mean.rds"))
# # types <- unique(c(names(Kathryn_pancancer_blueprint),names(Baolin_pancancer_blueprint),names(Moshe_pancancer_blueprint),names(Caushi_pancancer_blueprint)))
# types_pancancer_blueprint <- unique(intersect(intersect(names(Kathryn_pancancer_blueprint),names(Baolin_pancancer_blueprint)),intersect(names(Moshe_pancancer_blueprint),names(Caushi_pancancer_blueprint))))
# output_pancancer_blueprint <- NULL
# for (i in 1:length(types_pancancer_blueprint)) {
#   tmp1 <- merge(as.data.frame(Kathryn_pancancer_blueprint[[types_pancancer_blueprint[i]]]),as.data.frame(Baolin_pancancer_blueprint[[types_pancancer_blueprint[i]]]),by="row.names")
#   tmp2 <- merge(tmp1,as.data.frame(Moshe_pancancer_blueprint[[types_pancancer_blueprint[i]]]),by.x="Row.names",by.y="row.names")
#   tmp3 <- merge(tmp2,as.data.frame(Caushi_pancancer_blueprint[[types_pancancer_blueprint[i]]]),by.x="Row.names",by.y="row.names")
#   tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
#   tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
#   output_pancancer_blueprint <- c(output_pancancer_blueprint,tmp_cor_mean)
# }
# mean(output_pancancer_blueprint)
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# Baolin_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# Moshe_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# Caushi_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_cluster_mean.rds"))
# # types <- unique(c(names(Kathryn_ABTC),names(Baolin_ABTC),names(Moshe_ABTC),names(Caushi_ABTC)))
# types_ABTC <- unique(intersect(intersect(names(Kathryn_ABTC),names(Baolin_ABTC)),intersect(names(Moshe_ABTC),names(Caushi_ABTC))))
# output_ABTC <- NULL
# for (i in 1:length(types_ABTC)) {
#   tmp1 <- merge(as.data.frame(Kathryn_ABTC[[types_ABTC[i]]]),as.data.frame(Baolin_ABTC[[types_ABTC[i]]]),by="row.names")
#   tmp2 <- merge(tmp1,as.data.frame(Moshe_ABTC[[types_ABTC[i]]]),by.x="Row.names",by.y="row.names")
#   tmp3 <- merge(tmp2,as.data.frame(Caushi_ABTC[[types_ABTC[i]]]),by.x="Row.names",by.y="row.names")
#   tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
#   tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
#   output_ABTC <- c(output_ABTC,tmp_cor_mean)
# }
# mean(output_ABTC)
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells_cluster_mean.rds"))
# Baolin_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells_cluster_mean.rds"))
# Moshe_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells_cluster_mean.rds"))
# Caushi_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells_cluster_mean.rds"))
# # types <- unique(c(names(Kathryn_HCC),names(Baolin_HCC),names(Moshe_HCC),names(Caushi_HCC)))
# types_HCC <- unique(intersect(intersect(names(Kathryn_HCC),names(Baolin_HCC)),intersect(names(Moshe_HCC),names(Caushi_HCC))))
# output_HCC <- NULL
# for (i in 1:length(types_HCC)) {
#   tmp1 <- merge(as.data.frame(Kathryn_HCC[[types_HCC[i]]]),as.data.frame(Baolin_HCC[[types_HCC[i]]]),by="row.names")
#   tmp2 <- merge(tmp1,as.data.frame(Moshe_HCC[[types_HCC[i]]]),by.x="Row.names",by.y="row.names")
#   tmp3 <- merge(tmp2,as.data.frame(Caushi_HCC[[types_HCC[i]]]),by.x="Row.names",by.y="row.names")
#   tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
#   tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
#   output_HCC <- c(output_HCC,tmp_cor_mean)
# }
# mean(output_HCC)
# 
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells_cluster_mean.rds"))
# Baolin_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells_cluster_mean.rds"))
# Moshe_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells_cluster_mean.rds"))
# Caushi_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells_cluster_mean.rds"))
# # types <- unique(c(names(Kathryn_STAD),names(Baolin_STAD),names(Moshe_STAD),names(Caushi_STAD)))
# types_STAD <- unique(intersect(intersect(names(Kathryn_STAD),names(Baolin_STAD)),intersect(names(Moshe_STAD),names(Caushi_STAD))))
# output_STAD <- NULL
# for (i in 1:length(types_STAD)) {
#   tmp1 <- merge(as.data.frame(Kathryn_STAD[[types_STAD[i]]]),as.data.frame(Baolin_STAD[[types_STAD[i]]]),by="row.names")
#   tmp2 <- merge(tmp1,as.data.frame(Moshe_STAD[[types_STAD[i]]]),by.x="Row.names",by.y="row.names")
#   tmp3 <- merge(tmp2,as.data.frame(Caushi_STAD[[types_STAD[i]]]),by.x="Row.names",by.y="row.names")
#   tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
#   tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
#   output_STAD <- c(output_STAD,tmp_cor_mean)
# }
# mean(output_STAD)



#####################################################
#### 2023.01.19 correlations within mapped populations, use original exprs, not predicted reference exprs
#### change "TICA_T_NK_cells" to others
#####################################################

# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_pancancer_T_cells_CD8_CD4
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_pancancer_T_cells_CD8_CD4
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_pancancer_T_cells_CD8_CD4
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got[1:7]
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_pancancer_T_cells_CD8_CD4
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_pancancer_blueprint_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_ABTC_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_HCC_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_HCC_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_HCC_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_HCC_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
# rm(Caushi_T_NK_cells1)
# 
# library(Seurat)
# library(SeuratObject)
# library(lisi)
# library(mclust)
# library(cluster)
# workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Kathryn_T_NK_cells1) <- Kathryn_T_NK_cells1$predicted.to_STAD_T_NK_cells
# types_got <- as.character(unique(Idents(Kathryn_T_NK_cells1)))
# Kathryn_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Kathryn_sub_exprs <- subset(Kathryn_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Kathryn_sub_exprs <- Kathryn_sub_exprs[apply(Kathryn_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Kathryn_sub_sum <- apply(Kathryn_sub_exprs,1,mean)
#   rm(Kathryn_sub_exprs)
#   Kathryn_sub_exprs_mean[[i]] <- Kathryn_sub_sum
# }
# names(Kathryn_sub_exprs_mean) <- types_got
# saveRDS(Kathryn_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
# rm(Kathryn_T_NK_cells1)
# 
# 
# Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Baolin_T_NK_cells1) <- Baolin_T_NK_cells1$predicted.to_STAD_T_NK_cells
# # Baolin_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Baolin_T_NK_cells1)))
# Baolin_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Baolin_sub_exprs <- subset(Baolin_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Baolin_sub_exprs <- Baolin_sub_exprs[apply(Baolin_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Baolin_sub_sum <- apply(Baolin_sub_exprs,1,mean)
#   rm(Baolin_sub_exprs)
#   Baolin_sub_exprs_mean[[i]] <- Baolin_sub_sum
# }
# names(Baolin_sub_exprs_mean) <- types_got
# saveRDS(Baolin_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
# rm(Baolin_T_NK_cells1)
# 
# Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Moshe_T_NK_cells1) <- Moshe_T_NK_cells1$predicted.to_STAD_T_NK_cells
# # Moshe_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Moshe_T_NK_cells1)))
# Moshe_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Moshe_sub_exprs <- subset(Moshe_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Moshe_sub_exprs <- Moshe_sub_exprs[apply(Moshe_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Moshe_sub_sum <- apply(Moshe_sub_exprs,1,mean)
#   rm(Moshe_sub_exprs)
#   Moshe_sub_exprs_mean[[i]] <- Moshe_sub_sum
# }
# names(Moshe_sub_exprs_mean) <- types_got
# saveRDS(Moshe_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
# rm(Moshe_T_NK_cells1)
# 
# 
# Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells.rds"))
# Idents(Caushi_T_NK_cells1) <- Caushi_T_NK_cells1$predicted.to_STAD_T_NK_cells
# # Caushi_T_NK_cells1_markers <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells_markers.rds"))
# types_got <- as.character(unique(Idents(Caushi_T_NK_cells1)))
# Caushi_sub_exprs_mean <- list()
# for (i in 1:length(types_got)) {
#   Caushi_sub_exprs <- subset(Caushi_T_NK_cells1,idents=types_got[i])@assays$RNA@data
#   Caushi_sub_exprs <- Caushi_sub_exprs[apply(Caushi_sub_exprs,1,function(x) sum(x==0)<length(x)/2),]
#   Caushi_sub_sum <- apply(Caushi_sub_exprs,1,mean)
#   rm(Caushi_sub_exprs)
#   Caushi_sub_exprs_mean[[i]] <- Caushi_sub_sum
# }
# names(Caushi_sub_exprs_mean) <- types_got
# saveRDS(Caushi_sub_exprs_mean,file=paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
# rm(Caushi_T_NK_cells1)


##
library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
Kathryn_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
Baolin_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
Moshe_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
Caushi_TICA <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells_cluster_original_mean.rds"))
types_TICA <- unique(c(names(Kathryn_TICA),names(Baolin_TICA),names(Moshe_TICA),names(Caushi_TICA)))
types_TICAbubian <- types_TICA
# types_TICA <- unique(intersect(intersect(names(Kathryn_TICA),names(Baolin_TICA)),intersect(names(Moshe_TICA),names(Caushi_TICA))))
output_TICA <- NULL
for (i in 1:length(types_TICA)) {
  if (types_TICA[i] %in% names(Kathryn_TICA)) {
    tmp_Kathryn <- Kathryn_TICA[[types_TICA[i]]]
  } else {
    tmp_Kathryn <- NA
  }
  if (types_TICA[i] %in% names(Baolin_TICA)) {
    tmp_Baolin <- Baolin_TICA[[types_TICA[i]]]
  } else {
    tmp_Baolin <- NA
  }
  if (types_TICA[i] %in% names(Moshe_TICA)) {
    tmp_Moshe <- Moshe_TICA[[types_TICA[i]]]
  } else {
    tmp_Moshe <- NA
  }
  if (types_TICA[i] %in% names(Caushi_TICA)) {
    tmp_Caushi <- Caushi_TICA[[types_TICA[i]]]
  } else {
    tmp_Caushi <- NA
  }
  tmp1 <-  merge(as.data.frame(tmp_Kathryn),as.data.frame(tmp_Baolin),by="row.names",all=T)
  tmp2 <- merge(tmp1,as.data.frame(tmp_Moshe),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- merge(tmp2,as.data.frame(tmp_Caushi),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- tmp3[,2:ncol(tmp3)][,apply(tmp3[,2:ncol(tmp3)],2,function(x) sum(is.na(x)) < length(x))]
  if (!is.null(ncol(tmp3))) {
    if (ncol(tmp3) >1) {
      tmp3 <- tmp3[apply(tmp3,1,function(x) sum(!is.na(x)) == length(x)),]
      tmp_cor <- cor(tmp3)
      tmp_cor_mean <- mean(tmp_cor)#((sum(tmp_cor)-ncol(tmp_cor))/2)/ncol(tmp_cor)
      output_TICA <- c(output_TICA,tmp_cor_mean)
    }
  } else {
    types_TICAbubian <- types_TICAbubian[-i]
  }
  # tmp1 <- merge(as.data.frame(Kathryn_TICA[[types_TICA[i]]]),as.data.frame(Baolin_TICA[[types_TICA[i]]]),by="row.names")
  # tmp2 <- merge(tmp1,as.data.frame(Moshe_TICA[[types_TICA[i]]]),by.x="Row.names",by.y="row.names")
  # tmp3 <- merge(tmp2,as.data.frame(Caushi_TICA[[types_TICA[i]]]),by.x="Row.names",by.y="row.names")
  # tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
  # tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
  # output_TICA <- c(output_TICA,tmp_cor_mean)
}
names(output_TICA) <- paste0("TICA: ",types_TICAbubian)


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
Kathryn_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
Baolin_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
Moshe_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
Caushi_pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4_cluster_original_mean.rds"))
types_pancancer_T_cells_CD8_CD4 <- unique(c(names(Kathryn_pancancer_T_cells_CD8_CD4),names(Baolin_pancancer_T_cells_CD8_CD4),names(Moshe_pancancer_T_cells_CD8_CD4),names(Caushi_pancancer_T_cells_CD8_CD4)))
types_pancancer_T_cells_CD8_CD4bubian <- types_pancancer_T_cells_CD8_CD4
# types_pancancer_T_cells_CD8_CD4 <- unique(intersect(intersect(names(Kathryn_pancancer_T_cells_CD8_CD4),names(Baolin_pancancer_T_cells_CD8_CD4)),intersect(names(Moshe_pancancer_T_cells_CD8_CD4),names(Caushi_pancancer_T_cells_CD8_CD4))))
output_pancancer_T_cells_CD8_CD4 <- NULL
for (i in 1:length(types_pancancer_T_cells_CD8_CD4)) {
  if (types_pancancer_T_cells_CD8_CD4[i] %in% names(Kathryn_pancancer_T_cells_CD8_CD4)) {
    tmp_Kathryn <- Kathryn_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]
  } else {
    tmp_Kathryn <- NA
  }
  if (types_pancancer_T_cells_CD8_CD4[i] %in% names(Baolin_pancancer_T_cells_CD8_CD4)) {
    tmp_Baolin <- Baolin_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]
  } else {
    tmp_Baolin <- NA
  }
  if (types_pancancer_T_cells_CD8_CD4[i] %in% names(Moshe_pancancer_T_cells_CD8_CD4)) {
    tmp_Moshe <- Moshe_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]
  } else {
    tmp_Moshe <- NA
  }
  if (types_pancancer_T_cells_CD8_CD4[i] %in% names(Caushi_pancancer_T_cells_CD8_CD4)) {
    tmp_Caushi <- Caushi_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]
  } else {
    tmp_Caushi <- NA
  }
  tmp1 <-  merge(as.data.frame(tmp_Kathryn),as.data.frame(tmp_Baolin),by="row.names",all=T)
  tmp2 <- merge(tmp1,as.data.frame(tmp_Moshe),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- merge(tmp2,as.data.frame(tmp_Caushi),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- tmp3[,2:ncol(tmp3)][,apply(tmp3[,2:ncol(tmp3)],2,function(x) sum(is.na(x)) < length(x))]
  if (!is.null(ncol(tmp3))) {
    if (ncol(tmp3) >1) {
      tmp3 <- tmp3[apply(tmp3,1,function(x) sum(!is.na(x)) == length(x)),]
      tmp_cor <- cor(tmp3)
      tmp_cor_mean <- mean(tmp_cor)#((sum(tmp_cor)-ncol(tmp_cor))/2)/ncol(tmp_cor)
      output_pancancer_T_cells_CD8_CD4 <- c(output_pancancer_T_cells_CD8_CD4,tmp_cor_mean)
    }
  } else {
    types_pancancer_T_cells_CD8_CD4bubian <- types_pancancer_T_cells_CD8_CD4bubian[-i]
  }
  # tmp1 <- merge(as.data.frame(Kathryn_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]),as.data.frame(Baolin_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]),by="row.names")
  # tmp2 <- merge(tmp1,as.data.frame(Moshe_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]),by.x="Row.names",by.y="row.names")
  # tmp3 <- merge(tmp2,as.data.frame(Caushi_pancancer_T_cells_CD8_CD4[[types_pancancer_T_cells_CD8_CD4[i]]]),by.x="Row.names",by.y="row.names")
  # tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
  # tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
  # output_pancancer_T_cells_CD8_CD4 <- c(output_pancancer_T_cells_CD8_CD4,tmp_cor_mean)
}
names(output_pancancer_T_cells_CD8_CD4) <- types_pancancer_T_cells_CD8_CD4bubian


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
Kathryn_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
Baolin_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
Moshe_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
Caushi_pancancer_blueprint <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells_cluster_original_mean.rds"))
types_pancancer_blueprint <- unique(c(names(Kathryn_pancancer_blueprint),names(Baolin_pancancer_blueprint),names(Moshe_pancancer_blueprint),names(Caushi_pancancer_blueprint)))
types_pancancer_blueprintbubian <- types_pancancer_blueprint
# types_pancancer_blueprint <- unique(intersect(intersect(names(Kathryn_pancancer_blueprint),names(Baolin_pancancer_blueprint)),intersect(names(Moshe_pancancer_blueprint),names(Caushi_pancancer_blueprint))))
output_pancancer_blueprint <- NULL
for (i in 1:length(types_pancancer_blueprint)) {
  if (types_pancancer_blueprint[i] %in% names(Kathryn_pancancer_blueprint)) {
    tmp_Kathryn <- Kathryn_pancancer_blueprint[[types_pancancer_blueprint[i]]]
  } else {
    tmp_Kathryn <- NA
  }
  if (types_pancancer_blueprint[i] %in% names(Baolin_pancancer_blueprint)) {
    tmp_Baolin <- Baolin_pancancer_blueprint[[types_pancancer_blueprint[i]]]
  } else {
    tmp_Baolin <- NA
  }
  if (types_pancancer_blueprint[i] %in% names(Moshe_pancancer_blueprint)) {
    tmp_Moshe <- Moshe_pancancer_blueprint[[types_pancancer_blueprint[i]]]
  } else {
    tmp_Moshe <- NA
  }
  if (types_pancancer_blueprint[i] %in% names(Caushi_pancancer_blueprint)) {
    tmp_Caushi <- Caushi_pancancer_blueprint[[types_pancancer_blueprint[i]]]
  } else {
    tmp_Caushi <- NA
  }
  tmp1 <-  merge(as.data.frame(tmp_Kathryn),as.data.frame(tmp_Baolin),by="row.names",all=T)
  tmp2 <- merge(tmp1,as.data.frame(tmp_Moshe),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- merge(tmp2,as.data.frame(tmp_Caushi),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- tmp3[,2:ncol(tmp3)][,apply(tmp3[,2:ncol(tmp3)],2,function(x) sum(is.na(x)) < length(x))]
  if (!is.null(ncol(tmp3))) {
    if (ncol(tmp3) >1) {
      tmp3 <- tmp3[apply(tmp3,1,function(x) sum(!is.na(x)) == length(x)),]
      tmp_cor <- cor(tmp3)
      tmp_cor_mean <- mean(tmp_cor)#((sum(tmp_cor)-ncol(tmp_cor))/2)/ncol(tmp_cor)
      output_pancancer_blueprint <- c(output_pancancer_blueprint,tmp_cor_mean)
    }
  } else {
    types_pancancer_blueprintbubian <- types_pancancer_blueprintbubian[-i]
  }
  # tmp1 <- merge(as.data.frame(Kathryn_pancancer_blueprint[[types_pancancer_blueprint[i]]]),as.data.frame(Baolin_pancancer_blueprint[[types_pancancer_blueprint[i]]]),by="row.names")
  # tmp2 <- merge(tmp1,as.data.frame(Moshe_pancancer_blueprint[[types_pancancer_blueprint[i]]]),by.x="Row.names",by.y="row.names")
  # tmp3 <- merge(tmp2,as.data.frame(Caushi_pancancer_blueprint[[types_pancancer_blueprint[i]]]),by.x="Row.names",by.y="row.names")
  # tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
  # tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
  # output_pancancer_blueprint <- c(output_pancancer_blueprint,tmp_cor_mean)
}
names(output_pancancer_blueprint) <- paste0("pan- blueprint: ",types_pancancer_blueprintbubian)

library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
Kathryn_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
Baolin_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
Moshe_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
Caushi_ABTC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells_cluster_original_mean.rds"))
types_ABTC <- unique(c(names(Kathryn_ABTC),names(Baolin_ABTC),names(Moshe_ABTC),names(Caushi_ABTC)))
types_ABTCbubian <- types_ABTC
# types_ABTC <- unique(intersect(intersect(names(Kathryn_ABTC),names(Baolin_ABTC)),intersect(names(Moshe_ABTC),names(Caushi_ABTC))))
output_ABTC <- NULL
for (i in 1:length(types_ABTC)) {
  if (types_ABTC[i] %in% names(Kathryn_ABTC)) {
    tmp_Kathryn <- Kathryn_ABTC[[types_ABTC[i]]]
  } else {
    tmp_Kathryn <- NA
  }
  if (types_ABTC[i] %in% names(Baolin_ABTC)) {
    tmp_Baolin <- Baolin_ABTC[[types_ABTC[i]]]
  } else {
    tmp_Baolin <- NA
  }
  if (types_ABTC[i] %in% names(Moshe_ABTC)) {
    tmp_Moshe <- Moshe_ABTC[[types_ABTC[i]]]
  } else {
    tmp_Moshe <- NA
  }
  if (types_ABTC[i] %in% names(Caushi_ABTC)) {
    tmp_Caushi <- Caushi_ABTC[[types_ABTC[i]]]
  } else {
    tmp_Caushi <- NA
  }
  tmp1 <-  merge(as.data.frame(tmp_Kathryn),as.data.frame(tmp_Baolin),by="row.names",all=T)
  tmp2 <- merge(tmp1,as.data.frame(tmp_Moshe),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- merge(tmp2,as.data.frame(tmp_Caushi),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- tmp3[,2:ncol(tmp3)][,apply(tmp3[,2:ncol(tmp3)],2,function(x) sum(is.na(x)) < length(x))]
  if (!is.null(ncol(tmp3))) {
    if (ncol(tmp3) >1) {
      tmp3 <- tmp3[apply(tmp3,1,function(x) sum(!is.na(x)) == length(x)),]
      tmp_cor <- cor(tmp3)
      tmp_cor_mean <- mean(tmp_cor)#((sum(tmp_cor)-ncol(tmp_cor))/2)/ncol(tmp_cor)
      output_ABTC <- c(output_ABTC,tmp_cor_mean)
    }
  } else {
    types_ABTCbubian <- types_ABTCbubian[-i]
  }
  # tmp1 <- merge(as.data.frame(Kathryn_ABTC[[types_ABTC[i]]]),as.data.frame(Baolin_ABTC[[types_ABTC[i]]]),by="row.names")
  # tmp2 <- merge(tmp1,as.data.frame(Moshe_ABTC[[types_ABTC[i]]]),by.x="Row.names",by.y="row.names")
  # tmp3 <- merge(tmp2,as.data.frame(Caushi_ABTC[[types_ABTC[i]]]),by.x="Row.names",by.y="row.names")
  # tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
  # tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
  # output_ABTC <- c(output_ABTC,tmp_cor_mean)
}
names(output_ABTC) <- paste0("ABTC: ",types_ABTCbubian)

library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
Kathryn_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
Baolin_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
Moshe_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
Caushi_HCC <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells_cluster_original_mean.rds"))
types_HCC <- unique(c(names(Kathryn_HCC),names(Baolin_HCC),names(Moshe_HCC),names(Caushi_HCC)))
types_HCCbubian <- types_HCC
# types_HCC <- unique(intersect(intersect(names(Kathryn_HCC),names(Baolin_HCC)),intersect(names(Moshe_HCC),names(Caushi_HCC))))
output_HCC <- NULL
for (i in 1:length(types_HCC)) {
  if (types_HCC[i] %in% names(Kathryn_HCC)) {
    tmp_Kathryn <- Kathryn_HCC[[types_HCC[i]]]
  } else {
    tmp_Kathryn <- NA
  }
  if (types_HCC[i] %in% names(Baolin_HCC)) {
    tmp_Baolin <- Baolin_HCC[[types_HCC[i]]]
  } else {
    tmp_Baolin <- NA
  }
  if (types_HCC[i] %in% names(Moshe_HCC)) {
    tmp_Moshe <- Moshe_HCC[[types_HCC[i]]]
  } else {
    tmp_Moshe <- NA
  }
  if (types_HCC[i] %in% names(Caushi_HCC)) {
    tmp_Caushi <- Caushi_HCC[[types_HCC[i]]]
  } else {
    tmp_Caushi <- NA
  }
  tmp1 <-  merge(as.data.frame(tmp_Kathryn),as.data.frame(tmp_Baolin),by="row.names",all=T)
  tmp2 <- merge(tmp1,as.data.frame(tmp_Moshe),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- merge(tmp2,as.data.frame(tmp_Caushi),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- tmp3[,2:ncol(tmp3)][,apply(tmp3[,2:ncol(tmp3)],2,function(x) sum(is.na(x)) < length(x))]
  if (!is.null(ncol(tmp3))) {
    if (ncol(tmp3) >1) {
      tmp3 <- tmp3[apply(tmp3,1,function(x) sum(!is.na(x)) == length(x)),]
      tmp_cor <- cor(tmp3)
      tmp_cor_mean <- mean(tmp_cor)#((sum(tmp_cor)-ncol(tmp_cor))/2)/ncol(tmp_cor)
      output_HCC <- c(output_HCC,tmp_cor_mean)
    }
  } else {
    types_HCCbubian <- types_HCCbubian[-i]
  }
  # tmp1 <- merge(as.data.frame(Kathryn_HCC[[types_HCC[i]]]),as.data.frame(Baolin_HCC[[types_HCC[i]]]),by="row.names")
  # tmp2 <- merge(tmp1,as.data.frame(Moshe_HCC[[types_HCC[i]]]),by.x="Row.names",by.y="row.names")
  # tmp3 <- merge(tmp2,as.data.frame(Caushi_HCC[[types_HCC[i]]]),by.x="Row.names",by.y="row.names")
  # tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
  # tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
  # output_HCC <- c(output_HCC,tmp_cor_mean)
}
names(output_HCC) <- paste0("HCC: ",types_HCCbubian)


library(Seurat)
library(SeuratObject)
library(lisi)
library(mclust)
library(cluster)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
Kathryn_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
Baolin_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
Moshe_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
Caushi_STAD <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells_cluster_original_mean.rds"))
types_STAD <- unique(c(names(Kathryn_STAD),names(Baolin_STAD),names(Moshe_STAD),names(Caushi_STAD)))
types_STADbubian <- types_STAD
# types_STAD <- unique(intersect(intersect(names(Kathryn_STAD),names(Baolin_STAD)),intersect(names(Moshe_STAD),names(Caushi_STAD))))
output_STAD <- NULL
for (i in 1:length(types_STAD)) {
  if (types_STAD[i] %in% names(Kathryn_STAD)) {
    tmp_Kathryn <- Kathryn_STAD[[types_STAD[i]]]
  } else {
    tmp_Kathryn <- NA
  }
  if (types_STAD[i] %in% names(Baolin_STAD)) {
    tmp_Baolin <- Baolin_STAD[[types_STAD[i]]]
  } else {
    tmp_Baolin <- NA
  }
  if (types_STAD[i] %in% names(Moshe_STAD)) {
    tmp_Moshe <- Moshe_STAD[[types_STAD[i]]]
  } else {
    tmp_Moshe <- NA
  }
  if (types_STAD[i] %in% names(Caushi_STAD)) {
    tmp_Caushi <- Caushi_STAD[[types_STAD[i]]]
  } else {
    tmp_Caushi <- NA
  }
  tmp1 <-  merge(as.data.frame(tmp_Kathryn),as.data.frame(tmp_Baolin),by="row.names",all=T)
  tmp2 <- merge(tmp1,as.data.frame(tmp_Moshe),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- merge(tmp2,as.data.frame(tmp_Caushi),by.x="Row.names",by.y="row.names",all=T)
  tmp3 <- tmp3[,2:ncol(tmp3)][,apply(tmp3[,2:ncol(tmp3)],2,function(x) sum(is.na(x)) < length(x))]
  if (!is.null(ncol(tmp3))) {
    if (ncol(tmp3) >1) {
      tmp3 <- tmp3[apply(tmp3,1,function(x) sum(!is.na(x)) == length(x)),]
      tmp_cor <- cor(tmp3)
      tmp_cor_mean <- mean(tmp_cor)#((sum(tmp_cor)-ncol(tmp_cor))/2)/ncol(tmp_cor)
      output_STAD <- c(output_STAD,tmp_cor_mean)
    }
  } else {
    types_STADbubian <- types_STADbubian[-i]
  }
  # tmp1 <- merge(as.data.frame(Kathryn_STAD[[types_STAD[i]]]),as.data.frame(Baolin_STAD[[types_STAD[i]]]),by="row.names")
  # tmp2 <- merge(tmp1,as.data.frame(Moshe_STAD[[types_STAD[i]]]),by.x="Row.names",by.y="row.names")
  # tmp3 <- merge(tmp2,as.data.frame(Caushi_STAD[[types_STAD[i]]]),by.x="Row.names",by.y="row.names")
  # tmp_cor <- cor(tmp3[,2:ncol(tmp3)])
  # tmp_cor_mean <- ((sum(tmp_cor)-ncol(tmp_cor))/2)/6
  # output_STAD <- c(output_STAD,tmp_cor_mean)
}
names(output_STAD) <- paste0("STAD: ",types_STADbubian)

output_correlation_within_population <- c(output_TICA,
                                          # output_pancancer_T_cells_CD8_CD4,
                                          output_pancancer_blueprint,output_ABTC,output_HCC,output_STAD)
output <-rbind.fill(output_redetection,as.data.frame(t(as.matrix(output_correlation_within_population))))
output <- output[-1,]
colnames(output) <- as.character(new_cluster_names[colnames(output)])
write.table(output,file=paste0(workdir,"3_mapping_query_datasets_criteria_redetections_and_correlations_within_populations_20240319.txt"),quote = F,sep="\t",row.names = F)
