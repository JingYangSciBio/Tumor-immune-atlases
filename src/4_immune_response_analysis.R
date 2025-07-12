############################################
#### 2025.07.05 Jing Yang 
####  Caushi, Moshepre and Moshepost data, adjust covairant by age and gender
############################################
### for p values

Caushi<- read.table(paste0(workdir,"4_immune_response_Caushi_ageAndgender_to_reference.txt"),header=T,sep="\t",row.names = 1);dim(Moshepre)
colnames(Caushi) <- paste0("Caushi_",colnames(Caushi))
rownames(Caushi) <- as.character(as.data.frame(strsplit(rownames(Caushi),".predicted.to_"))[2,])

Moshepre <- read.table(paste0(workdir,"4_immune_response_Moshepre_ageAndgender_to_reference.txt"),header=T,sep="\t",row.names = 1);dim(Moshepre)
colnames(Moshepre) <- paste0("Moshepre_",colnames(Moshepre))
rownames(Moshepre) <- as.character(as.data.frame(strsplit(rownames(Moshepre),".predicted.to_"))[2,])
Moshepost <- read.table(paste0(workdir,"4_immune_response_Moshepost_ageAndgender_to_reference.txt"),header=T,sep="\t",row.names = 1);dim(Moshepost)
colnames(Moshepost) <- paste0("Moshepost_",colnames(Moshepost))
rownames(Moshepost) <- as.character(as.data.frame(strsplit(rownames(Moshepost),".predicted.to_"))[2,])
# #
# #
all <- merge(merge(Moshepre,Moshepost,by="row.names",all=T),Caushi,by.x="Row.names",by.y="row.names",all=T);dim(all)
all <- data.frame(all[,c("Row.names","Caushi_proportions","Caushi_odds.ratio","Caushi_NK_T_activation.Score_log2FC","Caushi_NK_T_activation.Score_p")],Caushi_NK_T_activation.Score_p_adj =p.adjust(all[,"Caushi_NK_T_activation.Score_p"],method="BH"),
                  all[,c("Moshepre_proportions","Moshepre_odds.ratio","Moshepre_NK_T_activation.Score_log2FC","Moshepre_NK_T_activation.Score_p")],Moshepre_NK_T_activation.Score_p_adj =p.adjust(all[,"Moshepre_NK_T_activation.Score_p"],method="BH"),
                  all[,c("Moshepost_proportions","Moshepost_NK_T_activation.Score_log2FC","Moshepost_NK_T_activation.Score_p")],Moshepost_NK_T_activation.Score_p_adj =p.adjust(all[,"Moshepost_NK_T_activation.Score_p"],method="BH"))
write.table(all,file=paste0(workdir,file="4_immune_response_Moshepre_Moshepost_Caushi_ageAndgender_to_reference_p_20250506.txt"),quote=F,sep="\t",row.names = F)



############################################
#### 2025.07.05 Jing Yang 
#### dot plot for Baolinpost (no Baolinpre plot because no responders in Baolinpre) and Moshepre and Moshepost data
############################################
### for p values

Baolinpost <- read.table(paste0(workdir,"4_immune_response_Baolinpost_to_reference.txt"),header=T,sep="\t",row.names = 1);dim(Baolinpost)
colnames(Baolinpost) <- paste0("Baolinpost_",colnames(Baolinpost))
rownames(Baolinpost) <- as.character(as.data.frame(strsplit(rownames(Baolinpost),".predicted.to_"))[2,])
Moshepre <- read.table(paste0(workdir,"4_immune_response_Moshepre_to_reference.txt"),header=T,sep="\t",row.names = 1);dim(Moshepre)
colnames(Moshepre) <- paste0("Moshepre_",colnames(Moshepre))
rownames(Moshepre) <- as.character(as.data.frame(strsplit(rownames(Moshepre),".predicted.to_"))[2,])
Moshepost <- read.table(paste0(workdir,"4_immune_response_Moshepost_to_reference.txt"),header=T,sep="\t",row.names = 1);dim(Moshepost)
colnames(Moshepost) <- paste0("Moshepost_",colnames(Moshepost))
rownames(Moshepost) <- as.character(as.data.frame(strsplit(rownames(Moshepost),".predicted.to_"))[2,])
# #
# #
all <- merge(merge(Baolinpost,Moshepre,by="row.names",all=T),Moshepost,by.x="Row.names",by.y="row.names",all=T);dim(all)
all <- data.frame(all[,c("Row.names","Baolinpost_proportions","Baolinpost_odds.ratio","Baolinpost_NK_T_activation.Score_log2FC","Baolinpost_NK_T_activation.Score_p")],Baolinpost_NK_T_activation.Score_p_adj =p.adjust(all[,"Baolinpost_NK_T_activation.Score_p"],method="BH"),
                  all[,c("Moshepre_proportions","Moshepre_odds.ratio","Moshepre_NK_T_activation.Score_log2FC","Moshepre_NK_T_activation.Score_p")],Moshepre_NK_T_activation.Score_p_adj =p.adjust(all[,"Moshepre_NK_T_activation.Score_p"],method="BH"),
                  all[,c("Moshepost_proportions","Moshepost_NK_T_activation.Score_log2FC","Moshepost_NK_T_activation.Score_p")],Moshepost_NK_T_activation.Score_p_adj =p.adjust(all[,"Moshepost_NK_T_activation.Score_p"],method="BH"))
write.table(all,file=paste0(workdir,file="4_immune_response_Baolinpost_Moshepre_Moshepost_to_reference_p_20250506.txt"),quote=F,sep="\t",row.names = F)


library(Seurat)
library(SeuratObject)
library(mclust)
library(cluster)
library(ggplot2)
library(plyr)
library(ggpubr)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"


# SampleInfo_Kathryn <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Clonal replacement of tumor-specific T cells following PD-1 blockade/clinical_meta.txt",header=T,sep="\t")
SampleInfo_Baolin <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Temporal single-cell tracing reveals clonal revival and expansion of precursor exhausted T cells during anti-PD-1 therapy in lung cancer/GSE179994_Tcell.metadata_with_annotation.txt",header=T,sep="\t");dim(SampleInfo_Baolin)

Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
Baolin_T_NK_cells1 <- subset(Baolin_T_NK_cells1,subset=prepost=="post")
activated_T_proliferation <- c("ABL1","AGER","ARG1","BTN2A2","BTN3A1","CADM1","CASP3","CD24","CD274","CLC","CRTAM","EPO","FADD","FOXP3","FYN","GPAM","HHLA2","HMGB1","ICOSLG","IDO1","IGF1","IGF2","IGFBP2","IL12B","IL12RB1","IL18","IL2","IL23A","IL23R","IL27RA","IL2RA","JAK3","LGALS9","LILRB4","LRRC32","MIR181C","MIR21","MIR30B","PDCD1LG2","PPP3CA","PRKAR1A","PRNP","PYCARD","RC3H1","RIPK3","RPS3","SCRIB","SLAMF1","STAT5B","TMIGD2","TNFSF9")
NK_T_activation <- c("CD300A","ELF4","HSPH1","IL12A","IL12B","IL15","IL18","IL23A","IL23R","JAK2","RASAL3","TYK2","ZBTB7B")
NK_T_proliferation <- c("ELF4","IL12B","IL15","IL18","IL23A","JAK2","RASAL3","TYK2","ZBTB7B")
NK_T_differentiation <- c("AP3B1","AP3D1","ATF2","ITK","PRDM1","TGFBR2","TOX","ZBTB16","ZBTB7B","ZNF683")

Baolin_T_NK_cells1 <- CellCycleScoring(Baolin_T_NK_cells1,s.features = activated_T_proliferation, g2m.features = NK_T_activation)
names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="S.Score")] <- "activated_T_proliferation.Score"
names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_activation.Score"
Baolin_T_NK_cells1 <- CellCycleScoring(Baolin_T_NK_cells1,s.features = NK_T_proliferation, g2m.features = NK_T_differentiation)
names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
Baolin_T_NK_cells1 <- CellCycleScoring(Baolin_T_NK_cells1,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
# names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
# names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"

Baolin_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Baolin_T_NK_cells2 <- subset(Baolin_T_NK_cells2,subset=prepost=="post")
Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4 <- Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
rm(Baolin_T_NK_cells2)
Baolin_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Baolin_T_NK_cells3 <- subset(Baolin_T_NK_cells3,subset=prepost=="post")
Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells <- Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
rm(Baolin_T_NK_cells3)
Baolin_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Baolin_T_NK_cells4 <- subset(Baolin_T_NK_cells4,subset=prepost=="post")
Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells <- Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells
rm(Baolin_T_NK_cells4)
Baolin_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"))
Baolin_T_NK_cells5 <- subset(Baolin_T_NK_cells5,subset=prepost=="post")
Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells <- Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells
rm(Baolin_T_NK_cells5)
Baolin_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"))
Baolin_T_NK_cells6 <- subset(Baolin_T_NK_cells6,subset=prepost=="post")
Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells <- Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells
rm(Baolin_T_NK_cells6)


Baolin_pre_output <- data.frame(cells = names(Baolin_T_NK_cells1$SampleId),
                                Baolin_T_NK_cells1$SampleId,
                                # Baolin_T_NK_cells1$condition,
                                Baolin_T_NK_cells1$RNA_snn_res.0.1,
                                Baolin_T_NK_cells1$RNA_snn_res.0.15,
                                Baolin_T_NK_cells1$RNA_snn_res.0.2,
                                Baolin_T_NK_cells1$RNA_snn_res.0.21,
                                Baolin_T_NK_cells1$RNA_snn_res.0.22,
                                Baolin_T_NK_cells1$RNA_snn_res.0.3,
                                Baolin_T_NK_cells1$RNA_snn_res.0.4,
                                Baolin_T_NK_cells1$RNA_snn_res.0.5,
                                Baolin_T_NK_cells1$RNA_snn_res.0.6,
                                Baolin_T_NK_cells1$RNA_snn_res.0.7,
                                Baolin_T_NK_cells1$RNA_snn_res.0.8,
                                Baolin_T_NK_cells1$RNA_snn_res.0.9,
                                Baolin_T_NK_cells1$RNA_snn_res.1,
                                Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,
                                Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4,
                                Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells,
                                Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells,
                                Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells,
                                Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells,
                                Baolin_T_NK_cells1$G2M.Score,
                                Baolin_T_NK_cells1$activated_T_proliferation.Score,
                                Baolin_T_NK_cells1$NK_T_activation.Score,
                                Baolin_T_NK_cells1$NK_T_differentiation.Score,
                                Baolin_T_NK_cells1$NK_T_proliferation.Score);dim(Baolin_pre_output)
Baolin_pre_output <- merge(SampleInfo_Baolin[,c("SampleId","condition")],Baolin_pre_output,by.x="SampleId",by.y="Baolin_T_NK_cells1.SampleId");dim(Baolin_pre_output)

Baolin_pre_output <- Baolin_pre_output[!is.na(Baolin_pre_output$condition),];dim(Baolin_pre_output)
col_names <- colnames(Baolin_pre_output)
Baolinpost_output <- NULL
for (eachcol in 4:22) { #colnames(Baolin_pre_output)[4:22]
  tmp_proportion <- as.matrix(table(Baolin_pre_output[,c("condition",col_names[eachcol])]))
  tmp_proportion_sum <- apply(tmp_proportion,1,sum)
  tmp_proportion <- as.data.frame(tmp_proportion/tmp_proportion_sum)
  tmp_proportion_NK_T_activation_socre <- data.frame(tmp_proportion,NK_T_activation.Score=NA)
  for (i in 1:nrow(tmp_proportion_NK_T_activation_socre)) {
    tmp_proportion_NK_T_activation_socre[i,4] <- mean(Baolin_pre_output[Baolin_pre_output$condition == tmp_proportion_NK_T_activation_socre[i,1] & Baolin_pre_output[,col_names[eachcol]] == tmp_proportion_NK_T_activation_socre[i,2],"Baolin_T_NK_cells1.NK_T_activation.Score"])
  }
  tmp_proportion_NK_T_activation_socre$condition <- paste0(col_names[eachcol],"post_",tmp_proportion_NK_T_activation_socre$condition)
  colnames(tmp_proportion_NK_T_activation_socre)[1:2] <- c("res_condition","clusters")
  Baolinpost_output <- rbind(Baolinpost_output,tmp_proportion_NK_T_activation_socre)
}
rm(Baolin_T_NK_cells1)




SampleInfo_Moshe <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From Board/Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma/cells_all_scp.txt",header=T,sep="\t")
SampleInfo_Moshe <- SampleInfo_Moshe[-1,]
SampleInfo_Moshe <- unique(SampleInfo_Moshe[,c(2:3,8,10:13)])
SampleInfo_Moshe <- SampleInfo_Moshe[SampleInfo_Moshe$response %in% c("R","NR"),]
SampleInfo_Moshe <- data.frame(SampleId = SampleInfo_Moshe$patient,condition=SampleInfo_Moshe$response,treatment=SampleInfo_Moshe$therapy,prepost=SampleInfo_Moshe$prepost,gender=SampleInfo_Moshe$gender,age=SampleInfo_Moshe$age,os=SampleInfo_Moshe$survival_days,stringsAsFactors=F)
SampleInfo_Moshe$prepost <- tolower(SampleInfo_Moshe$prepost)
SampleInfo_Moshe <- SampleInfo_Moshe[c(1:10,12:48),]
SampleInfo_Moshe <- subset(SampleInfo_Moshe,prepost=="pre")
# SampleInfo_Caushi <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Transcriptional programs of neoantigen-specific TIL in anti-PD-1-treated lung cancers/clinical_meta4_simplified.txt",header=T,sep="\t")

Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
Moshe_T_NK_cells1 <- subset(Moshe_T_NK_cells1,subset=prepost=="Pre")
# Naive_genes <- c("TCF7", "SELL", "LEF1", "CCR7", "IL7R", "CD27", "CD28", "S1PR1")
# cytotoxic_genes <- c("CX3CR1", "PRF1", "GZMB", "GMZH", "GNLY", "FGFBP2", "FCGR3A", "KLRG1", "LYAR", "GZMK", "GZMM", "TXNIP", "FCRL6", "NKG7", "CCL4", "CST7", "GZMA", "IFNG", "CCL3", "KLRD1", "GZMH", "TBX21", "EOMES", "S1PR1", "S1PR5", "CXCR4", "CXCR3", "CD44")
# preexhausted_genes <- c("GZMK","PDCD1","ZNF683","ITGAE","CD28")
# rm_genes <- c("ITGAE", "HAVCR2", "GZMA", "GZMB", "IFNG", "ENTPD1", "CXCL13", "TNFRSF9", "PDCD1", "CCL3", "CTLA4", "TIGIT", "LAG3", "PRF1", "CD6", "CXCL1", "XCL2", "MYADM", "CAPG", "RORA", "NR4A1", "NR4A2", "NR4A3", "CD69", "CD39", "MKI67", "TOP2A", "CCNA2", "KIF2C ", "HMGB2", "TUBA1B", "TUBB", "H2AFZ", "CKS1B", "STMN1")
# exhausted_genes <- c("LAYN", "ITGAE", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "CXCL13", "CD38", "ENTPD1", "TOX", "IFNG", "GZMB", "MIR155HG", "TNFRSF9", "CDK1", "CCNB1", "MKI67", "CDK4", "RB1", "HSPH1", "HSPB1")
# memory_genes <- c("TCF7 ", "CCR7", "SELL ", "IL7R", "SELL", "LTB", "LEF1", "EOMES", "GZMK", "CXCR5", "GPR183", "CD27", "CD28", "GZMA", "CCL5", "S1PR1", "MYADM", "VIM", "ANKRD28", "ATP2B1") 
activated_T_proliferation <- c("ABL1","AGER","ARG1","BTN2A2","BTN3A1","CADM1","CASP3","CD24","CD274","CLC","CRTAM","EPO","FADD","FOXP3","FYN","GPAM","HHLA2","HMGB1","ICOSLG","IDO1","IGF1","IGF2","IGFBP2","IL12B","IL12RB1","IL18","IL2","IL23A","IL23R","IL27RA","IL2RA","JAK3","LGALS9","LILRB4","LRRC32","MIR181C","MIR21","MIR30B","PDCD1LG2","PPP3CA","PRKAR1A","PRNP","PYCARD","RC3H1","RIPK3","RPS3","SCRIB","SLAMF1","STAT5B","TMIGD2","TNFSF9")
NK_T_activation <- c("CD300A","ELF4","HSPH1","IL12A","IL12B","IL15","IL18","IL23A","IL23R","JAK2","RASAL3","TYK2","ZBTB7B")
NK_T_proliferation <- c("ELF4","IL12B","IL15","IL18","IL23A","JAK2","RASAL3","TYK2","ZBTB7B")
NK_T_differentiation <- c("AP3B1","AP3D1","ATF2","ITK","PRDM1","TGFBR2","TOX","ZBTB16","ZBTB7B","ZNF683")

# Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = Naive_genes, g2m.features = cytotoxic_genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "Naive.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "Cytotoxic.Score"
# Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = preexhausted_genes, g2m.features = rm_genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "Preexhausted.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "Rm.Score"
# Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = exhausted_genes, g2m.features = memory_genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "Exhausted.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "Memory.Score"
Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = activated_T_proliferation, g2m.features = NK_T_activation)
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "activated_T_proliferation.Score"
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_activation.Score"
Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = NK_T_proliferation, g2m.features = NK_T_differentiation)
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"

Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Moshe_T_NK_cells2 <- subset(Moshe_T_NK_cells2,subset=prepost=="Pre")
Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4 <- Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
rm(Moshe_T_NK_cells2)
Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Moshe_T_NK_cells3 <- subset(Moshe_T_NK_cells3,subset=prepost=="Pre")
Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells <- Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
rm(Moshe_T_NK_cells3)
Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Moshe_T_NK_cells4 <- subset(Moshe_T_NK_cells4,subset=prepost=="Pre")
Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells <- Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells
rm(Moshe_T_NK_cells4)
Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
Moshe_T_NK_cells5 <- subset(Moshe_T_NK_cells5,subset=prepost=="Pre")
Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells <- Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells
rm(Moshe_T_NK_cells5)
Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
Moshe_T_NK_cells6 <- subset(Moshe_T_NK_cells6,subset=prepost=="Pre")
Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells <- Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells
rm(Moshe_T_NK_cells6)


Moshe_pre_output <- data.frame(cells = names(Moshe_T_NK_cells1$SampleId),
                               Moshe_T_NK_cells1$SampleId,
                               # Moshe_T_NK_cells1$condition,
                               Moshe_T_NK_cells1$RNA_snn_res.0.1,
                               Moshe_T_NK_cells1$RNA_snn_res.0.2,
                               Moshe_T_NK_cells1$RNA_snn_res.0.3,
                               Moshe_T_NK_cells1$RNA_snn_res.0.35,
                               Moshe_T_NK_cells1$RNA_snn_res.0.4,
                               Moshe_T_NK_cells1$RNA_snn_res.0.5,
                               Moshe_T_NK_cells1$RNA_snn_res.0.52,
                               Moshe_T_NK_cells1$RNA_snn_res.0.6,
                               Moshe_T_NK_cells1$RNA_snn_res.0.7,
                               Moshe_T_NK_cells1$RNA_snn_res.0.8,
                               Moshe_T_NK_cells1$RNA_snn_res.0.82,
                               Moshe_T_NK_cells1$RNA_snn_res.0.8245,
                               Moshe_T_NK_cells1$RNA_snn_res.0.9,
                               Moshe_T_NK_cells1$RNA_snn_res.1,
                               Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,
                               Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4,
                               Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells,
                               Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells,
                               Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells,
                               Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells,
                               Moshe_T_NK_cells1$G2M.Score,
                               Moshe_T_NK_cells1$activated_T_proliferation.Score,
                               Moshe_T_NK_cells1$NK_T_activation.Score,
                               Moshe_T_NK_cells1$NK_T_differentiation.Score,
                               Moshe_T_NK_cells1$NK_T_proliferation.Score);dim(Moshe_pre_output)
Moshe_pre_output <- merge(SampleInfo_Moshe[,c("SampleId","condition")],Moshe_pre_output,by.x="SampleId",by.y="Moshe_T_NK_cells1.SampleId");dim(Moshe_pre_output)

Moshe_pre_output <- Moshe_pre_output[!is.na(Moshe_pre_output$condition),];dim(Moshe_pre_output)
col_names <- colnames(Moshe_pre_output)
Moshepre_output <- NULL
for (eachcol in 4:23) { #colnames(Moshe_pre_output)[4:23]
  tmp_proportion <- as.matrix(table(Moshe_pre_output[,c("condition",col_names[eachcol])]))
  tmp_proportion_sum <- apply(tmp_proportion,1,sum)
  tmp_proportion <- as.data.frame(tmp_proportion/tmp_proportion_sum)
  tmp_proportion_NK_T_activation_socre <- data.frame(tmp_proportion,NK_T_activation.Score=NA)
  for (i in 1:nrow(tmp_proportion_NK_T_activation_socre)) {
    tmp_proportion_NK_T_activation_socre[i,4] <- mean(Moshe_pre_output[Moshe_pre_output$condition == tmp_proportion_NK_T_activation_socre[i,1] & Moshe_pre_output[,col_names[eachcol]] == tmp_proportion_NK_T_activation_socre[i,2],"Moshe_T_NK_cells1.NK_T_activation.Score"])
  }
  tmp_proportion_NK_T_activation_socre$condition <- paste0(col_names[eachcol],"pre_",tmp_proportion_NK_T_activation_socre$condition)
  colnames(tmp_proportion_NK_T_activation_socre)[1:2] <- c("res_condition","clusters")
  Moshepre_output <- rbind(Moshepre_output,tmp_proportion_NK_T_activation_socre)
}
rm(Moshe_T_NK_cells1)


SampleInfo_Moshe <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From Board/Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma/cells_all_scp.txt",header=T,sep="\t")
SampleInfo_Moshe <- SampleInfo_Moshe[-1,]
SampleInfo_Moshe <- unique(SampleInfo_Moshe[,c(2:3,8,10:13)])
SampleInfo_Moshe <- SampleInfo_Moshe[SampleInfo_Moshe$response %in% c("R","NR"),]
SampleInfo_Moshe <- data.frame(SampleId = SampleInfo_Moshe$patient,condition=SampleInfo_Moshe$response,treatment=SampleInfo_Moshe$therapy,prepost=SampleInfo_Moshe$prepost,gender=SampleInfo_Moshe$gender,age=SampleInfo_Moshe$age,os=SampleInfo_Moshe$survival_days,stringsAsFactors=F)
SampleInfo_Moshe$prepost <- tolower(SampleInfo_Moshe$prepost)
SampleInfo_Moshe <- SampleInfo_Moshe[c(1:10,12:48),]
SampleInfo_Moshe <- subset(SampleInfo_Moshe,prepost=="post")
# SampleInfo_Caushi <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Transcriptional programs of neoantigen-specific TIL in anti-PD-1-treated lung cancers/clinical_meta4_simplified.txt",header=T,sep="\t")

Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
Moshe_T_NK_cells1 <- subset(Moshe_T_NK_cells1,subset=prepost=="Post")
# Naive_genes <- c("TCF7", "SELL", "LEF1", "CCR7", "IL7R", "CD27", "CD28", "S1PR1")
# cytotoxic_genes <- c("CX3CR1", "PRF1", "GZMB", "GMZH", "GNLY", "FGFBP2", "FCGR3A", "KLRG1", "LYAR", "GZMK", "GZMM", "TXNIP", "FCRL6", "NKG7", "CCL4", "CST7", "GZMA", "IFNG", "CCL3", "KLRD1", "GZMH", "TBX21", "EOMES", "S1PR1", "S1PR5", "CXCR4", "CXCR3", "CD44")
# preexhausted_genes <- c("GZMK","PDCD1","ZNF683","ITGAE","CD28")
# rm_genes <- c("ITGAE", "HAVCR2", "GZMA", "GZMB", "IFNG", "ENTPD1", "CXCL13", "TNFRSF9", "PDCD1", "CCL3", "CTLA4", "TIGIT", "LAG3", "PRF1", "CD6", "CXCL1", "XCL2", "MYADM", "CAPG", "RORA", "NR4A1", "NR4A2", "NR4A3", "CD69", "CD39", "MKI67", "TOP2A", "CCNA2", "KIF2C ", "HMGB2", "TUBA1B", "TUBB", "H2AFZ", "CKS1B", "STMN1")
# exhausted_genes <- c("LAYN", "ITGAE", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "CXCL13", "CD38", "ENTPD1", "TOX", "IFNG", "GZMB", "MIR155HG", "TNFRSF9", "CDK1", "CCNB1", "MKI67", "CDK4", "RB1", "HSPH1", "HSPB1")
# memory_genes <- c("TCF7 ", "CCR7", "SELL ", "IL7R", "SELL", "LTB", "LEF1", "EOMES", "GZMK", "CXCR5", "GPR183", "CD27", "CD28", "GZMA", "CCL5", "S1PR1", "MYADM", "VIM", "ANKRD28", "ATP2B1") 
activated_T_proliferation <- c("ABL1","AGER","ARG1","BTN2A2","BTN3A1","CADM1","CASP3","CD24","CD274","CLC","CRTAM","EPO","FADD","FOXP3","FYN","GPAM","HHLA2","HMGB1","ICOSLG","IDO1","IGF1","IGF2","IGFBP2","IL12B","IL12RB1","IL18","IL2","IL23A","IL23R","IL27RA","IL2RA","JAK3","LGALS9","LILRB4","LRRC32","MIR181C","MIR21","MIR30B","PDCD1LG2","PPP3CA","PRKAR1A","PRNP","PYCARD","RC3H1","RIPK3","RPS3","SCRIB","SLAMF1","STAT5B","TMIGD2","TNFSF9")
NK_T_activation <- c("CD300A","ELF4","HSPH1","IL12A","IL12B","IL15","IL18","IL23A","IL23R","JAK2","RASAL3","TYK2","ZBTB7B")
NK_T_proliferation <- c("ELF4","IL12B","IL15","IL18","IL23A","JAK2","RASAL3","TYK2","ZBTB7B")
NK_T_differentiation <- c("AP3B1","AP3D1","ATF2","ITK","PRDM1","TGFBR2","TOX","ZBTB16","ZBTB7B","ZNF683")

# Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = Naive_genes, g2m.features = cytotoxic_genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "Naive.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "Cytotoxic.Score"
# Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = preexhausted_genes, g2m.features = rm_genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "Preexhausted.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "Rm.Score"
# Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = exhausted_genes, g2m.features = memory_genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "Exhausted.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "Memory.Score"
Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = activated_T_proliferation, g2m.features = NK_T_activation)
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "activated_T_proliferation.Score"
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_activation.Score"
Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = NK_T_proliferation, g2m.features = NK_T_differentiation)
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
# names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"

Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
Moshe_T_NK_cells2 <- subset(Moshe_T_NK_cells2,subset=prepost=="Post")
Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4 <- Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
rm(Moshe_T_NK_cells2)
Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
Moshe_T_NK_cells3 <- subset(Moshe_T_NK_cells3,subset=prepost=="Post")
Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells <- Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
rm(Moshe_T_NK_cells3)
Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
Moshe_T_NK_cells4 <- subset(Moshe_T_NK_cells4,subset=prepost=="Post")
Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells <- Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells
rm(Moshe_T_NK_cells4)
Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
Moshe_T_NK_cells5 <- subset(Moshe_T_NK_cells5,subset=prepost=="Post")
Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells <- Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells
rm(Moshe_T_NK_cells5)
Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
Moshe_T_NK_cells6 <- subset(Moshe_T_NK_cells6,subset=prepost=="Post")
Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells <- Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells
rm(Moshe_T_NK_cells6)


Moshe_pre_output <- data.frame(cells = names(Moshe_T_NK_cells1$SampleId),
                               Moshe_T_NK_cells1$SampleId,
                               # Moshe_T_NK_cells1$condition,
                               Moshe_T_NK_cells1$RNA_snn_res.0.1,
                               Moshe_T_NK_cells1$RNA_snn_res.0.2,
                               Moshe_T_NK_cells1$RNA_snn_res.0.3,
                               Moshe_T_NK_cells1$RNA_snn_res.0.35,
                               Moshe_T_NK_cells1$RNA_snn_res.0.4,
                               Moshe_T_NK_cells1$RNA_snn_res.0.5,
                               Moshe_T_NK_cells1$RNA_snn_res.0.52,
                               Moshe_T_NK_cells1$RNA_snn_res.0.6,
                               Moshe_T_NK_cells1$RNA_snn_res.0.7,
                               Moshe_T_NK_cells1$RNA_snn_res.0.8,
                               Moshe_T_NK_cells1$RNA_snn_res.0.82,
                               Moshe_T_NK_cells1$RNA_snn_res.0.8245,
                               Moshe_T_NK_cells1$RNA_snn_res.0.9,
                               Moshe_T_NK_cells1$RNA_snn_res.1,
                               Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,
                               Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4,
                               Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells,
                               Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells,
                               Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells,
                               Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells,
                               Moshe_T_NK_cells1$G2M.Score,
                               Moshe_T_NK_cells1$activated_T_proliferation.Score,
                               Moshe_T_NK_cells1$NK_T_activation.Score,
                               Moshe_T_NK_cells1$NK_T_differentiation.Score,
                               Moshe_T_NK_cells1$NK_T_proliferation.Score);dim(Moshe_pre_output)
Moshe_pre_output <- merge(SampleInfo_Moshe[,c("SampleId","condition")],Moshe_pre_output,by.x="SampleId",by.y="Moshe_T_NK_cells1.SampleId");dim(Moshe_pre_output)

Moshe_pre_output <- Moshe_pre_output[!is.na(Moshe_pre_output$condition),];dim(Moshe_pre_output)
col_names <- colnames(Moshe_pre_output)
Moshepost_output <- NULL
for (eachcol in 4:23) { #colnames(Moshe_pre_output)[4:23]
  tmp_proportion <- as.matrix(table(Moshe_pre_output[,c("condition",col_names[eachcol])]))
  tmp_proportion_sum <- apply(tmp_proportion,1,sum)
  tmp_proportion <- as.data.frame(tmp_proportion/tmp_proportion_sum)
  tmp_proportion_NK_T_activation_socre <- data.frame(tmp_proportion,NK_T_activation.Score=NA)
  for (i in 1:nrow(tmp_proportion_NK_T_activation_socre)) {
    tmp_proportion_NK_T_activation_socre[i,4] <- mean(Moshe_pre_output[Moshe_pre_output$condition == tmp_proportion_NK_T_activation_socre[i,1] & Moshe_pre_output[,col_names[eachcol]] == tmp_proportion_NK_T_activation_socre[i,2],"Moshe_T_NK_cells1.NK_T_activation.Score"])
  }
  tmp_proportion_NK_T_activation_socre$condition <- paste0(col_names[eachcol],"post_",tmp_proportion_NK_T_activation_socre$condition)
  colnames(tmp_proportion_NK_T_activation_socre)[1:2] <- c("res_condition","clusters")
  Moshepost_output <- rbind(Moshepost_output,tmp_proportion_NK_T_activation_socre)
}
rm(Moshe_T_NK_cells1)

Baolinpost_Moshepre_Moshepost_output <- rbind(Baolinpost_output,Moshepre_output,Moshepost_output)
write.table(Baolinpost_Moshepre_Moshepost_output,file=paste0(workdir,"4_immune_response_Baolinpost_Moshepre_Moshepost_to_reference_plot.txt"),quote = F,sep="\t")

##### TICA plot
TICA_plot <- Baolinpost_Moshepre_Moshepost_output[Baolinpost_Moshepre_Moshepost_output$res_condition %in% 
                                                  c("Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_NR","Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_R",
                                                    "Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspre_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspre_R",
                                                    "Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_R"
                                                    ),]

Baolinpost_TICA_plot <- TICA_plot[TICA_plot$res_condition %in% 
                                c("Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_NR","Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_R"
                                ),]
if (length(unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Baolinpost_TICA_plot$clusters)])>0) {
  
  add_needed <- unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Baolinpost_TICA_plot$clusters)]
  add_needed <- data.frame(rep(c("Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_NR","Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Baolinpost_TICA_plot)
  Baolinpost_TICA_plot <- rbind(Baolinpost_TICA_plot,add_needed)
}

Moshepre_TICA_plot <- TICA_plot[TICA_plot$res_condition %in% 
                               c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspre_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspre_R"
                               ),]
if (length(unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Moshepre_TICA_plot$clusters)])>0) {
  
  add_needed <- unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Moshepre_TICA_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspre_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspre_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepre_TICA_plot)
  Moshepre_TICA_plot <- rbind(Moshepre_TICA_plot,add_needed)
}


Moshepost_TICA_plot <- TICA_plot[TICA_plot$res_condition %in% 
                                  c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_R"
                                  ),]
if (length(unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Moshepost_TICA_plot$clusters)])>0) {
  
  add_needed <- unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Moshepost_TICA_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepost_TICA_plot)
  Moshepost_TICA_plot <- rbind(Moshepost_TICA_plot,add_needed)
}



# pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_TICA_dotplot.pdf"),width = 4.88, height = 5.49)
p1 <- ggplot(data=Baolinpost_TICA_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepre_TICA_dotplot.pdf"),width = 4.88, height = 5.49)
p2 <- ggplot(data=Moshepre_TICA_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepost_TICA_dotplot.pdf"),width = 4.88, height = 5.49)
p3 <- ggplot(data=Moshepost_TICA_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_Moshepre_Moshepost_TICA_dotplot.pdf"),width = 3.33, height = 5.49)
CombinePlots(list(p1,p2,p3),ncol=3)
dev.off()

##### pancancer_blueprint plot
pancancer_blueprint_plot <- Baolinpost_Moshepre_Moshepost_output[Baolinpost_Moshepre_Moshepost_output$res_condition %in% 
                                                                   c("Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_NR","Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_R",
                                                                     "Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspre_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspre_R",
                                                                     "Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_R"
                                                                   ),]

Baolinpost_pancancer_blueprint_plot <- pancancer_blueprint_plot[pancancer_blueprint_plot$res_condition %in% 
                                                                  c("Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_NR","Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_R"
                                                                  ),]
if (length(unique(pancancer_blueprint_plot$clusters)[!(unique(pancancer_blueprint_plot$clusters) %in% Baolinpost_pancancer_blueprint_plot$clusters)])>0) {
  
  add_needed <- unique(pancancer_blueprint_plot$clusters)[!(unique(pancancer_blueprint_plot$clusters) %in% Baolinpost_pancancer_blueprint_plot$clusters)]
  add_needed <- data.frame(rep(c("Baolin_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cellspost_NR","Baolin_T_NK_cells3.predicted.to_pancancer_blueprint_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Baolinpost_pancancer_blueprint_plot)
  Baolinpost_pancancer_blueprint_plot <- rbind(Baolinpost_pancancer_blueprint_plot,add_needed)
}

Moshepre_pancancer_blueprint_plot <- pancancer_blueprint_plot[pancancer_blueprint_plot$res_condition %in% 
                                                                c("Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspre_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspre_R"
                                                                ),]
if (length(unique(pancancer_blueprint_plot$clusters)[!(unique(pancancer_blueprint_plot$clusters) %in% Moshepre_pancancer_blueprint_plot$clusters)])>0) {
  
  add_needed <- unique(pancancer_blueprint_plot$clusters)[!(unique(pancancer_blueprint_plot$clusters) %in% Moshepre_pancancer_blueprint_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspre_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspre_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepre_pancancer_blueprint_plot)
  Moshepre_pancancer_blueprint_plot <- rbind(Moshepre_pancancer_blueprint_plot,add_needed)
}


Moshepost_pancancer_blueprint_plot <- pancancer_blueprint_plot[pancancer_blueprint_plot$res_condition %in% 
                                                                 c("Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_R"
                                                                 ),]
if (length(unique(pancancer_blueprint_plot$clusters)[!(unique(pancancer_blueprint_plot$clusters) %in% Moshepost_pancancer_blueprint_plot$clusters)])>0) {
  
  add_needed <- unique(pancancer_blueprint_plot$clusters)[!(unique(pancancer_blueprint_plot$clusters) %in% Moshepost_pancancer_blueprint_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepost_pancancer_blueprint_plot)
  Moshepost_pancancer_blueprint_plot <- rbind(Moshepost_pancancer_blueprint_plot,add_needed)
}



# pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_pancancer_blueprint_dotplot.pdf"),width = 4.1, height = 5.49)
p1 <- ggplot(data=Baolinpost_pancancer_blueprint_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepre_pancancer_blueprint_dotplot.pdf"),width = 4.1, height = 5.49)
p2 <- ggplot(data=Moshepre_pancancer_blueprint_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepost_pancancer_blueprint_dotplot.pdf"),width = 4.1, height = 5.49)
p3 <- ggplot(data=Moshepost_pancancer_blueprint_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_Moshepre_Moshepost_pancancer_blueprint_dotplot.pdf"),width = 3.33, height = 5.49)
CombinePlots(list(p1,p2,p3),ncol=3)
dev.off()

##### ABTC plot
ABTC_plot <- Baolinpost_Moshepre_Moshepost_output[Baolinpost_Moshepre_Moshepost_output$res_condition %in% 
                                                    c("Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_NR","Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_R",
                                                      "Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspre_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspre_R",
                                                      "Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_R"
                                                    ),]

Baolinpost_ABTC_plot <- ABTC_plot[ABTC_plot$res_condition %in% 
                                    c("Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_NR","Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_R"
                                    ),]
if (length(unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Baolinpost_ABTC_plot$clusters)])>0) {
  
  add_needed <- unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Baolinpost_ABTC_plot$clusters)]
  add_needed <- data.frame(rep(c("Baolin_T_NK_cells4.predicted.to_ABTC_T_NK_cellspost_NR","Baolin_T_NK_cells4.predicted.to_ABTC_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Baolinpost_ABTC_plot)
  Baolinpost_ABTC_plot <- rbind(Baolinpost_ABTC_plot,add_needed)
}

Moshepre_ABTC_plot <- ABTC_plot[ABTC_plot$res_condition %in% 
                                  c("Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspre_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspre_R"
                                  ),]
if (length(unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Moshepre_ABTC_plot$clusters)])>0) {
  
  add_needed <- unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Moshepre_ABTC_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspre_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspre_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepre_ABTC_plot)
  Moshepre_ABTC_plot <- rbind(Moshepre_ABTC_plot,add_needed)
}


Moshepost_ABTC_plot <- ABTC_plot[ABTC_plot$res_condition %in% 
                                   c("Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_R"
                                   ),]
if (length(unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Moshepost_ABTC_plot$clusters)])>0) {
  
  add_needed <- unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Moshepost_ABTC_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepost_ABTC_plot)
  Moshepost_ABTC_plot <- rbind(Moshepost_ABTC_plot,add_needed)
}



p1 <- ggplot(data=Baolinpost_ABTC_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") + 
  theme_classic() + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                          axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p2 <- ggplot(data=Moshepre_ABTC_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p3 <- ggplot(data=Moshepost_ABTC_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_Moshepre_Moshepost_ABTC_dotplot.pdf"),width = 3.33, height = 5.49)
CombinePlots(list(p1,p2,p3),ncol=3)
dev.off()

##### HCC plot
HCC_plot <- Baolinpost_Moshepre_Moshepost_output[Baolinpost_Moshepre_Moshepost_output$res_condition %in% 
                                                   c("Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_NR","Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_R",
                                                     "Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspre_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspre_R",
                                                     "Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_R"
                                                   ),]

Baolinpost_HCC_plot <- HCC_plot[HCC_plot$res_condition %in% 
                                  c("Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_NR","Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_R"
                                  ),]
if (length(unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Baolinpost_HCC_plot$clusters)])>0) {
  
  add_needed <- unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Baolinpost_HCC_plot$clusters)]
  add_needed <- data.frame(rep(c("Baolin_T_NK_cells5.predicted.to_HCC_T_NK_cellspost_NR","Baolin_T_NK_cells5.predicted.to_HCC_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Baolinpost_HCC_plot)
  Baolinpost_HCC_plot <- rbind(Baolinpost_HCC_plot,add_needed)
}

Moshepre_HCC_plot <- HCC_plot[HCC_plot$res_condition %in% 
                                c("Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspre_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspre_R"
                                ),]
if (length(unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Moshepre_HCC_plot$clusters)])>0) {
  
  add_needed <- unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Moshepre_HCC_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspre_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspre_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepre_HCC_plot)
  Moshepre_HCC_plot <- rbind(Moshepre_HCC_plot,add_needed)
}


Moshepost_HCC_plot <- HCC_plot[HCC_plot$res_condition %in% 
                                 c("Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_R"
                                 ),]
if (length(unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Moshepost_HCC_plot$clusters)])>0) {
  
  add_needed <- unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Moshepost_HCC_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepost_HCC_plot)
  Moshepost_HCC_plot <- rbind(Moshepost_HCC_plot,add_needed)
}



# pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_HCC_dotplot.pdf"),width = 4.88, height = 5.49)
p1 <- ggplot(data=Baolinpost_HCC_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepre_HCC_dotplot.pdf"),width = 4.88, height = 5.49)
p2 <- ggplot(data=Moshepre_HCC_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepost_HCC_dotplot.pdf"),width = 4.88, height = 5.49)
p3 <- ggplot(data=Moshepost_HCC_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_Moshepre_Moshepost_HCC_dotplot.pdf"),width = 3.33, height = 5.49)
CombinePlots(list(p1,p2,p3),ncol=3)
dev.off()

##### STAD plot
STAD_plot <- Baolinpost_Moshepre_Moshepost_output[Baolinpost_Moshepre_Moshepost_output$res_condition %in% 
                                                    c("Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_NR","Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_R",
                                                      "Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspre_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspre_R",
                                                      "Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_R"
                                                    ),]

Baolinpost_STAD_plot <- STAD_plot[STAD_plot$res_condition %in% 
                                    c("Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_NR","Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_R"
                                    ),]
if (length(unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Baolinpost_STAD_plot$clusters)])>0) {
  
  add_needed <- unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Baolinpost_STAD_plot$clusters)]
  add_needed <- data.frame(rep(c("Baolin_T_NK_cells6.predicted.to_STAD_T_NK_cellspost_NR","Baolin_T_NK_cells6.predicted.to_STAD_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Baolinpost_STAD_plot)
  Baolinpost_STAD_plot <- rbind(Baolinpost_STAD_plot,add_needed)
}

Moshepre_STAD_plot <- STAD_plot[STAD_plot$res_condition %in% 
                                  c("Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspre_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspre_R"
                                  ),]
if (length(unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Moshepre_STAD_plot$clusters)])>0) {
  
  add_needed <- unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Moshepre_STAD_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspre_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspre_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepre_STAD_plot)
  Moshepre_STAD_plot <- rbind(Moshepre_STAD_plot,add_needed)
}


Moshepost_STAD_plot <- STAD_plot[STAD_plot$res_condition %in% 
                                   c("Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_R"
                                   ),]
if (length(unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Moshepost_STAD_plot$clusters)])>0) {
  
  add_needed <- unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Moshepost_STAD_plot$clusters)]
  add_needed <- data.frame(rep(c("Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cellspost_R"),length(add_needed)),
                           unlist(lapply(add_needed,function(x) rep(x,2))),
                           rep(NA,2*length(add_needed)),
                           rep(NA,2*length(add_needed)))
  colnames(add_needed) <- colnames(Moshepost_STAD_plot)
  Moshepost_STAD_plot <- rbind(Moshepost_STAD_plot,add_needed)
}



# pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_STAD_dotplot.pdf"),width = 4.88, height = 5.49)
p1 <- ggplot(data=Baolinpost_STAD_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepre_STAD_dotplot.pdf"),width = 4.88, height = 5.49)
p2 <- ggplot(data=Moshepre_STAD_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
# pdf(paste0(workdir,"4_immune_response_plots_Moshepost_STAD_dotplot.pdf"),width = 4.88, height = 5.49)
p3 <- ggplot(data=Moshepost_STAD_plot,
       mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
  geom_point(mapping=aes_string(size='Freq')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic()+ theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
                         axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# dev.off()
pdf(paste0(workdir,"4_immune_response_plots_Baolinpost_Moshepre_Moshepost_STAD_dotplot.pdf"),width = 3.33, height = 5.49)
CombinePlots(list(p1,p2,p3),ncol=3)
dev.off()


#### prepare plots 
  
  SampleInfo_Kathryn <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Clonal replacement of tumor-specific T cells following PD-1 blockade/clinical_meta.txt",header=T,sep="\t")
  
  Kathryn_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_TICA_T_NK_cells.rds"))
  activated_T_proliferation <- c("ABL1","AGER","ARG1","BTN2A2","BTN3A1","CADM1","CASP3","CD24","CD274","CLC","CRTAM","EPO","FADD","FOXP3","FYN","GPAM","HHLA2","HMGB1","ICOSLG","IDO1","IGF1","IGF2","IGFBP2","IL12B","IL12RB1","IL18","IL2","IL23A","IL23R","IL27RA","IL2RA","JAK3","LGALS9","LILRB4","LRRC32","MIR181C","MIR21","MIR30B","PDCD1LG2","PPP3CA","PRKAR1A","PRNP","PYCARD","RC3H1","RIPK3","RPS3","SCRIB","SLAMF1","STAT5B","TMIGD2","TNFSF9")
  NK_T_activation <- c("CD300A","ELF4","HSPH1","IL12A","IL12B","IL15","IL18","IL23A","IL23R","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_proliferation <- c("ELF4","IL12B","IL15","IL18","IL23A","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_differentiation <- c("AP3B1","AP3D1","ATF2","ITK","PRDM1","TGFBR2","TOX","ZBTB16","ZBTB7B","ZNF683")
  
  Kathryn_T_NK_cells1 <- CellCycleScoring(Kathryn_T_NK_cells1,s.features = activated_T_proliferation, g2m.features = NK_T_activation)
  names(Kathryn_T_NK_cells1@meta.data)[which(names(Kathryn_T_NK_cells1@meta.data)=="S.Score")] <- "activated_T_proliferation.Score"
  names(Kathryn_T_NK_cells1@meta.data)[which(names(Kathryn_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_activation.Score"
  Kathryn_T_NK_cells1 <- CellCycleScoring(Kathryn_T_NK_cells1,s.features = NK_T_proliferation, g2m.features = NK_T_differentiation)
  names(Kathryn_T_NK_cells1@meta.data)[which(names(Kathryn_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
  names(Kathryn_T_NK_cells1@meta.data)[which(names(Kathryn_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
  Kathryn_T_NK_cells1 <- CellCycleScoring(Kathryn_T_NK_cells1,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  # names(Kathryn_T_NK_cells1@meta.data)[which(names(Kathryn_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
  # names(Kathryn_T_NK_cells1@meta.data)[which(names(Kathryn_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
  
  Kathryn_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
  Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4 <- Kathryn_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
  rm(Kathryn_T_NK_cells2)
  Kathryn_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
  Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells <- Kathryn_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
  rm(Kathryn_T_NK_cells3)
  Kathryn_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_ABTC_T_NK_cells.rds"))
  Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells <- Kathryn_T_NK_cells4$predicted.to_ABTC_T_NK_cells
  rm(Kathryn_T_NK_cells4)
  Kathryn_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_HCC_T_NK_cells.rds"))
  Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells <- Kathryn_T_NK_cells5$predicted.to_HCC_T_NK_cells
  rm(Kathryn_T_NK_cells5)
  Kathryn_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Kathryn_T_NK_cells_to_STAD_T_NK_cells.rds"))
  Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells <- Kathryn_T_NK_cells6$predicted.to_STAD_T_NK_cells
  rm(Kathryn_T_NK_cells6)
  
  
  Kathryn_pre_output <- data.frame(cells = names(Kathryn_T_NK_cells1$SampleId),
                           Kathryn_T_NK_cells1$SampleId,
                           # Kathryn_T_NK_cells1$condition,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.1,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.15,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.16,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.2,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.217,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.22,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.25,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.28,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.3,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.4,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.5,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.6,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.7,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.8,
                           Kathryn_T_NK_cells1$RNA_snn_res.0.9,
                           Kathryn_T_NK_cells1$RNA_snn_res.1,
                           Kathryn_T_NK_cells1$predicted.to_TICA_T_NK_cells,
                           Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4,
                           Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells,
                           Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells,
                           Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells,
                           Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells,
                           Kathryn_T_NK_cells1$G2M.Score,
                           Kathryn_T_NK_cells1$activated_T_proliferation.Score,
                           Kathryn_T_NK_cells1$NK_T_activation.Score,
                           Kathryn_T_NK_cells1$NK_T_differentiation.Score,
                           Kathryn_T_NK_cells1$NK_T_proliferation.Score);dim(Kathryn_pre_output)
  Kathryn_pre_output <- merge(SampleInfo_Kathryn[,c("SampleId","condition")],Kathryn_pre_output,by.x="SampleId",by.y="Kathryn_T_NK_cells1.SampleId");dim(Kathryn_pre_output)
  
  Kathryn_pre_output <- Kathryn_pre_output[!is.na(Kathryn_pre_output$condition),];dim(Kathryn_pre_output)
  col_names <- colnames(Kathryn_pre_output)
  Kathryn_output <- NULL
  for (eachcol in 4:25) { #colnames(Kathryn_pre_output)[4:25]
    tmp_proportion <- as.matrix(table(Kathryn_pre_output[,c("condition",col_names[eachcol])]))
    tmp_proportion_sum <- apply(tmp_proportion,1,sum)
    tmp_proportion <- as.data.frame(tmp_proportion/tmp_proportion_sum)
    tmp_proportion_NK_T_activation_socre <- data.frame(tmp_proportion,NK_T_activation.Score=NA)
    for (i in 1:nrow(tmp_proportion_NK_T_activation_socre)) {
      tmp_proportion_NK_T_activation_socre[i,4] <- mean(Kathryn_pre_output[Kathryn_pre_output$condition == tmp_proportion_NK_T_activation_socre[i,1] & Kathryn_pre_output[,col_names[eachcol]] == tmp_proportion_NK_T_activation_socre[i,2],"Kathryn_T_NK_cells1.NK_T_activation.Score"])
    }
    tmp_proportion_NK_T_activation_socre$condition <- paste0(col_names[eachcol],"_",tmp_proportion_NK_T_activation_socre$condition)
    colnames(tmp_proportion_NK_T_activation_socre)[1:2] <- c("res_condition","clusters")
    Kathryn_output <- rbind(Kathryn_output,tmp_proportion_NK_T_activation_socre)
  }
  rm(Kathryn_T_NK_cells1)
  
  # SampleInfo_Kathryn <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Clonal replacement of tumor-specific T cells following PD-1 blockade/clinical_meta.txt",header=T,sep="\t")
  SampleInfo_Baolin <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Temporal single-cell tracing reveals clonal revival and expansion of precursor exhausted T cells during anti-PD-1 therapy in lung cancer/GSE179994_Tcell.metadata_with_annotation.txt",header=T,sep="\t");dim(SampleInfo_Baolin)
    
  Baolin_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_TICA_T_NK_cells.rds"))
    activated_T_proliferation <- c("ABL1","AGER","ARG1","BTN2A2","BTN3A1","CADM1","CASP3","CD24","CD274","CLC","CRTAM","EPO","FADD","FOXP3","FYN","GPAM","HHLA2","HMGB1","ICOSLG","IDO1","IGF1","IGF2","IGFBP2","IL12B","IL12RB1","IL18","IL2","IL23A","IL23R","IL27RA","IL2RA","JAK3","LGALS9","LILRB4","LRRC32","MIR181C","MIR21","MIR30B","PDCD1LG2","PPP3CA","PRKAR1A","PRNP","PYCARD","RC3H1","RIPK3","RPS3","SCRIB","SLAMF1","STAT5B","TMIGD2","TNFSF9")
  NK_T_activation <- c("CD300A","ELF4","HSPH1","IL12A","IL12B","IL15","IL18","IL23A","IL23R","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_proliferation <- c("ELF4","IL12B","IL15","IL18","IL23A","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_differentiation <- c("AP3B1","AP3D1","ATF2","ITK","PRDM1","TGFBR2","TOX","ZBTB16","ZBTB7B","ZNF683")
  
   Baolin_T_NK_cells1 <- CellCycleScoring(Baolin_T_NK_cells1,s.features = activated_T_proliferation, g2m.features = NK_T_activation)
  names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="S.Score")] <- "activated_T_proliferation.Score"
  names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_activation.Score"
  Baolin_T_NK_cells1 <- CellCycleScoring(Baolin_T_NK_cells1,s.features = NK_T_proliferation, g2m.features = NK_T_differentiation)
  names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
  names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
  Baolin_T_NK_cells1 <- CellCycleScoring(Baolin_T_NK_cells1,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  # names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
  # names(Baolin_T_NK_cells1@meta.data)[which(names(Baolin_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
  
  Baolin_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
  Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4 <- Baolin_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
  rm(Baolin_T_NK_cells2)
  Baolin_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
  Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells <- Baolin_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
  rm(Baolin_T_NK_cells3)
  Baolin_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_ABTC_T_NK_cells.rds"))
  Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells <- Baolin_T_NK_cells4$predicted.to_ABTC_T_NK_cells
  rm(Baolin_T_NK_cells4)
  Baolin_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_HCC_T_NK_cells.rds"))
  Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells <- Baolin_T_NK_cells5$predicted.to_HCC_T_NK_cells
  rm(Baolin_T_NK_cells5)
  Baolin_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Baolin_T_NK_cells_to_STAD_T_NK_cells.rds"))
  Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells <- Baolin_T_NK_cells6$predicted.to_STAD_T_NK_cells
  rm(Baolin_T_NK_cells6)
  
  
  Baolin_pre_output <- data.frame(cells = names(Baolin_T_NK_cells1$SampleId),
                                  Baolin_T_NK_cells1$SampleId,
                                  # Baolin_T_NK_cells1$condition,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.1,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.15,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.2,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.21,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.22,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.3,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.4,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.5,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.6,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.7,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.8,
                                  Baolin_T_NK_cells1$RNA_snn_res.0.9,
                                  Baolin_T_NK_cells1$RNA_snn_res.1,
                                  Baolin_T_NK_cells1$predicted.to_TICA_T_NK_cells,
                                  Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4,
                                  Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells,
                                  Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells,
                                  Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells,
                                  Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells,
                                  Baolin_T_NK_cells1$G2M.Score,
                                  Baolin_T_NK_cells1$activated_T_proliferation.Score,
                                  Baolin_T_NK_cells1$NK_T_activation.Score,
                                  Baolin_T_NK_cells1$NK_T_differentiation.Score,
                                  Baolin_T_NK_cells1$NK_T_proliferation.Score);dim(Baolin_pre_output)
  Baolin_pre_output <- merge(SampleInfo_Baolin[,c("SampleId","condition")],Baolin_pre_output,by.x="SampleId",by.y="Baolin_T_NK_cells1.SampleId");dim(Baolin_pre_output)
  
  Baolin_pre_output <- Baolin_pre_output[!is.na(Baolin_pre_output$condition),];dim(Baolin_pre_output)
  col_names <- colnames(Baolin_pre_output)
  Baolin_output <- NULL
  for (eachcol in 4:22) { #colnames(Baolin_pre_output)[4:22]
    tmp_proportion <- as.matrix(table(Baolin_pre_output[,c("condition",col_names[eachcol])]))
    tmp_proportion_sum <- apply(tmp_proportion,1,sum)
    tmp_proportion <- as.data.frame(tmp_proportion/tmp_proportion_sum)
    tmp_proportion_NK_T_activation_socre <- data.frame(tmp_proportion,NK_T_activation.Score=NA)
    for (i in 1:nrow(tmp_proportion_NK_T_activation_socre)) {
      tmp_proportion_NK_T_activation_socre[i,4] <- mean(Baolin_pre_output[Baolin_pre_output$condition == tmp_proportion_NK_T_activation_socre[i,1] & Baolin_pre_output[,col_names[eachcol]] == tmp_proportion_NK_T_activation_socre[i,2],"Baolin_T_NK_cells1.NK_T_activation.Score"])
    }
    tmp_proportion_NK_T_activation_socre$condition <- paste0(col_names[eachcol],"_",tmp_proportion_NK_T_activation_socre$condition)
    colnames(tmp_proportion_NK_T_activation_socre)[1:2] <- c("res_condition","clusters")
    Baolin_output <- rbind(Baolin_output,tmp_proportion_NK_T_activation_socre)
  }
  rm(Baolin_T_NK_cells1)
  
  # SampleInfo_Kathryn <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Clonal replacement of tumor-specific T cells following PD-1 blockade/clinical_meta.txt",header=T,sep="\t")
  # SampleInfo_Baolin <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Temporal single-cell tracing reveals clonal revival and expansion of precursor exhausted T cells during anti-PD-1 therapy in lung cancer/GSE179994_Tcell.metadata_with_annotation.txt",header=T,sep="\t");dim(SampleInfo_Moshe)
  SampleInfo_Moshe <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From Board/Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma/cells_all_scp.txt",header=T,sep="\t")
  SampleInfo_Moshe <- SampleInfo_Moshe[-1,]
  SampleInfo_Moshe <- unique(SampleInfo_Moshe[,c(2:3,8,10:13)])
  SampleInfo_Moshe <- SampleInfo_Moshe[SampleInfo_Moshe$response %in% c("R","NR"),]
  SampleInfo_Moshe <- data.frame(SampleId = SampleInfo_Moshe$patient,condition=SampleInfo_Moshe$response,treatment=SampleInfo_Moshe$therapy,prepost=SampleInfo_Moshe$prepost,gender=SampleInfo_Moshe$gender,age=SampleInfo_Moshe$age,os=SampleInfo_Moshe$survival_days,stringsAsFactors=F)
  SampleInfo_Moshe$prepost <- tolower(SampleInfo_Moshe$prepost)
  SampleInfo_Moshe <- SampleInfo_Moshe[c(1:10,12:48),]
  # SampleInfo_Caushi <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Transcriptional programs of neoantigen-specific TIL in anti-PD-1-treated lung cancers/clinical_meta4_simplified.txt",header=T,sep="\t")
  
  Moshe_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_TICA_T_NK_cells.rds"))
  activated_T_proliferation <- c("ABL1","AGER","ARG1","BTN2A2","BTN3A1","CADM1","CASP3","CD24","CD274","CLC","CRTAM","EPO","FADD","FOXP3","FYN","GPAM","HHLA2","HMGB1","ICOSLG","IDO1","IGF1","IGF2","IGFBP2","IL12B","IL12RB1","IL18","IL2","IL23A","IL23R","IL27RA","IL2RA","JAK3","LGALS9","LILRB4","LRRC32","MIR181C","MIR21","MIR30B","PDCD1LG2","PPP3CA","PRKAR1A","PRNP","PYCARD","RC3H1","RIPK3","RPS3","SCRIB","SLAMF1","STAT5B","TMIGD2","TNFSF9")
  NK_T_activation <- c("CD300A","ELF4","HSPH1","IL12A","IL12B","IL15","IL18","IL23A","IL23R","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_proliferation <- c("ELF4","IL12B","IL15","IL18","IL23A","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_differentiation <- c("AP3B1","AP3D1","ATF2","ITK","PRDM1","TGFBR2","TOX","ZBTB16","ZBTB7B","ZNF683")
  
  Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = activated_T_proliferation, g2m.features = NK_T_activation)
  names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "activated_T_proliferation.Score"
  names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_activation.Score"
  Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = NK_T_proliferation, g2m.features = NK_T_differentiation)
  names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
  names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
  Moshe_T_NK_cells1 <- CellCycleScoring(Moshe_T_NK_cells1,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  # names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
  # names(Moshe_T_NK_cells1@meta.data)[which(names(Moshe_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
  
  Moshe_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
  Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4 <- Moshe_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
  rm(Moshe_T_NK_cells2)
  Moshe_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
  Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells <- Moshe_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
  rm(Moshe_T_NK_cells3)
  Moshe_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_ABTC_T_NK_cells.rds"))
  Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells <- Moshe_T_NK_cells4$predicted.to_ABTC_T_NK_cells
  rm(Moshe_T_NK_cells4)
  Moshe_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_HCC_T_NK_cells.rds"))
  Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells <- Moshe_T_NK_cells5$predicted.to_HCC_T_NK_cells
  rm(Moshe_T_NK_cells5)
  Moshe_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Moshe_T_NK_cells_to_STAD_T_NK_cells.rds"))
  Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells <- Moshe_T_NK_cells6$predicted.to_STAD_T_NK_cells
  rm(Moshe_T_NK_cells6)
  
  
  Moshe_pre_output <- data.frame(cells = names(Moshe_T_NK_cells1$SampleId),
                                 Moshe_T_NK_cells1$SampleId,
                                 # Moshe_T_NK_cells1$condition,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.1,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.2,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.3,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.35,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.4,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.5,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.52,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.6,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.7,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.8,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.82,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.8245,
                                 Moshe_T_NK_cells1$RNA_snn_res.0.9,
                                 Moshe_T_NK_cells1$RNA_snn_res.1,
                                 Moshe_T_NK_cells1$predicted.to_TICA_T_NK_cells,
                                 Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4,
                                 Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells,
                                 Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells,
                                 Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells,
                                 Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells,
                                 Moshe_T_NK_cells1$G2M.Score,
                                 Moshe_T_NK_cells1$activated_T_proliferation.Score,
                                 Moshe_T_NK_cells1$NK_T_activation.Score,
                                 Moshe_T_NK_cells1$NK_T_differentiation.Score,
                                 Moshe_T_NK_cells1$NK_T_proliferation.Score);dim(Moshe_pre_output)
  Moshe_pre_output <- merge(SampleInfo_Moshe[,c("SampleId","condition")],Moshe_pre_output,by.x="SampleId",by.y="Moshe_T_NK_cells1.SampleId");dim(Moshe_pre_output)
  
  Moshe_pre_output <- Moshe_pre_output[!is.na(Moshe_pre_output$condition),];dim(Moshe_pre_output)
  col_names <- colnames(Moshe_pre_output)
  Moshe_output <- NULL
  for (eachcol in 4:23) { #colnames(Moshe_pre_output)[4:23]
    tmp_proportion <- as.matrix(table(Moshe_pre_output[,c("condition",col_names[eachcol])]))
    tmp_proportion_sum <- apply(tmp_proportion,1,sum)
    tmp_proportion <- as.data.frame(tmp_proportion/tmp_proportion_sum)
    tmp_proportion_NK_T_activation_socre <- data.frame(tmp_proportion,NK_T_activation.Score=NA)
    for (i in 1:nrow(tmp_proportion_NK_T_activation_socre)) {
      tmp_proportion_NK_T_activation_socre[i,4] <- mean(Moshe_pre_output[Moshe_pre_output$condition == tmp_proportion_NK_T_activation_socre[i,1] & Moshe_pre_output[,col_names[eachcol]] == tmp_proportion_NK_T_activation_socre[i,2],"Moshe_T_NK_cells1.NK_T_activation.Score"])
    }
    tmp_proportion_NK_T_activation_socre$condition <- paste0(col_names[eachcol],"_",tmp_proportion_NK_T_activation_socre$condition)
    colnames(tmp_proportion_NK_T_activation_socre)[1:2] <- c("res_condition","clusters")
    Moshe_output <- rbind(Moshe_output,tmp_proportion_NK_T_activation_socre)
  }
  rm(Moshe_T_NK_cells1)
  
  
  SampleInfo_Caushi <- read.table("/Users/jingyang/Dropbox/work/Single-cell meta-analysis/data sets/From GEO/Transcriptional programs of neoantigen-specific TIL in anti-PD-1-treated lung cancers/clinical_meta4_simplified.txt",header=T,sep="\t")
  
  Caushi_T_NK_cells1 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_TICA_T_NK_cells.rds"))
  activated_T_proliferation <- c("ABL1","AGER","ARG1","BTN2A2","BTN3A1","CADM1","CASP3","CD24","CD274","CLC","CRTAM","EPO","FADD","FOXP3","FYN","GPAM","HHLA2","HMGB1","ICOSLG","IDO1","IGF1","IGF2","IGFBP2","IL12B","IL12RB1","IL18","IL2","IL23A","IL23R","IL27RA","IL2RA","JAK3","LGALS9","LILRB4","LRRC32","MIR181C","MIR21","MIR30B","PDCD1LG2","PPP3CA","PRKAR1A","PRNP","PYCARD","RC3H1","RIPK3","RPS3","SCRIB","SLAMF1","STAT5B","TMIGD2","TNFSF9")
  NK_T_activation <- c("CD300A","ELF4","HSPH1","IL12A","IL12B","IL15","IL18","IL23A","IL23R","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_proliferation <- c("ELF4","IL12B","IL15","IL18","IL23A","JAK2","RASAL3","TYK2","ZBTB7B")
  NK_T_differentiation <- c("AP3B1","AP3D1","ATF2","ITK","PRDM1","TGFBR2","TOX","ZBTB16","ZBTB7B","ZNF683")
  
  Caushi_T_NK_cells1 <- CellCycleScoring(Caushi_T_NK_cells1,s.features = activated_T_proliferation, g2m.features = NK_T_activation)
  names(Caushi_T_NK_cells1@meta.data)[which(names(Caushi_T_NK_cells1@meta.data)=="S.Score")] <- "activated_T_proliferation.Score"
  names(Caushi_T_NK_cells1@meta.data)[which(names(Caushi_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_activation.Score"
  Caushi_T_NK_cells1 <- CellCycleScoring(Caushi_T_NK_cells1,s.features = NK_T_proliferation, g2m.features = NK_T_differentiation)
  names(Caushi_T_NK_cells1@meta.data)[which(names(Caushi_T_NK_cells1@meta.data)=="S.Score")] <- "NK_T_proliferation.Score"
  names(Caushi_T_NK_cells1@meta.data)[which(names(Caushi_T_NK_cells1@meta.data)=="G2M.Score")] <- "NK_T_differentiation.Score"
  Caushi_T_NK_cells1 <- CellCycleScoring(Caushi_T_NK_cells1,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  
  Caushi_T_NK_cells2 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_T_cells_CD8_CD4.rds"))
  Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4 <- Caushi_T_NK_cells2$predicted.to_pancancer_T_cells_CD8_CD4
  rm(Caushi_T_NK_cells2)
  Caushi_T_NK_cells3 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_pancancer_blueprint_T_NK_cells.rds"))
  Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells <- Caushi_T_NK_cells3$predicted.to_pancancer_blueprint_T_NK_cells
  rm(Caushi_T_NK_cells3)
  Caushi_T_NK_cells4 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_ABTC_T_NK_cells.rds"))
  Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells <- Caushi_T_NK_cells4$predicted.to_ABTC_T_NK_cells
  rm(Caushi_T_NK_cells4)
  Caushi_T_NK_cells5 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_HCC_T_NK_cells.rds"))
  Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells <- Caushi_T_NK_cells5$predicted.to_HCC_T_NK_cells
  rm(Caushi_T_NK_cells5)
  Caushi_T_NK_cells6 <- readRDS(paste0(workdir,"3_mapping_query_dataset_Caushi_T_NK_cells_to_STAD_T_NK_cells.rds"))
  Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells <- Caushi_T_NK_cells6$predicted.to_STAD_T_NK_cells
  rm(Caushi_T_NK_cells6)
  
  
  Caushi_pre_output <- data.frame(cells = names(Caushi_T_NK_cells1$SampleId),
                                  Caushi_T_NK_cells1$SampleId,
                                  # Caushi_T_NK_cells1$condition,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.1,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.2,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.208,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.21,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.3,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.31,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.35,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.4,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.5,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.6,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.7,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.8,
                                  Caushi_T_NK_cells1$RNA_snn_res.0.9,
                                  Caushi_T_NK_cells1$RNA_snn_res.1,
                                  Caushi_T_NK_cells1$predicted.to_TICA_T_NK_cells,
                                  Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4,
                                  Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells,
                                  Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells,
                                  Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells,
                                  Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells,
                                  Caushi_T_NK_cells1$G2M.Score,
                                  Caushi_T_NK_cells1$activated_T_proliferation.Score,
                                  Caushi_T_NK_cells1$NK_T_activation.Score,
                                  Caushi_T_NK_cells1$NK_T_differentiation.Score,
                                  Caushi_T_NK_cells1$NK_T_proliferation.Score);dim(Caushi_pre_output)
  Caushi_pre_output <- merge(SampleInfo_Caushi[,c("SampleId","condition")],Caushi_pre_output,by.x="SampleId",by.y="Caushi_T_NK_cells1.SampleId");dim(Caushi_pre_output)
  
  Caushi_pre_output <- Caushi_pre_output[!is.na(Caushi_pre_output$condition),];dim(Caushi_pre_output)
  col_names <- colnames(Caushi_pre_output)
  Caushi_output <- NULL
  for (eachcol in 4:23) { #colnames(Caushi_pre_output)[4:23]
    tmp_proportion <- as.matrix(table(Caushi_pre_output[,c("condition",col_names[eachcol])]))
    tmp_proportion_sum <- apply(tmp_proportion,1,sum)
    tmp_proportion <- as.data.frame(tmp_proportion/tmp_proportion_sum)
    tmp_proportion_NK_T_activation_socre <- data.frame(tmp_proportion,NK_T_activation.Score=NA)
    for (i in 1:nrow(tmp_proportion_NK_T_activation_socre)) {
      tmp_proportion_NK_T_activation_socre[i,4] <- mean(Caushi_pre_output[Caushi_pre_output$condition == tmp_proportion_NK_T_activation_socre[i,1] & Caushi_pre_output[,col_names[eachcol]] == tmp_proportion_NK_T_activation_socre[i,2],"Caushi_T_NK_cells1.NK_T_activation.Score"])
    }
    tmp_proportion_NK_T_activation_socre$condition <- paste0(col_names[eachcol],"_",tmp_proportion_NK_T_activation_socre$condition)
    colnames(tmp_proportion_NK_T_activation_socre)[1:2] <- c("res_condition","clusters")
    Caushi_output <- rbind(Caushi_output,tmp_proportion_NK_T_activation_socre)
  }
  rm(Caushi_T_NK_cells1)
  
  Kathryn_Baolin_Moshe_Caushi_output <- rbind(Kathryn_output,Baolin_output,Moshe_output,Caushi_output)
  write.table(Kathryn_Baolin_Moshe_Caushi_output,file=paste0(workdir,"4_immune_response_Kathryn_Baolin_Moshe_Caushi_to_reference_plot.txt"),quote = F,sep="\t")
  
  
##### TICA plot
  TICA_plot <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% 
                                                    c("Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells_R",
                                                      "Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells_R",
                                                      "Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells_R",
                                                      "Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells_non-MPR","Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells_MPR"
  ),]
  Kathryn_TICA_plot <- TICA_plot[TICA_plot$res_condition %in% 
                                   c("Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells_R"
                                   ),]
  if (length(unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Kathryn_TICA_plot$clusters)])>0) {
    
    add_needed <- unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Kathryn_TICA_plot$clusters)]
    add_needed <- data.frame(rep(c("Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Kathryn_T_NK_cells1.predicted.to_TICA_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Kathryn_TICA_plot)
    Kathryn_TICA_plot <- rbind(Kathryn_TICA_plot,add_needed)
  }
  Baolin_TICA_plot <- TICA_plot[TICA_plot$res_condition %in% 
                                  c("Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells_R"
                                  ),]
  if (length(unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Baolin_TICA_plot$clusters)])>0) {
    
    add_needed <- unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Baolin_TICA_plot$clusters)]
    add_needed <- data.frame(rep(c("Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Baolin_T_NK_cells1.predicted.to_TICA_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Baolin_TICA_plot)
    Baolin_TICA_plot <- rbind(Baolin_TICA_plot,add_needed)
  }
  
  Moshe_TICA_plot <- TICA_plot[TICA_plot$res_condition %in% 
                                 c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells_R"
                                 ),]
  if (length(unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Moshe_TICA_plot$clusters)])>0) {
    
    add_needed <- unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Moshe_TICA_plot$clusters)]
    add_needed <- data.frame(rep(c("Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells_NR","Moshe_T_NK_cells1.predicted.to_TICA_T_NK_cells_R"),length(add_needed)),
                                  unlist(lapply(add_needed,function(x) rep(x,2))),
                                  rep(NA,2*length(add_needed)),
                                  rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Moshe_TICA_plot)
    Moshe_TICA_plot <- rbind(Moshe_TICA_plot,add_needed)
  }
  Caushi_TICA_plot <- TICA_plot[TICA_plot$res_condition %in% 
                                  c("Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells_non-MPR","Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells_MPR"
                                  ),]
  if (length(unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Caushi_TICA_plot$clusters)])>0) {
    
    add_needed <- unique(TICA_plot$clusters)[!(unique(TICA_plot$clusters) %in% Caushi_TICA_plot$clusters)]
    add_needed <- data.frame(rep(c("Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells_non-MPR","Caushi_T_NK_cells1.predicted.to_TICA_T_NK_cells_MPR"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Caushi_TICA_plot)
    Caushi_TICA_plot <- rbind(Caushi_TICA_plot,add_needed)
  }
  pdf(paste0(workdir,"4_immune_response_plots_Kathryn_TICA_dotplot.pdf"),width = 4.88, height = 5.49)
  ggplot(data=Kathryn_TICA_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Baolin_TICA_dotplot.pdf"),width = 4.88, height = 5.49)
  ggplot(data=Baolin_TICA_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Moshe_TICA_dotplot.pdf"),width = 4.88, height = 5.49)
  ggplot(data=Moshe_TICA_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Caushi_TICA_dotplot.pdf"),width = 4.88, height = 5.49)
  ggplot(data=Caushi_TICA_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()

#### plot panT
  panT_plot <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% 
                                                    c("Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R",
                                                      "Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R",
                                                      "Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R",
                                                      "Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_non-MPR","Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_MPR"
                                                    ),];dim(panT_plot)
  Kathryn_panT_plot <- panT_plot[panT_plot$res_condition %in% 
                                   c("Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R"
                                   ),]
  if (length(unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Kathryn_panT_plot$clusters)])>0) {
    
    add_needed <- unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Kathryn_panT_plot$clusters)]
    add_needed <- data.frame(rep(c("Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Kathryn_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Kathryn_panT_plot)
    Kathryn_panT_plot <- rbind(Kathryn_panT_plot,add_needed)
  }
  Baolin_panT_plot <- panT_plot[panT_plot$res_condition %in% 
                                  c("Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R"
                                  ),]
  if (length(unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Baolin_panT_plot$clusters)])>0) {
    
    add_needed <- unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Baolin_panT_plot$clusters)]
    add_needed <- data.frame(rep(c("Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Baolin_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Baolin_panT_plot)
    Baolin_panT_plot <- rbind(Baolin_panT_plot,add_needed)
  }
  
  Moshe_panT_plot <- panT_plot[panT_plot$res_condition %in% 
                                 c("Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R"
                                 ),]
  if (length(unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Moshe_panT_plot$clusters)])>0) {
    
    add_needed <- unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Moshe_panT_plot$clusters)]
    add_needed <- data.frame(rep(c("Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_NR","Moshe_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Moshe_panT_plot)
    Moshe_panT_plot <- rbind(Moshe_panT_plot,add_needed)
  }
  Caushi_panT_plot <- panT_plot[panT_plot$res_condition %in% 
                                  c("Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_non-MPR","Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_MPR"
                                  ),]
  if (length(unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Caushi_panT_plot$clusters)])>0) {
    
    add_needed <- unique(panT_plot$clusters)[!(unique(panT_plot$clusters) %in% Caushi_panT_plot$clusters)]
    add_needed <- data.frame(rep(c("Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_non-MPR","Caushi_T_NK_cells2_predicted.to_pancancer_T_cells_CD8_CD4_MPR"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Caushi_panT_plot)
    Caushi_panT_plot <- rbind(Caushi_panT_plot,add_needed)
  }
  pdf(paste0(workdir,"4_immune_response_plots_Kathryn_panT_dotplot.pdf"),width = 10.86, height = 4.71)
  ggplot(data=Kathryn_panT_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Baolin_panT_dotplot.pdf"),width = 10.86, height = 4.71)
  ggplot(data=Baolin_panT_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Moshe_panT_dotplot.pdf"),width = 10.86, height = 4.71)
  ggplot(data=Moshe_panT_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Caushi_panT_dotplot.pdf"),width = 10.86, height = 4.71)
  ggplot(data=Caushi_panT_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  
##### panblueprint plot
  panblueprint_plot <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% 
                                                            c("Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R",
                                                              "Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R",
                                                              "Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R",
                                                              "Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_non-MPR","Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_MPR"
                                                            ),];dim(panblueprint_plot)
  Kathryn_panblueprint_plot <- panblueprint_plot[panblueprint_plot$res_condition %in% 
                                                   c("Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R"
                                                   ),]
  if (length(unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Kathryn_panblueprint_plot$clusters)])>0) {
    
    add_needed <- unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Kathryn_panblueprint_plot$clusters)]
    add_needed <- data.frame(rep(c("Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Kathryn_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Kathryn_panblueprint_plot)
    Kathryn_panblueprint_plot <- rbind(Kathryn_panblueprint_plot,add_needed)
  }
  Baolin_panblueprint_plot <- panblueprint_plot[panblueprint_plot$res_condition %in% 
                                                  c("Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R"
                                                  ),]
  if (length(unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Baolin_panblueprint_plot$clusters)])>0) {
    
    add_needed <- unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Baolin_panblueprint_plot$clusters)]
    add_needed <- data.frame(rep(c("Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Baolin_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Baolin_panblueprint_plot)
    Baolin_panblueprint_plot <- rbind(Baolin_panblueprint_plot,add_needed)
  }
  
  Moshe_panblueprint_plot <- panblueprint_plot[panblueprint_plot$res_condition %in% 
                                                 c("Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R"
                                                 ),]
  if (length(unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Moshe_panblueprint_plot$clusters)])>0) {
    
    add_needed <- unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Moshe_panblueprint_plot$clusters)]
    add_needed <- data.frame(rep(c("Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_NR","Moshe_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Moshe_panblueprint_plot)
    Moshe_panblueprint_plot <- rbind(Moshe_panblueprint_plot,add_needed)
  }
  Caushi_panblueprint_plot <- panblueprint_plot[panblueprint_plot$res_condition %in% 
                                                  c("Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_non-MPR","Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_MPR"
                                                  ),]
  if (length(unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Caushi_panblueprint_plot$clusters)])>0) {
    
    add_needed <- unique(panblueprint_plot$clusters)[!(unique(panblueprint_plot$clusters) %in% Caushi_panblueprint_plot$clusters)]
    add_needed <- data.frame(rep(c("Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_non-MPR","Caushi_T_NK_cells3_predicted.to_pancancer_blueprint_T_NK_cells_MPR"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Caushi_panblueprint_plot)
    Caushi_panblueprint_plot <- rbind(Caushi_panblueprint_plot,add_needed)
  }
  pdf(paste0(workdir,"4_immune_response_plots_Kathryn_panblueprint_dotplot.pdf"),width = 4.05, height = 4.71)
  ggplot(data=Kathryn_panblueprint_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Baolin_panblueprint_dotplot.pdf"),width = 4.05, height = 4.71)
  ggplot(data=Baolin_panblueprint_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Moshe_panblueprint_dotplot.pdf"),width = 4.05, height = 4.71)
  ggplot(data=Moshe_panblueprint_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Caushi_panblueprint_dotplot.pdf"),width = 4.05, height = 4.71)
  ggplot(data=Caushi_panblueprint_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  
##### ABTC plots
  ABTC_plot <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% 
                                                    c("Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R",
                                                      "Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R",
                                                      "Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R",
                                                      "Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells_non-MPR","Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells_MPR"
                                                    ),];dim(ABTC_plot)
  Kathryn_ABTC_plot <- ABTC_plot[ABTC_plot$res_condition %in% 
                                   c("Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R"
                                   ),]
  if (length(unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Kathryn_ABTC_plot$clusters)])>0) {
    
    add_needed <- unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Kathryn_ABTC_plot$clusters)]
    add_needed <- data.frame(rep(c("Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Kathryn_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Kathryn_ABTC_plot)
    Kathryn_ABTC_plot <- rbind(Kathryn_ABTC_plot,add_needed)
  }
  Baolin_ABTC_plot <- ABTC_plot[ABTC_plot$res_condition %in% 
                                  c("Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R"
                                  ),]
  if (length(unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Baolin_ABTC_plot$clusters)])>0) {
    
    add_needed <- unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Baolin_ABTC_plot$clusters)]
    add_needed <- data.frame(rep(c("Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Baolin_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Baolin_ABTC_plot)
    Baolin_ABTC_plot <- rbind(Baolin_ABTC_plot,add_needed)
  }
  
  Moshe_ABTC_plot <- ABTC_plot[ABTC_plot$res_condition %in% 
                                 c("Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R"
                                 ),]
  if (length(unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Moshe_ABTC_plot$clusters)])>0) {
    
    add_needed <- unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Moshe_ABTC_plot$clusters)]
    add_needed <- data.frame(rep(c("Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells_NR","Moshe_T_NK_cells4_predicted.to_ABTC_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Moshe_ABTC_plot)
    Moshe_ABTC_plot <- rbind(Moshe_ABTC_plot,add_needed)
  }
  Caushi_ABTC_plot <- ABTC_plot[ABTC_plot$res_condition %in% 
                                  c("Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells_non-MPR","Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells_MPR"
                                  ),]
  if (length(unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Caushi_ABTC_plot$clusters)])>0) {
    
    add_needed <- unique(ABTC_plot$clusters)[!(unique(ABTC_plot$clusters) %in% Caushi_ABTC_plot$clusters)]
    add_needed <- data.frame(rep(c("Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells_non-MPR","Caushi_T_NK_cells4_predicted.to_ABTC_T_NK_cells_MPR"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Caushi_ABTC_plot)
    Caushi_ABTC_plot <- rbind(Caushi_ABTC_plot,add_needed)
  }
  pdf(paste0(workdir,"4_immune_response_plots_Kathryn_ABTC_dotplot.pdf"),width = 4.28, height = 4.71)
  ggplot(data=Kathryn_ABTC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Baolin_ABTC_dotplot.pdf"),width = 4.28, height = 4.71)
  ggplot(data=Baolin_ABTC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Moshe_ABTC_dotplot.pdf"),width = 4.28, height = 4.71)
  ggplot(data=Moshe_ABTC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Caushi_ABTC_dotplot.pdf"),width = 4.28, height = 4.71)
  ggplot(data=Caushi_ABTC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  
##### HCC plots
  HCC_plot <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% 
                                                   c("Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells_R",
                                                     "Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells_R",
                                                     "Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells_R",
                                                     "Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells_non-MPR","Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells_MPR"
                                                   ),];dim(HCC_plot)
  Kathryn_HCC_plot <- HCC_plot[HCC_plot$res_condition %in% 
                                 c("Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells_R"
                                 ),]
  if (length(unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Kathryn_HCC_plot$clusters)])>0) {
    
    add_needed <- unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Kathryn_HCC_plot$clusters)]
    add_needed <- data.frame(rep(c("Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Kathryn_T_NK_cells5_predicted.to_HCC_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Kathryn_HCC_plot)
    Kathryn_HCC_plot <- rbind(Kathryn_HCC_plot,add_needed)
  }
  Baolin_HCC_plot <- HCC_plot[HCC_plot$res_condition %in% 
                                c("Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells_R"
                                ),]
  if (length(unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Baolin_HCC_plot$clusters)])>0) {
    
    add_needed <- unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Baolin_HCC_plot$clusters)]
    add_needed <- data.frame(rep(c("Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Baolin_T_NK_cells5_predicted.to_HCC_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Baolin_HCC_plot)
    Baolin_HCC_plot <- rbind(Baolin_HCC_plot,add_needed)
  }
  
  Moshe_HCC_plot <- HCC_plot[HCC_plot$res_condition %in% 
                               c("Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells_R"
                               ),]
  if (length(unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Moshe_HCC_plot$clusters)])>0) {
    
    add_needed <- unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Moshe_HCC_plot$clusters)]
    add_needed <- data.frame(rep(c("Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells_NR","Moshe_T_NK_cells5_predicted.to_HCC_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Moshe_HCC_plot)
    Moshe_HCC_plot <- rbind(Moshe_HCC_plot,add_needed)
  }
  Caushi_HCC_plot <- HCC_plot[HCC_plot$res_condition %in% 
                                c("Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells_non-MPR","Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells_MPR"
                                ),]
  if (length(unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Caushi_HCC_plot$clusters)])>0) {
    
    add_needed <- unique(HCC_plot$clusters)[!(unique(HCC_plot$clusters) %in% Caushi_HCC_plot$clusters)]
    add_needed <- data.frame(rep(c("Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells_non-MPR","Caushi_T_NK_cells5_predicted.to_HCC_T_NK_cells_MPR"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Caushi_HCC_plot)
    Caushi_HCC_plot <- rbind(Caushi_HCC_plot,add_needed)
  }
  pdf(paste0(workdir,"4_immune_response_plots_Kathryn_HCC_dotplot.pdf"),width = 3.53, height = 5.23)
  ggplot(data=Kathryn_HCC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Baolin_HCC_dotplot.pdf"),width = 3.53, height = 5.23)
  ggplot(data=Baolin_HCC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Moshe_HCC_dotplot.pdf"),width = 3.53, height = 5.23)
  ggplot(data=Moshe_HCC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Caushi_HCC_dotplot.pdf"),width = 3.53, height = 5.23)
  ggplot(data=Caushi_HCC_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  
##### STAD plot
  STAD_plot <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% 
                                                    c("Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells_R",
                                                      "Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells_R",
                                                      "Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells_R",
                                                      "Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells_non-MPR","Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells_MPR"
                                                    ),];dim(STAD_plot)
  Kathryn_STAD_plot <- STAD_plot[STAD_plot$res_condition %in% 
                                   c("Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells_R"
                                   ),]
  if (length(unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Kathryn_STAD_plot$clusters)])>0) {
    
    add_needed <- unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Kathryn_STAD_plot$clusters)]
    add_needed <- data.frame(rep(c("Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Kathryn_T_NK_cells6_predicted.to_STAD_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Kathryn_STAD_plot)
    Kathryn_STAD_plot <- rbind(Kathryn_STAD_plot,add_needed)
  }
  Baolin_STAD_plot <- STAD_plot[STAD_plot$res_condition %in% 
                                  c("Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells_R"
                                  ),]
  if (length(unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Baolin_STAD_plot$clusters)])>0) {
    
    add_needed <- unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Baolin_STAD_plot$clusters)]
    add_needed <- data.frame(rep(c("Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Baolin_T_NK_cells6_predicted.to_STAD_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Baolin_STAD_plot)
    Baolin_STAD_plot <- rbind(Baolin_STAD_plot,add_needed)
  }
  
  Moshe_STAD_plot <- STAD_plot[STAD_plot$res_condition %in% 
                                 c("Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells_R"
                                 ),]
  if (length(unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Moshe_STAD_plot$clusters)])>0) {
    
    add_needed <- unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Moshe_STAD_plot$clusters)]
    add_needed <- data.frame(rep(c("Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells_NR","Moshe_T_NK_cells6_predicted.to_STAD_T_NK_cells_R"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Moshe_STAD_plot)
    Moshe_STAD_plot <- rbind(Moshe_STAD_plot,add_needed)
  }
  Caushi_STAD_plot <- STAD_plot[STAD_plot$res_condition %in% 
                                  c("Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells_non-MPR","Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells_MPR"
                                  ),]
  if (length(unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Caushi_STAD_plot$clusters)])>0) {
    
    add_needed <- unique(STAD_plot$clusters)[!(unique(STAD_plot$clusters) %in% Caushi_STAD_plot$clusters)]
    add_needed <- data.frame(rep(c("Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells_non-MPR","Caushi_T_NK_cells6_predicted.to_STAD_T_NK_cells_MPR"),length(add_needed)),
                             unlist(lapply(add_needed,function(x) rep(x,2))),
                             rep(NA,2*length(add_needed)),
                             rep(NA,2*length(add_needed)))
    colnames(add_needed) <- colnames(Caushi_STAD_plot)
    Caushi_STAD_plot <- rbind(Caushi_STAD_plot,add_needed)
  }
  pdf(paste0(workdir,"4_immune_response_plots_Kathryn_STAD_dotplot.pdf"),width = 3.32, height = 4.59)
  ggplot(data=Kathryn_STAD_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Baolin_STAD_dotplot.pdf"),width = 3.32, height = 4.59)
  ggplot(data=Baolin_STAD_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Moshe_STAD_dotplot.pdf"),width = 3.32, height = 4.59)
  ggplot(data=Moshe_STAD_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  pdf(paste0(workdir,"4_immune_response_plots_Caushi_STAD_dotplot.pdf"),width = 3.32, height = 4.59)
  ggplot(data=Caushi_STAD_plot,
         mapping=aes_string(x='res_condition',y='clusters',color='NK_T_activation.Score')) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    theme_classic()
  dev.off()
  
  
  all[all$Kathryn_NK_T_activation.Score_p < 0.05 & all$Kathryn_NK_T_activation.Score_log2FC > 0.05,c(1,2,9,22,35,48)]
  all[all$Baolin_NK_T_activation.Score_p < 0.05 & all$Baolin_NK_T_activation.Score_log2FC > 0.05,c(1,15)]
  all[all$Moshe_NK_T_activation.Score_p < 0.05 & all$Moshe_NK_T_activation.Score_log2FC > 0.05,c(1,28)]
  all[all$Caushi_NK_T_activation.Score_p < 0.05 & all$Caushi_NK_T_activation.Score_log2FC > 0.05,c(1,41)]
  
  
#### plot for different resolution
  library(Seurat)
  library(SeuratObject)
  library(lisi)
  library(mclust)
  library(cluster)
  library(ggplot2)
  library(plyr)
  library(ggpubr)
  workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
  Kathryn_Baolin_Moshe_Caushi_output<- read.table(file=paste0(workdir,"4_immune_response_Kathryn_Baolin_Moshe_Caushi_to_reference_plot.txt"),header=T,sep="\t")
  Kathryn_resolution <- c(0.1,0.15,0.16,0.2,0.217,0.22,0.25,0.28)
  multi_plot <- list()
  for (i in 1:length(Kathryn_resolution)) {
  Kathryn_plot <- c(paste0(paste0("Kathryn_T_NK_cells1.RNA_snn_res.",Kathryn_resolution[i]),"_R"),paste0(paste0("Kathryn_T_NK_cells1.RNA_snn_res.",Kathryn_resolution[i]),"_NR"))
Kathryn_plot_data <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% Kathryn_plot,]
  Kathryn_plot_data$clusters <- factor(Kathryn_plot_data$clusters,levels = 0:26)
  # pdf(paste0(workdir,"4_immune_response_plots_Kathryn_resolutions",Kathryn_resolution[i],"_dotplot.pdf"),width = 3, height = 5.74)
  multi_plot[[i]] <- ggplot(data=Kathryn_plot_data,
         mapping=aes_string(x='res_condition',y='clusters',color="NK_T_activation.Score")) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    labs (y="") +
    theme_classic() +
    theme(text = element_text(size=20, colour = "black"))
  # dev.off()
  }
  pdf(paste0(workdir,"4_immune_response_plots_Kathryn_resolutions_dotplot_combined.pdf"),width = 10, height = 5.74)
  ggarrange(multi_plot[[1]],multi_plot[[2]],multi_plot[[3]],multi_plot[[4]],multi_plot[[5]],multi_plot[[6]],multi_plot[[7]],ncol=7,nrow=1,legend = F)
  dev.off()
  Kathryn_plot_data_p <- read.table(file=paste0(workdir,"4_immune_response_Kathryn.txt"),header=T,sep="\t")
  Kathryn_plot_data_p <- Kathryn_plot_data_p[Kathryn_plot_data_p$X %in% paste0("Kathryn_T_NK_cells1.RNA_snn_res.",Kathryn_resolution),]
  Kathryn_plot_data_p[Kathryn_plot_data_p$NK_T_activation.Score_p < 0.05 & Kathryn_plot_data_p$NK_T_activation.Score_log2FC > 0.05,c(1:5,10:11)]
  
  
  Baolin_resolution <- c(0.1,0.15,0.2,0.21,0.22,0.3,0.4)
  multi_plot <- list()
  for (i in 1:length(Baolin_resolution)) {
  Baolin_plot <- c(paste0(paste0("Baolin_T_NK_cells1.RNA_snn_res.",Baolin_resolution[i]),"_R"),paste0(paste0("Baolin_T_NK_cells1.RNA_snn_res.",Baolin_resolution[i]),"_NR"))
  Baolin_plot_data <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% Baolin_plot,]
  Baolin_plot_data$clusters <- factor(Baolin_plot_data$clusters,levels = 0:28)
  # pdf(paste0(workdir,"4_immune_response_plots_Baolin_resolutions_dotplot.pdf"),width = 5.5, height = 5.74)
  multi_plot[[i]] <-ggplot(data=Baolin_plot_data,
         mapping=aes_string(x='res_condition',y='clusters',color="NK_T_activation.Score")) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    labs (y="") +
    theme_classic() +
    theme(text = element_text(size=20, colour = "black"))
  # dev.off()
  }
  pdf(paste0(workdir,"4_immune_response_plots_Baolin_resolutions_dotplot_combined.pdf"),width = 10, height = 5.74)
  ggarrange(multi_plot[[1]],multi_plot[[2]],multi_plot[[3]],multi_plot[[4]],multi_plot[[5]],multi_plot[[6]],multi_plot[[7]],ncol=7,nrow=1,legend = F)
  dev.off()
  
  Baolin_plot_data_p <- read.table(file=paste0(workdir,"4_immune_response_Baolin.txt"),header=T,sep="\t")
  Baolin_plot_data_p <- Baolin_plot_data_p[Baolin_plot_data_p$X %in% paste0("Baolin_T_NK_cells1.RNA_snn_res.",Baolin_resolution),]
  Baolin_plot_data_p[Baolin_plot_data_p$NK_T_activation.Score_p < 0.05 & Baolin_plot_data_p$NK_T_activation.Score_log2FC > 0.05,c(1:5,10:11)]
  
  
  
  Moshe_resolution <- c(0.35,0.4,0.52,0.6,0.8,0.82,0.8245)
  multi_plot <- list()
  for (i in 1:length(Moshe_resolution)) {
  Moshe_plot <- c(paste0(paste0("Moshe_T_NK_cells1.RNA_snn_res.",Moshe_resolution[i]),"_R"),paste0(paste0("Moshe_T_NK_cells1.RNA_snn_res.",Moshe_resolution[i]),"_NR"))
  Moshe_plot_data <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% Moshe_plot,]
  Moshe_plot_data$clusters <- factor(Moshe_plot_data$clusters,levels = 0:28)
  # pdf(paste0(workdir,"4_immune_response_plots_Moshe_resolutions_dotplot.pdf"),width = 5.5, height = 5.74)
  multi_plot[[i]] <- ggplot(data=Moshe_plot_data,
         mapping=aes_string(x='res_condition',y='clusters',color="NK_T_activation.Score")) + 
    geom_point(mapping=aes_string(size='Freq')) +
    labs (y="") +
    scale_color_gradient(low="blue",high="red") +
    theme_classic() + 
    theme(text = element_text(size=20, colour = "black"))
  # dev.off()
  }
  pdf(paste0(workdir,"4_immune_response_plots_Moshe_resolutions_dotplot_combined.pdf"),width = 10, height = 5.74)
  ggarrange(multi_plot[[1]],multi_plot[[2]],multi_plot[[3]],multi_plot[[4]],multi_plot[[5]],multi_plot[[6]],multi_plot[[7]],ncol=7,nrow=1,legend = F)
  dev.off()
  
  Moshe_plot_data_p <- read.table(file=paste0(workdir,"4_immune_response_Moshe.txt"),header=T,sep="\t")
  Moshe_plot_data_p <- Moshe_plot_data_p[Moshe_plot_data_p$X %in% paste0("Moshe_T_NK_cells1.RNA_snn_res.",Moshe_resolution),]
  Moshe_plot_data_p[Moshe_plot_data_p$NK_T_activation.Score_p < 0.05 & Moshe_plot_data_p$NK_T_activation.Score_log2FC > 0.05,c(1:5,10:11)]
  
  
  
  Caushi_resolution <- c(0.2,0.208,0.21,0.3,0.31,0.35,0.4)
  multi_plot <- list()
  for (i in 1:length(Caushi_resolution)) {
  Caushi_plot <- c(paste0(paste0("Caushi_T_NK_cells1.RNA_snn_res.",Caushi_resolution[i]),"_MPR"),paste0(paste0("Caushi_T_NK_cells1.RNA_snn_res.",Caushi_resolution[i]),"_non-MPR"))
  Caushi_plot_data <- Kathryn_Baolin_Moshe_Caushi_output[Kathryn_Baolin_Moshe_Caushi_output$res_condition %in% Caushi_plot,]
  Caushi_plot_data$clusters <- factor(Caushi_plot_data$clusters,levels = 0:27)
  Caushi_plot_data$res_condition <- factor(Caushi_plot_data$res_condition,levels = c("Caushi_T_NK_cells1.RNA_snn_res.0.2_non-MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.2_MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.208_non-MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.208_MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.21_non-MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.21_MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.3_non-MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.3_MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.31_non-MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.31_MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.35_non-MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.35_MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.4_non-MPR",
                                                                                     "Caushi_T_NK_cells1.RNA_snn_res.0.4_MPR"#,
                                                                                     ))
  # pdf(paste0(workdir,"4_immune_response_plots_Caushi_resolutions_dotplot.pdf"),width = 5.5, height = 5.74)
  multi_plot[[i]] <- ggplot(data=Caushi_plot_data,
         mapping=aes_string(x='res_condition',y='clusters',color="NK_T_activation.Score")) + 
    geom_point(mapping=aes_string(size='Freq')) +
    scale_color_gradient(low="blue",high="red") +
    labs (y="") +
    theme_classic() +
    theme(text = element_text(size=20, colour = "black"))
  # dev.off()
  }
  pdf(paste0(workdir,"4_immune_response_plots_Caushi_resolutions_dotplot_combined.pdf"),width = 10, height = 5.74)
  ggarrange(multi_plot[[1]],multi_plot[[2]],multi_plot[[3]],multi_plot[[4]],multi_plot[[5]],multi_plot[[6]],multi_plot[[7]],ncol=7,nrow=1,legend = F)
  dev.off()
  
  Caushi_plot_data_p <- read.table(file=paste0(workdir,"4_immune_response_Caushi.txt"),header=T,sep="\t")
  Caushi_plot_data_p <- Caushi_plot_data_p[Caushi_plot_data_p$X %in% paste0("Caushi_T_NK_cells1.RNA_snn_res.",Caushi_resolution),]
  Caushi_plot_data_p[Caushi_plot_data_p$NK_T_activation.Score_p < 0.05 & Caushi_plot_data_p$NK_T_activation.Score_log2FC > 0.05,c(1:5,10:11)]
  
