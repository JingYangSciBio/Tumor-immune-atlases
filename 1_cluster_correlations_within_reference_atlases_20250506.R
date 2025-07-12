#########################################
#### Jing Yang
#### 2024.02.14 Figure 1 part 1. #reorder plots on output_cor, output_ove and output_gsea_es 
#########################################
library(purrr)
library(tidyverse)
library(scales)
library(Seurat)
library(infotheo)
library(pheatmap)
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"

output_ove <- readRDS(paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_ove_plot_20230606.rds")) #20240214

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

names(new_cluster_names) <- c("TICA_T_NK_cells__Recently activated CD4 T cells", "TICA_T_NK_cells__Naive-memory CD4 T cells", "TICA_T_NK_cells__Transitional memory CD4 T cells","TICA_T_NK_cells__Cytotoxic CD8 T cells",
                       "TICA_T_NK_cells__Effector memory CD8 T cells", "TICA_T_NK_cells__Th17 cells", "TICA_T_NK_cells__NK", "TICA_T_NK_cells__Terminally exhausted CD8 T cells",
                       "TICA_T_NK_cells__Naive T cells", "TICA_T_NK_cells__Regulatory T cells", "TICA_T_NK_cells__T helper cells", "TICA_T_NK_cells__Proliferative T cells",
                       "TICA_T_NK_cells__Pre-exhausted CD8 T cells",
                       
                       
                       "pancancer_T_cells_CD8_CD4__CD8.c07.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD8.c09.Tk.KIR2DL4", "pancancer_T_cells_CD8_CD4__CD8.c05.Tem.CXCR5" ,
                       "pancancer_T_cells_CD8_CD4__CD8.c15.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD8.c01.Tn.MAL", "pancancer_T_cells_CD8_CD4__CD8.c17.Tm.NME1", "pancancer_T_cells_CD8_CD4__CD8.c08.Tk.TYROBP" ,
                       "pancancer_T_cells_CD8_CD4__CD8.c12.Tex.CXCL13", "pancancer_T_cells_CD8_CD4__CD8.c02.Tm.IL7R", "pancancer_T_cells_CD8_CD4__CD8.c10.Trm.ZNF683", "pancancer_T_cells_CD8_CD4__CD8.c03.Tm.RPS12" ,
                       "pancancer_T_cells_CD8_CD4__CD8.c14.Tex.TCF7", "pancancer_T_cells_CD8_CD4__CD8.c06.Tem.GZMK", "pancancer_T_cells_CD8_CD4__CD8.c11.Tex.PDCD1", "pancancer_T_cells_CD8_CD4__CD8.c16.MAIT.SLC4A10" ,
                       "pancancer_T_cells_CD8_CD4__CD8.c04.Tm.CD52", "pancancer_T_cells_CD8_CD4__CD8.c13.Tex.myl12a", "pancancer_T_cells_CD8_CD4__CD4.c20.Treg.TNFRSF9", "pancancer_T_cells_CD8_CD4__CD4.c18.Treg.RTKN2" ,
                       "pancancer_T_cells_CD8_CD4__CD4.c06.Tm.ANXA1", "pancancer_T_cells_CD8_CD4__CD4.c22.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD4.c01.Tn.TCF7", "pancancer_T_cells_CD8_CD4__CD4.c12.Tem.GZMK" ,
                       "pancancer_T_cells_CD8_CD4__CD4.c11.Tm.GZMA", "pancancer_T_cells_CD8_CD4__CD4.c13.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD4.c16.Tfh.CXCR5", "pancancer_T_cells_CD8_CD4__CD4.c07.Tm.ANXA2" ,
                       "pancancer_T_cells_CD8_CD4__CD4.c17.TfhTh1.CXCL13", "pancancer_T_cells_CD8_CD4__CD4.c19.Treg.S1PR1", "pancancer_T_cells_CD8_CD4__CD4.c10.Tm.CAPG", "pancancer_T_cells_CD8_CD4__CD4.c24.Mix.NME2" ,
                       "pancancer_T_cells_CD8_CD4__CD4.c21.Treg.OAS1", "pancancer_T_cells_CD8_CD4__CD4.c14.Th17.SLC4A10", "pancancer_T_cells_CD8_CD4__CD4.c09.Tm.CCL5", "pancancer_T_cells_CD8_CD4__CD4.c02.Tn.PASK" ,
                       "pancancer_T_cells_CD8_CD4__CD4.c03.Tn.ADSL", "pancancer_T_cells_CD8_CD4__CD4.c05.Tm.TNF", "pancancer_T_cells_CD8_CD4__CD4.c08.Tm.CREM", "pancancer_T_cells_CD8_CD4__CD4.c04.Tn.il7r" ,
                       "pancancer_T_cells_CD8_CD4__CD4.c15.Th17.IL23R", "pancancer_T_cells_CD8_CD4__CD4.c23.Mix.NME1",
                       
                       
                       "pancancer_blueprint_T_NK_cells__CD8 memory T", "pancancer_blueprint_T_NK_cells__CD8 exhausted cytotoxic T",
                       "pancancer_blueprint_T_NK_cells__CD8 pre-effector T",   "pancancer_blueprint_T_NK_cells__CD4 memory/effector T",
                       "pancancer_blueprint_T_NK_cells__NK_XCL1",   "pancancer_blueprint_T_NK_cells__CD4 regulatory T",
                       "pancancer_blueprint_T_NK_cells__CD4 naive T",   "pancancer_blueprint_T_NK_cells__CD8 effector T",
                       "pancancer_blueprint_T_NK_cells__CD4 exhausted effector T", "pancancer_blueprint_T_NK_cells__Cytotoxic NK",
                       
                       "ABTC_T_NK_cells__Naive CD4 T cell",     "ABTC_T_NK_cells__Exhausted CD8 T cell",
                       "ABTC_T_NK_cells__Memory CD8 T cell",    "ABTC_T_NK_cells__Treg" ,               
                       "ABTC_T_NK_cells__DNT",                  "ABTC_T_NK_cells__Effector CD8 T cell" ,
                       "ABTC_T_NK_cells__NCAM1nFCGR3Ap NK",     "ABTC_T_NK_cells__Naive CD8 T cell" ,   
                       "ABTC_T_NK_cells__NKT and NK cell",      "ABTC_T_NK_cells__NCAM1pFCGR3Ap NK"   ,"ABTC_T_NK_cells__NCAM1pFCGR3An NK",
                       
                       "HCC_T_NK_cells__13_CD8_effector memory T", "HCC_T_NK_cells__2_NK", "HCC_T_NK_cells__2_CD8_MAIT",  "HCC_T_NK_cells__20_CD4_naive T",
                       "HCC_T_NK_cells__9_CD4_regulatory T",  "HCC_T_NK_cells__8_NK",  "HCC_T_NK_cells__35_CD4_central memory T", "HCC_T_NK_cells__18_NK",
                       "HCC_T_NK_cells__11_CD8_inter-states T", "HCC_T_NK_cells__14_CD8_TR_memory T", "HCC_T_NK_cells__32_CD4_other T", "HCC_T_NK_cells__1_CD8_cytotoxic T" ,
                       "HCC_T_NK_cells__28_NK", "HCC_T_NK_cells__51_CD8_other T",
                       
                       "STAD_T_NK_cells__LNK1", "STAD_T_NK_cells__LNK3", "STAD_T_NK_cells__LNK4", "STAD_T_NK_cells__CD8 naive T" ,
                       "STAD_T_NK_cells__CD4 naive T",  "STAD_T_NK_cells__LNK2", "STAD_T_NK_cells__CD4 helper T",  "STAD_T_NK_cells__Regulatory T" ,
                       "STAD_T_NK_cells__CD8 effector T",  "STAD_T_NK_cells__Proliferative T")

# rownames(output_ove) <- as.character(new_cluster_names[rownames(output_ove)])
# colnames(output_ove) <- as.character(new_cluster_names[colnames(output_ove)])
# rownames(output_ove)[which(rownames(output_ove)=="HCC_T_NK_cells: 14_TR_memory T")] <- "HCC_T_NK_cells: 14_CD8_TR_memory T"
# colnames(output_ove)[which(colnames(output_ove)=="HCC_T_NK_cells: 14_TR_memory T")] <- "HCC_T_NK_cells: 14_CD8_TR_memory T"
# rownames(output_ove)[which(rownames(output_ove)=="HCC_T_NK_cells: 35_central memory T")] <- "HCC_T_NK_cells: 35_CD4_central memory T"
# colnames(output_ove)[which(colnames(output_ove)=="HCC_T_NK_cells: 35_central memory T")] <- "HCC_T_NK_cells: 35_CD4_central memory T"
# output_ove <- output_ove[c(1:32,34:58),c(1:32,34:58)]
# ##### CD4
# # output_ove <- output_ove[rownames(output_ove) %in% c("pancancer_blueprint_T_NK_cells__CD4 exhausted effector T", "pancancer_T_cells_CD8_CD4__CD4.c16.Tfh.CXCR5", "pancancer_T_cells_CD8_CD4__CD4.c17.TfhTh1.CXCL13", "pancancer_T_cells_CD8_CD4__CD4.c02.Tn.PASK", "pancancer_T_cells_CD8_CD4__CD4.c03.Tn.ADSL", "pancancer_T_cells_CD8_CD4__CD4.c04.Tn.il7r", "TICA_T_NK_cells__Naive-memory CD4 T cells", "STAD_T_NK_cells__CD4 naive T", "pancancer_blueprint_T_NK_cells__CD4 naive T", "HCC_T_NK_cells__20_CD4_naive T", "pancancer_T_cells_CD8_CD4__CD4.c01.Tn.TCF7", "ABTC_T_NK_cells__Naive CD4 T cell", "TICA_T_NK_cells__Naive T cells", "pancancer_T_cells_CD8_CD4__CD4.c21.Treg.OAS1", "pancancer_T_cells_CD8_CD4__CD4.c22.ISG.IFIT1", "STAD_T_NK_cells__Proliferative T", "TICA_T_NK_cells__Proliferative T cells", "HCC_T_NK_cells__32_CD4_other T", "TICA_T_NK_cells__T helper cells", "STAD_T_NK_cells__Regulatory T", "HCC_T_NK_cells__9_CD4_regulatory T", "pancancer_blueprint_T_NK_cells__CD4 regulatory T", "pancancer_T_cells_CD8_CD4__CD4.c20.Treg.TNFRSF9", "TICA_T_NK_cells__Regulatory T cells", "ABTC_T_NK_cells__Treg", "pancancer_T_cells_CD8_CD4__CD4.c12.Tem.GZMK", "pancancer_T_cells_CD8_CD4__CD4.c13.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD4.c10.Tm.CAPG", "pancancer_T_cells_CD8_CD4__CD4.c08.Tm.CREM", "STAD_T_NK_cells__CD4 helper T", "pancancer_blueprint_T_NK_cells__CD4 memory/effector T", "TICA_T_NK_cells__Transitional memory CD4 T cells", "pancancer_T_cells_CD8_CD4__CD4.c18.Treg.RTKN2", "pancancer_T_cells_CD8_CD4__CD4.c19.Treg.S1PR1", "TICA_T_NK_cells__Recently activated CD4 T cells", "HCC_T_NK_cells__35_CD4_central memory T"),
# #                          colnames(output_ove) %in% c("pancancer_blueprint_T_NK_cells__CD4 exhausted effector T", "pancancer_T_cells_CD8_CD4__CD4.c16.Tfh.CXCR5", "pancancer_T_cells_CD8_CD4__CD4.c17.TfhTh1.CXCL13", "pancancer_T_cells_CD8_CD4__CD4.c02.Tn.PASK", "pancancer_T_cells_CD8_CD4__CD4.c03.Tn.ADSL", "pancancer_T_cells_CD8_CD4__CD4.c04.Tn.il7r", "TICA_T_NK_cells__Naive-memory CD4 T cells", "STAD_T_NK_cells__CD4 naive T", "pancancer_blueprint_T_NK_cells__CD4 naive T", "HCC_T_NK_cells__20_CD4_naive T", "pancancer_T_cells_CD8_CD4__CD4.c01.Tn.TCF7", "ABTC_T_NK_cells__Naive CD4 T cell", "TICA_T_NK_cells__Naive T cells", "pancancer_T_cells_CD8_CD4__CD4.c21.Treg.OAS1", "pancancer_T_cells_CD8_CD4__CD4.c22.ISG.IFIT1", "STAD_T_NK_cells__Proliferative T", "TICA_T_NK_cells__Proliferative T cells", "HCC_T_NK_cells__32_CD4_other T", "TICA_T_NK_cells__T helper cells", "STAD_T_NK_cells__Regulatory T", "HCC_T_NK_cells__9_CD4_regulatory T", "pancancer_blueprint_T_NK_cells__CD4 regulatory T", "pancancer_T_cells_CD8_CD4__CD4.c20.Treg.TNFRSF9", "TICA_T_NK_cells__Regulatory T cells", "ABTC_T_NK_cells__Treg", "pancancer_T_cells_CD8_CD4__CD4.c12.Tem.GZMK", "pancancer_T_cells_CD8_CD4__CD4.c13.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD4.c10.Tm.CAPG", "pancancer_T_cells_CD8_CD4__CD4.c08.Tm.CREM", "STAD_T_NK_cells__CD4 helper T", "pancancer_blueprint_T_NK_cells__CD4 memory/effector T", "TICA_T_NK_cells__Transitional memory CD4 T cells", "pancancer_T_cells_CD8_CD4__CD4.c18.Treg.RTKN2", "pancancer_T_cells_CD8_CD4__CD4.c19.Treg.S1PR1", "TICA_T_NK_cells__Recently activated CD4 T cells", "HCC_T_NK_cells__35_CD4_central memory T")]
# output_ove <- output_ove[rownames(output_ove) %in% c("TICA_T_NK_cells__Recently activated CD4 T cells", "TICA_T_NK_cells__Naive-memory CD4 T cells", "TICA_T_NK_cells__Transitional memory CD4 T cells","TICA_T_NK_cells__Naive T cells", "TICA_T_NK_cells__Regulatory T cells", "TICA_T_NK_cells__T helper cells", "TICA_T_NK_cells__Proliferative T cells","TICA_T_NK_cells__Th17 cells",
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c06.Tm.ANXA1", "pancancer_T_cells_CD8_CD4__CD4.c22.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD4.c01.Tn.TCF7", "pancancer_T_cells_CD8_CD4__CD4.c12.Tem.GZMK" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c11.Tm.GZMA", "pancancer_T_cells_CD8_CD4__CD4.c13.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD4.c16.Tfh.CXCR5", "pancancer_T_cells_CD8_CD4__CD4.c07.Tm.ANXA2" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c17.TfhTh1.CXCL13", "pancancer_T_cells_CD8_CD4__CD4.c19.Treg.S1PR1", "pancancer_T_cells_CD8_CD4__CD4.c10.Tm.CAPG", "pancancer_T_cells_CD8_CD4__CD4.c24.Mix.NME2" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c21.Treg.OAS1", "pancancer_T_cells_CD8_CD4__CD4.c14.Th17.SLC4A10", "pancancer_T_cells_CD8_CD4__CD4.c09.Tm.CCL5", "pancancer_T_cells_CD8_CD4__CD4.c02.Tn.PASK" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c03.Tn.ADSL", "pancancer_T_cells_CD8_CD4__CD4.c05.Tm.TNF", "pancancer_T_cells_CD8_CD4__CD4.c08.Tm.CREM", "pancancer_T_cells_CD8_CD4__CD4.c04.Tn.il7r" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c15.Th17.IL23R", "pancancer_T_cells_CD8_CD4__CD4.c23.Mix.NME1","pancancer_T_cells_CD8_CD4__CD4.c20.Treg.TNFRSF9", "pancancer_T_cells_CD8_CD4__CD4.c18.Treg.RTKN2" ,
#                                                      "pancancer_blueprint_T_NK_cells__CD4 memory/effector T","pancancer_blueprint_T_NK_cells__CD4 naive T",  "pancancer_blueprint_T_NK_cells__CD4 regulatory T", "pancancer_blueprint_T_NK_cells__CD4 exhausted effector T",
#                                                      "ABTC_T_NK_cells__Naive CD4 T cell","ABTC_T_NK_cells__Treg" ,
#                                                      "HCC_T_NK_cells__20_CD4_naive T","HCC_T_NK_cells__9_CD4_regulatory T","HCC_T_NK_cells__35_CD4_central memory T","HCC_T_NK_cells__32_CD4_other T",
#                                                      "STAD_T_NK_cells__CD4 naive T", "STAD_T_NK_cells__CD4 helper T",  "STAD_T_NK_cells__Regulatory T" ,"STAD_T_NK_cells__Proliferative T"),
#                          colnames(output_ove) %in% c("TICA_T_NK_cells__Recently activated CD4 T cells", "TICA_T_NK_cells__Naive-memory CD4 T cells", "TICA_T_NK_cells__Transitional memory CD4 T cells","TICA_T_NK_cells__Naive T cells", "TICA_T_NK_cells__Regulatory T cells", "TICA_T_NK_cells__T helper cells", "TICA_T_NK_cells__Proliferative T cells","TICA_T_NK_cells__Th17 cells",
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c06.Tm.ANXA1", "pancancer_T_cells_CD8_CD4__CD4.c22.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD4.c01.Tn.TCF7", "pancancer_T_cells_CD8_CD4__CD4.c12.Tem.GZMK" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c11.Tm.GZMA", "pancancer_T_cells_CD8_CD4__CD4.c13.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD4.c16.Tfh.CXCR5", "pancancer_T_cells_CD8_CD4__CD4.c07.Tm.ANXA2" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c17.TfhTh1.CXCL13", "pancancer_T_cells_CD8_CD4__CD4.c19.Treg.S1PR1", "pancancer_T_cells_CD8_CD4__CD4.c10.Tm.CAPG", "pancancer_T_cells_CD8_CD4__CD4.c24.Mix.NME2" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c21.Treg.OAS1", "pancancer_T_cells_CD8_CD4__CD4.c14.Th17.SLC4A10", "pancancer_T_cells_CD8_CD4__CD4.c09.Tm.CCL5", "pancancer_T_cells_CD8_CD4__CD4.c02.Tn.PASK" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c03.Tn.ADSL", "pancancer_T_cells_CD8_CD4__CD4.c05.Tm.TNF", "pancancer_T_cells_CD8_CD4__CD4.c08.Tm.CREM", "pancancer_T_cells_CD8_CD4__CD4.c04.Tn.il7r" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD4.c15.Th17.IL23R", "pancancer_T_cells_CD8_CD4__CD4.c23.Mix.NME1","pancancer_T_cells_CD8_CD4__CD4.c20.Treg.TNFRSF9", "pancancer_T_cells_CD8_CD4__CD4.c18.Treg.RTKN2" ,
#                                                      "pancancer_blueprint_T_NK_cells__CD4 memory/effector T","pancancer_blueprint_T_NK_cells__CD4 naive T",  "pancancer_blueprint_T_NK_cells__CD4 regulatory T", "pancancer_blueprint_T_NK_cells__CD4 exhausted effector T",
#                                                      "ABTC_T_NK_cells__Naive CD4 T cell","ABTC_T_NK_cells__Treg" ,
#                                                      "HCC_T_NK_cells__20_CD4_naive T","HCC_T_NK_cells__9_CD4_regulatory T","HCC_T_NK_cells__35_CD4_central memory T","HCC_T_NK_cells__32_CD4_other T",
#                                                      "STAD_T_NK_cells__CD4 naive T", "STAD_T_NK_cells__CD4 helper T",  "STAD_T_NK_cells__Regulatory T" ,"STAD_T_NK_cells__Proliferative T")]
# dim(output_ove)
# ##### CD8
# # output_ove <- output_ove[rownames(output_ove) %in% c("pancancer_T_cells_CD8_CD4__CD8.c14.Tex.TCF7", "pancancer_T_cells_CD8_CD4__CD8.c12.Tex.CXCL13", "pancancer_T_cells_CD8_CD4__CD8.c13.Tex.myl12a", "pancancer_T_cells_CD8_CD4__CD8.c11.Tex.PDCD1", "ABTC_T_NK_cells__Exhausted CD8 T cell", "TICA_T_NK_cells__Terminally exhausted CD8 T cells", "pancancer_blueprint_T_NK_cells__CD8 exhausted cytotoxic T", "pancancer_T_cells_CD8_CD4__CD8.c01.Tn.MAL", "ABTC_T_NK_cells__Naive CD8 T cell", "TICA_T_NK_cells__Pre-exhausted CD8 T cells", "pancancer_T_cells_CD8_CD4__CD8.c15.ISG.IFIT1", "HCC_T_NK_cells__51_CD8_other T", "HCC_T_NK_cells__11_CD8_inter-states T", "HCC_T_NK_cells__1_CD8_cytotoxic T", "pancancer_T_cells_CD8_CD4__CD8.c05.Tem.CXCR5", "pancancer_T_cells_CD8_CD4__CD8.c06.Tem.GZMK", "pancancer_blueprint_T_NK_cells__CD8 pre-effector T", "ABTC_T_NK_cells__Memory CD8 T cell", "pancancer_blueprint_T_NK_cells__CD8 effector T", "HCC_T_NK_cells__13_CD8_effector memory T", "pancancer_T_cells_CD8_CD4__CD8.c07.Temra.CX3CR1", "ABTC_T_NK_cells__Effector CD8 T cell", "pancancer_T_cells_CD8_CD4__CD8.c10.Trm.ZNF683", "pancancer_blueprint_T_NK_cells__CD8 memory T", "TICA_T_NK_cells__Cytotoxic CD8 T cells", "pancancer_T_cells_CD8_CD4__CD8.c02.Tm.IL7R", "STAD_T_NK_cells__CD8 naive T", "HCC_T_NK_cells__2_CD8_MAIT", "STAD_T_NK_cells__CD8 effector T", "TICA_T_NK_cells__Effector memory CD8 T cells", "HCC_T_NK_cells__14_CD8_TR_memory T"),
# #                          colnames(output_ove) %in% c("pancancer_T_cells_CD8_CD4__CD8.c14.Tex.TCF7", "pancancer_T_cells_CD8_CD4__CD8.c12.Tex.CXCL13", "pancancer_T_cells_CD8_CD4__CD8.c13.Tex.myl12a", "pancancer_T_cells_CD8_CD4__CD8.c11.Tex.PDCD1", "ABTC_T_NK_cells__Exhausted CD8 T cell", "TICA_T_NK_cells__Terminally exhausted CD8 T cells", "pancancer_blueprint_T_NK_cells__CD8 exhausted cytotoxic T", "pancancer_T_cells_CD8_CD4__CD8.c01.Tn.MAL", "ABTC_T_NK_cells__Naive CD8 T cell", "TICA_T_NK_cells__Pre-exhausted CD8 T cells", "pancancer_T_cells_CD8_CD4__CD8.c15.ISG.IFIT1", "HCC_T_NK_cells__51_CD8_other T", "HCC_T_NK_cells__11_CD8_inter-states T", "HCC_T_NK_cells__1_CD8_cytotoxic T", "pancancer_T_cells_CD8_CD4__CD8.c05.Tem.CXCR5", "pancancer_T_cells_CD8_CD4__CD8.c06.Tem.GZMK", "pancancer_blueprint_T_NK_cells__CD8 pre-effector T", "ABTC_T_NK_cells__Memory CD8 T cell", "pancancer_blueprint_T_NK_cells__CD8 effector T", "HCC_T_NK_cells__13_CD8_effector memory T", "pancancer_T_cells_CD8_CD4__CD8.c07.Temra.CX3CR1", "ABTC_T_NK_cells__Effector CD8 T cell", "pancancer_T_cells_CD8_CD4__CD8.c10.Trm.ZNF683", "pancancer_blueprint_T_NK_cells__CD8 memory T", "TICA_T_NK_cells__Cytotoxic CD8 T cells", "pancancer_T_cells_CD8_CD4__CD8.c02.Tm.IL7R", "STAD_T_NK_cells__CD8 naive T", "HCC_T_NK_cells__2_CD8_MAIT", "STAD_T_NK_cells__CD8 effector T", "TICA_T_NK_cells__Effector memory CD8 T cells", "HCC_T_NK_cells__14_CD8_TR_memory T")]
# output_ove <- output_ove[rownames(output_ove) %in% c("TICA_T_NK_cells__Cytotoxic CD8 T cells","TICA_T_NK_cells__Terminally exhausted CD8 T cells","TICA_T_NK_cells__Pre-exhausted CD8 T cells", "TICA_T_NK_cells__Effector memory CD8 T cells",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c07.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD8.c05.Tem.CXCR5" ,"pancancer_T_cells_CD8_CD4__CD8.c04.Tm.CD52", "pancancer_T_cells_CD8_CD4__CD8.c13.Tex.myl12a",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c15.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD8.c01.Tn.MAL", "pancancer_T_cells_CD8_CD4__CD8.c17.Tm.NME1",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c12.Tex.CXCL13", "pancancer_T_cells_CD8_CD4__CD8.c02.Tm.IL7R", "pancancer_T_cells_CD8_CD4__CD8.c10.Trm.ZNF683", "pancancer_T_cells_CD8_CD4__CD8.c03.Tm.RPS12" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c14.Tex.TCF7", "pancancer_T_cells_CD8_CD4__CD8.c06.Tem.GZMK", "pancancer_T_cells_CD8_CD4__CD8.c11.Tex.PDCD1", "pancancer_T_cells_CD8_CD4__CD8.c16.MAIT.SLC4A10" ,
#                                                      "pancancer_blueprint_T_NK_cells__CD8 memory T", "pancancer_blueprint_T_NK_cells__CD8 exhausted cytotoxic T","pancancer_blueprint_T_NK_cells__CD8 pre-effector T",  "pancancer_blueprint_T_NK_cells__CD8 effector T",
#                                                      "ABTC_T_NK_cells__Exhausted CD8 T cell","ABTC_T_NK_cells__Memory CD8 T cell", "ABTC_T_NK_cells__Effector CD8 T cell" , "ABTC_T_NK_cells__Naive CD8 T cell" ,
#                                                      "HCC_T_NK_cells__13_CD8_effector memory T",  "HCC_T_NK_cells__2_CD8_MAIT",  "HCC_T_NK_cells__11_CD8_inter-states T","HCC_T_NK_cells__14_CD8_TR_memory T","HCC_T_NK_cells__1_CD8_cytotoxic T" ,"HCC_T_NK_cells__51_CD8_other T",
#                                                      "STAD_T_NK_cells__CD8 naive T" ,"STAD_T_NK_cells__CD8 effector T"),
#                          colnames(output_ove) %in% c("TICA_T_NK_cells__Cytotoxic CD8 T cells","TICA_T_NK_cells__Terminally exhausted CD8 T cells","TICA_T_NK_cells__Pre-exhausted CD8 T cells", "TICA_T_NK_cells__Effector memory CD8 T cells",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c07.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD8.c05.Tem.CXCR5" ,"pancancer_T_cells_CD8_CD4__CD8.c04.Tm.CD52", "pancancer_T_cells_CD8_CD4__CD8.c13.Tex.myl12a",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c15.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD8.c01.Tn.MAL", "pancancer_T_cells_CD8_CD4__CD8.c17.Tm.NME1",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c12.Tex.CXCL13", "pancancer_T_cells_CD8_CD4__CD8.c02.Tm.IL7R", "pancancer_T_cells_CD8_CD4__CD8.c10.Trm.ZNF683", "pancancer_T_cells_CD8_CD4__CD8.c03.Tm.RPS12" ,
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c14.Tex.TCF7", "pancancer_T_cells_CD8_CD4__CD8.c06.Tem.GZMK", "pancancer_T_cells_CD8_CD4__CD8.c11.Tex.PDCD1", "pancancer_T_cells_CD8_CD4__CD8.c16.MAIT.SLC4A10" ,
#                                                      "pancancer_blueprint_T_NK_cells__CD8 memory T", "pancancer_blueprint_T_NK_cells__CD8 exhausted cytotoxic T","pancancer_blueprint_T_NK_cells__CD8 pre-effector T",  "pancancer_blueprint_T_NK_cells__CD8 effector T",
#                                                      "ABTC_T_NK_cells__Exhausted CD8 T cell","ABTC_T_NK_cells__Memory CD8 T cell", "ABTC_T_NK_cells__Effector CD8 T cell" , "ABTC_T_NK_cells__Naive CD8 T cell" ,
#                                                      "HCC_T_NK_cells__13_CD8_effector memory T",  "HCC_T_NK_cells__2_CD8_MAIT",  "HCC_T_NK_cells__11_CD8_inter-states T","HCC_T_NK_cells__14_CD8_TR_memory T","HCC_T_NK_cells__1_CD8_cytotoxic T" ,"HCC_T_NK_cells__51_CD8_other T",
#                                                      "STAD_T_NK_cells__CD8 naive T" ,"STAD_T_NK_cells__CD8 effector T")]
# dim(output_ove)
# ##### NK
# # output_ove <- output_ove[rownames(output_ove) %in% c("ABTC_T_NK_cells__NCAM1pFCGR3An NK", "HCC_T_NK_cells__8_NK", "pancancer_blueprint_T_NK_cells__NK_XCL1", "STAD_T_NK_cells__LNK1", "ABTC_T_NK_cells__NCAM1pFCGR3Ap NK", "ABTC_T_NK_cells__NKT and NK cell", "HCC_T_NK_cells__28_NK", "TICA_T_NK_cells__NK", "pancancer_blueprint_T_NK_cells__Cytotoxic NK", "STAD_T_NK_cells__LNK2", "ABTC_T_NK_cells__DNT", "HCC_T_NK_cells__18_NK", "STAD_T_NK_cells__LNK3", "STAD_T_NK_cells__LNK4", "HCC_T_NK_cells__2_NK", "pancancer_T_cells_CD8_CD4__CD8.c09.Tk.KIR2DL4", "pancancer_T_cells_CD8_CD4__CD8.c08.Tk.TYROBP"),
# #                          colnames(output_ove) %in% c("ABTC_T_NK_cells__NCAM1pFCGR3An NK", "HCC_T_NK_cells__8_NK", "pancancer_blueprint_T_NK_cells__NK_XCL1", "STAD_T_NK_cells__LNK1", "ABTC_T_NK_cells__NCAM1pFCGR3Ap NK", "ABTC_T_NK_cells__NKT and NK cell", "HCC_T_NK_cells__28_NK", "TICA_T_NK_cells__NK", "pancancer_blueprint_T_NK_cells__Cytotoxic NK", "STAD_T_NK_cells__LNK2", "ABTC_T_NK_cells__DNT", "HCC_T_NK_cells__18_NK", "STAD_T_NK_cells__LNK3", "STAD_T_NK_cells__LNK4", "HCC_T_NK_cells__2_NK", "pancancer_T_cells_CD8_CD4__CD8.c09.Tk.KIR2DL4", "pancancer_T_cells_CD8_CD4__CD8.c08.Tk.TYROBP")]
# output_ove <- output_ove[rownames(output_ove) %in% c("TICA_T_NK_cells__NK",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c08.Tk.TYROBP" ,"pancancer_T_cells_CD8_CD4__CD8.c09.Tk.KIR2DL4",
#                                                      "pancancer_blueprint_T_NK_cells__NK_XCL1", "pancancer_blueprint_T_NK_cells__Cytotoxic NK",
#                                                      "ABTC_T_NK_cells__DNT", "ABTC_T_NK_cells__NCAM1nFCGR3Ap NK", "ABTC_T_NK_cells__NKT and NK cell",      "ABTC_T_NK_cells__NCAM1pFCGR3Ap NK"   ,"ABTC_T_NK_cells__NCAM1pFCGR3An NK",
#                                                      "HCC_T_NK_cells__2_NK","HCC_T_NK_cells__8_NK", "HCC_T_NK_cells__18_NK","HCC_T_NK_cells__28_NK",
#                                                      "STAD_T_NK_cells__LNK1", "STAD_T_NK_cells__LNK3", "STAD_T_NK_cells__LNK4", "STAD_T_NK_cells__LNK2"),
#                          colnames(output_ove) %in% c("TICA_T_NK_cells__NK",
#                                                      "pancancer_T_cells_CD8_CD4__CD8.c08.Tk.TYROBP" ,"pancancer_T_cells_CD8_CD4__CD8.c09.Tk.KIR2DL4",
#                                                      "pancancer_blueprint_T_NK_cells__NK_XCL1", "pancancer_blueprint_T_NK_cells__Cytotoxic NK",
#                                                      "ABTC_T_NK_cells__DNT", "ABTC_T_NK_cells__NCAM1nFCGR3Ap NK", "ABTC_T_NK_cells__NKT and NK cell",      "ABTC_T_NK_cells__NCAM1pFCGR3Ap NK"   ,"ABTC_T_NK_cells__NCAM1pFCGR3An NK",
#                                                      "HCC_T_NK_cells__2_NK","HCC_T_NK_cells__8_NK", "HCC_T_NK_cells__18_NK","HCC_T_NK_cells__28_NK",
#                                                      "STAD_T_NK_cells__LNK1", "STAD_T_NK_cells__LNK3", "STAD_T_NK_cells__LNK4", "STAD_T_NK_cells__LNK2")]
# #### a sub cluster
# a_sub_cluster <- c("ABTC: Exhausted CD8 T cell", "TICA: Terminally exhausted CD8 T cells", "pan- blueprint: CD8 exhausted cytotoxic T", "pan- blueprint: CD8 memory T", "TICA: Cytotoxic CD8 T cells", "STAD: LNK2", "HCC: 11_CD8_inter-states T", "pan- blueprint: CD8 pre-effector T", "ABTC: Memory CD8 T cell", "ABTC: NCAM1pFCGR3An NK", "HCC: 8_NK", "pan- blueprint: NK_XCL1", "STAD: LNK1", "ABTC: NCAM1pFCGR3Ap NK", "HCC: 28_NK", "ABTC: Effector CD8 T cell", "ABTC: NCAM1nFCGR3Ap NK", "TICA: NK", "pan- blueprint: Cytotoxic NK", "ABTC: NKT and NK cell", "pan- blueprint: CD8 effector T", "HCC: 13_CD8_effector memory T")
# output_ove <- output_ove[rownames(output_ove) %in% a_sub_cluster,
#                          colnames(output_ove) %in% a_sub_cluster]
dim(output_ove)
dd <- as.dist(1-output_ove)
heatmap_output_ove <- pheatmap(output_ove,clustering_distance_cols = dd,clustering_distance_rows = dd,show_colnames = F)
heatmap_order_row <- rownames(output_ove)[heatmap_output_ove$tree_row$order]
heatmap_order_col <- colnames(output_ove)[heatmap_output_ove$tree_col$order]
# colnames(output_ove)
# output_ove <- output_ove[heatmap_order_row,heatmap_order_col]
# heatmap_order_row <- heatmap_order_col <- c( "HCC: 51_CD8_other T", "	STAD: Proliferative T", "	TICA: Proliferative T cells", "	HCC: 32_CD4_other T", "	TICA: T helper cells", "	pan- blueprint: CD4 exhausted effector T", "	HCC: 9_CD4_regulatory T", "	pan- blueprint: CD4 regulatory T", "	STAD: Regulatory T", "	TICA: Regulatory T cells", "	ABTC: Treg", "	STAD: CD4 naive T cell", "	ABTC: Naive CD4 T cell", "	ABTC: Naive CD8 T cell", "	pan- blueprint: CD4 naive T", "	HCC: 20_CD4_naive T", "	TICA: Naive T cells", "	STAD: CD4 helper T", "	pan- blueprint: CD4 memory/effector T", "	HCC: 2_CD8_MAIT", "	HCC: 35_CD4_central memory T", "	ABTC: Exhausted CD8 T cell", "	TICA: Terminally exhausted CD8 T cells", "	pan- blueprint: CD8 exhausted cytotoxic T", "	pan- blueprint: CD8 memory T", "	TICA: Cytotoxic CD8 T cells", "	STAD: LNK2", "	HCC: 11_CD8_inter-states T", "	pan- blueprint: CD8 pre-effector T", "	ABTC: Memory CD8 T cell", "	ABTC: NCAM1pFCGR3An NK", "	HCC: 8_NK", "	pan- blueprint: NK_XCL1", "	STAD: LNK1", "	ABTC: NCAM1pFCGR3Ap NK", "	HCC: 28_NK", "	ABTC: Effector CD8 T cell", "	ABTC: NCAM1nFCGR3Ap NK", "	TICA: NK", "	pan- blueprint: Cytotoxic NK", "	ABTC: NKT and NK cell", "	pan- blueprint: CD8 effector T", "	HCC: 13_CD8_effector memory T", "	HCC: 1_CD8_cytotoxic T", "	STAD: LNK3", "	HCC: 2_NK", "	TICA: Th17 cells", "	STAD: CD8 naive T", "	STAD: CD8 effector T cell", "	TICA: Transitional memory CD4 T cells", "	ABTC: DNT", "	TICA: Recently activated CD4 T cells", "	TICA: Pre-exhausted CD8 T cells", "	HCC: 18_NK", "	STAD: LNK4", "	HCC: 14_CD8_TR_memory T", "	TICA: Naive-memory CD4 T cells", "	TICA: Effector memory CD8 T cells")
# output_ove <- output_ove[heatmap_order_row,heatmap_order_col]
# ##### CD8
# heatmap_order_row <- heatmap_order_col <- c( "pancancer_T_cells_CD8_CD4__CD8.c14.Tex.TCF7", "pancancer_T_cells_CD8_CD4__CD8.c12.Tex.CXCL13", "pancancer_T_cells_CD8_CD4__CD8.c13.Tex.myl12a", "pancancer_T_cells_CD8_CD4__CD8.c11.Tex.PDCD1", "ABTC_T_NK_cells__Exhausted CD8 T cell", "TICA_T_NK_cells__Terminally exhausted CD8 T cells", "pancancer_blueprint_T_NK_cells__CD8 exhausted cytotoxic T", "pancancer_T_cells_CD8_CD4__CD8.c10.Trm.ZNF683", "HCC_T_NK_cells__51_CD8_other T", "STAD_T_NK_cells__CD8 effector T", "TICA_T_NK_cells__Pre-exhausted CD8 T cells", "pancancer_T_cells_CD8_CD4__CD8.c15.ISG.IFIT1", "TICA_T_NK_cells__Effector memory CD8 T cells", "HCC_T_NK_cells__14_CD8_TR_memory T", "pancancer_T_cells_CD8_CD4__CD8.c01.Tn.MAL", "ABTC_T_NK_cells__Naive CD8 T cell", "pancancer_T_cells_CD8_CD4__CD8.c02.Tm.IL7R", "HCC_T_NK_cells__2_CD8_MAIT", "pancancer_T_cells_CD8_CD4__CD8.c07.Temra.CX3CR1", "ABTC_T_NK_cells__Effector CD8 T cell", "pancancer_blueprint_T_NK_cells__CD8 effector T", "HCC_T_NK_cells__13_CD8_effector memory T", "TICA_T_NK_cells__Cytotoxic CD8 T cells", "pancancer_blueprint_T_NK_cells__CD8 memory T", "pancancer_T_cells_CD8_CD4__CD8.c05.Tem.CXCR5", "pancancer_T_cells_CD8_CD4__CD8.c06.Tem.GZMK", "pancancer_blueprint_T_NK_cells__CD8 pre-effector T", "ABTC_T_NK_cells__Memory CD8 T cell", "HCC_T_NK_cells__11_CD8_inter-states T", "HCC_T_NK_cells__1_CD8_cytotoxic T", "STAD_T_NK_cells__CD8 naive T")
# ##### CD4
# # heatmap_order_row <- heatmap_order_col <- c( "pancancer_T_cells_CD8_CD4__CD4.c12.Tem.GZMK", "pancancer_T_cells_CD8_CD4__CD4.c13.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD4.c22.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD4.c21.Treg.OAS1", "STAD_T_NK_cells__Proliferative T", "TICA_T_NK_cells__Proliferative T cells", "HCC_T_NK_cells__32_CD4_other T", "TICA_T_NK_cells__T helper cells", "STAD_T_NK_cells__Regulatory T", "HCC_T_NK_cells__9_CD4_regulatory T", "pancancer_blueprint_T_NK_cells__CD4 regulatory T", "pancancer_T_cells_CD8_CD4__CD4.c20.Treg.TNFRSF9", "TICA_T_NK_cells__Regulatory T cells", "ABTC_T_NK_cells__Treg", "pancancer_T_cells_CD8_CD4__CD4.c18.Treg.RTKN2", "pancancer_T_cells_CD8_CD4__CD4.c19.Treg.S1PR1", "STAD_T_NK_cells__CD4 naive T", "HCC_T_NK_cells__20_CD4_naive T", "HCC_T_NK_cells__35_CD4_central memory T", "pancancer_T_cells_CD8_CD4__CD4.c04.Tn.il7r", "pancancer_T_cells_CD8_CD4__CD4.c03.Tn.ADSL", "pancancer_T_cells_CD8_CD4__CD4.c01.Tn.TCF7", "TICA_T_NK_cells__Naive T cells", "pancancer_blueprint_T_NK_cells__CD4 naive T", "ABTC_T_NK_cells__Naive CD4 T cell", "pancancer_T_cells_CD8_CD4__CD4.c17.TfhTh1.CXCL13", "pancancer_T_cells_CD8_CD4__CD4.c16.Tfh.CXCR5", "pancancer_blueprint_T_NK_cells__CD4 exhausted effector T", "TICA_T_NK_cells__Recently activated CD4 T cells", "pancancer_T_cells_CD8_CD4__CD4.c10.Tm.CAPG", "TICA_T_NK_cells__Naive-memory CD4 T cells", "pancancer_blueprint_T_NK_cells__CD4 memory/effector T", "STAD_T_NK_cells__CD4 helper T", "TICA_T_NK_cells__Transitional memory CD4 T cells", "pancancer_T_cells_CD8_CD4__CD4.c02.Tn.PASK", "pancancer_T_cells_CD8_CD4__CD4.c08.Tm.CREM")  
# heatmap_order_row <- heatmap_order_col <- c( "pancancer_T_cells_CD8_CD4__CD4.c12.Tem.GZMK", "pancancer_T_cells_CD8_CD4__CD4.c13.Temra.CX3CR1", "pancancer_T_cells_CD8_CD4__CD4.c22.ISG.IFIT1", "pancancer_T_cells_CD8_CD4__CD4.c21.Treg.OAS1", "STAD_T_NK_cells__Proliferative T", "TICA_T_NK_cells__Proliferative T cells", "HCC_T_NK_cells__32_CD4_other T", "STAD_T_NK_cells__Regulatory T", "HCC_T_NK_cells__9_CD4_regulatory T", "pancancer_blueprint_T_NK_cells__CD4 regulatory T", "pancancer_T_cells_CD8_CD4__CD4.c20.Treg.TNFRSF9", "TICA_T_NK_cells__Regulatory T cells", "ABTC_T_NK_cells__Treg", "pancancer_T_cells_CD8_CD4__CD4.c18.Treg.RTKN2", "pancancer_T_cells_CD8_CD4__CD4.c19.Treg.S1PR1", "pancancer_T_cells_CD8_CD4__CD4.c04.Tn.il7r", "pancancer_T_cells_CD8_CD4__CD4.c03.Tn.ADSL", "pancancer_T_cells_CD8_CD4__CD4.c01.Tn.TCF7", "TICA_T_NK_cells__Naive T cells", "pancancer_blueprint_T_NK_cells__CD4 naive T", "ABTC_T_NK_cells__Naive CD4 T cell", "TICA_T_NK_cells__Naive-memory CD4 T cells", "STAD_T_NK_cells__CD4 naive T", "HCC_T_NK_cells__20_CD4_naive T", "HCC_T_NK_cells__35_CD4_central memory T", "pancancer_T_cells_CD8_CD4__CD4.c17.TfhTh1.CXCL13", "pancancer_T_cells_CD8_CD4__CD4.c16.Tfh.CXCR5", "pancancer_blueprint_T_NK_cells__CD4 exhausted effector T", "TICA_T_NK_cells__T helper cells", "STAD_T_NK_cells__CD4 helper T", "pancancer_T_cells_CD8_CD4__CD4.c06.Tm.ANXA1", "pancancer_T_cells_CD8_CD4__CD4.c14.Th17.SLC4A10", "pancancer_T_cells_CD8_CD4__CD4.c15.Th17.IL23R",  "pancancer_blueprint_T_NK_cells__CD4 memory/effector T", "TICA_T_NK_cells__Th17 cells", "pancancer_T_cells_CD8_CD4__CD4.c08.Tm.CREM", "TICA_T_NK_cells__Transitional memory CD4 T cells", "pancancer_T_cells_CD8_CD4__CD4.c02.Tn.PASK", "pancancer_T_cells_CD8_CD4__CD4.c24.Mix.NME2", "pancancer_T_cells_CD8_CD4__CD4.c23.Mix.NME1", "TICA_T_NK_cells__Recently activated CD4 T cells", "pancancer_T_cells_CD8_CD4__CD4.c05.Tm.TNF", "pancancer_T_cells_CD8_CD4__CD4.c09.Tm.CCL5","pancancer_T_cells_CD8_CD4__CD4.c07.Tm.ANXA2", "pancancer_T_cells_CD8_CD4__CD4.c11.Tm.GZMA", "pancancer_T_cells_CD8_CD4__CD4.c10.Tm.CAPG")  
# 
# ##### NK
# heatmap_order_row <- heatmap_order_col <- c( "ABTC_T_NK_cells__DNT", "STAD_T_NK_cells__LNK2", "pancancer_T_cells_CD8_CD4__CD8.c09.Tk.KIR2DL4", "pancancer_T_cells_CD8_CD4__CD8.c08.Tk.TYROBP", "ABTC_T_NK_cells__NCAM1pFCGR3An NK", "HCC_T_NK_cells__8_NK", "pancancer_blueprint_T_NK_cells__NK_XCL1", "STAD_T_NK_cells__LNK1", "ABTC_T_NK_cells__NCAM1pFCGR3Ap NK", "HCC_T_NK_cells__18_NK", "STAD_T_NK_cells__LNK3", "STAD_T_NK_cells__LNK4", "HCC_T_NK_cells__28_NK", "TICA_T_NK_cells__NK", "pancancer_blueprint_T_NK_cells__Cytotoxic NK", "ABTC_T_NK_cells__NKT and NK cell", "HCC_T_NK_cells__2_NK")  


heatmap_order_row <- heatmap_order_col <- c( "STAD: Proliferative T", "TICA: Proliferative T cells", "HCC: 32_CD4_other T", "HCC: 51_CD8_other T", 
                                             
                                             "HCC: 9_CD4_regulatory T", "pan- blueprint: CD4 regulatory T", "STAD: Regulatory T", "TICA: Regulatory T cells", "ABTC: Treg", 
                                             
                                             "STAD: CD4 naive T cell", "ABTC: Naive CD4 T cell",  "pan- blueprint: CD4 naive T", "HCC: 20_CD4_naive T", "TICA: Naive T cells",  "HCC: 35_CD4_central memory T",
                                             
                                             "pan- blueprint: CD4 memory/effector T",  "TICA: Naive-memory CD4 T cells", 
                                             
                                             "STAD: CD8 effector T cell", "TICA: Transitional memory CD4 T cells", "ABTC: DNT", 
                                            
                                             "TICA: T helper cells", "STAD: CD4 helper T","pan- blueprint: CD4 exhausted effector T", 
                                             
                                             
                                             
                                             "ABTC: Exhausted CD8 T cell", "pan- blueprint: CD8 exhausted cytotoxic T",  "TICA: Terminally exhausted CD8 T cells", "STAD: LNK2", "TICA: Pre-exhausted CD8 T cells",
                                             
                                             "pan- blueprint: CD8 effector T", "HCC: 13_CD8_effector memory T", "ABTC: NKT and NK cell", "STAD: LNK3","HCC: 28_NK", "ABTC: Effector CD8 T cell",
                                             
                                             "TICA: Cytotoxic CD8 T cells", "HCC: 11_CD8_inter-states T", "pan- blueprint: CD8 pre-effector T", "ABTC: Memory CD8 T cell", "HCC: 1_CD8_cytotoxic T", 
                                             
                                             "pan- blueprint: CD8 memory T", "TICA: Effector memory CD8 T cells","TICA: Recently activated CD4 T cells", 
                                             
                                             "ABTC: Naive CD8 T cell","HCC: 2_CD8_MAIT", "STAD: CD8 naive T","HCC: 14_CD8_TR_memory T","TICA: Th17 cells", 
                                             
                                             
                                             
                                             
                                             
                                             
                                             
                                             
                                            
                                             "ABTC: NCAM1pFCGR3An NK", "HCC: 8_NK", "pan- blueprint: NK_XCL1", 
                                             
                                             # "ABTC: NCAM1pFCGR3Ap NK", 
                                             
                                             "ABTC: NCAM1pFCGR3Ap NK", "HCC: 2_NK","STAD: LNK1","TICA: NK", "pan- blueprint: Cytotoxic NK",
                                             
                                             "ABTC: NCAM1nFCGR3Ap NK",   #"HCC: 2_NK",#"STAD: LNK3",
                                             
                                             "HCC: 18_NK",  "STAD: LNK4"
                                             
                                            
                                             
                                             
                                             
                                             )
output_ove <- output_ove[heatmap_order_row,heatmap_order_col]
pheatmap(output_ove,
         # clustering_distance_cols = dd,clustering_distance_rows = dd,
         cluster_rows = F,cluster_cols = F,
         show_colnames = F,fontsize = 6)
write.table(output_ove,file=paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_ove_plot_ordered_20240214.txt"),quote = F,sep="\t")
## Then go to excel

# output_cor <- readRDS(paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_cor_plot_20230606.rds"))
# # output_cor <- output_cor[c(1:32,34:58),c(1:32,34:58)]
# # a_sub_cluster <- c("ABTC: Exhausted CD8 T cell", "TICA: Terminally exhausted CD8 T cells", "pan- blueprint: CD8 exhausted cytotoxic T", "pan- blueprint: CD8 memory T", "TICA: Cytotoxic CD8 T cells", "STAD: LNK2", "HCC: 11_CD8_inter-states T", "pan- blueprint: CD8 pre-effector T", "ABTC: Memory CD8 T cell", "HCC: 28_NK", "ABTC: Effector CD8 T cell", "ABTC: NCAM1nFCGR3Ap NK", "TICA: NK", "pan- blueprint: Cytotoxic NK", "ABTC: NKT and NK cell", "pan- blueprint: CD8 effector T", "HCC: 13_CD8_effector memory T", "ABTC: NCAM1pFCGR3An NK", "HCC: 8_NK", "pan- blueprint: NK_XCL1", "STAD: LNK1", "HCC: 1_CD8_cytotoxic T", "TICA: Th17 cells", "STAD: CD8 naive T")
# # output_cor <- output_cor[rownames(output_cor) %in% a_sub_cluster,
# #                          colnames(output_cor) %in% a_sub_cluster];dim(output_cor)
# output_cor <- output_cor[heatmap_order_row,heatmap_order_col]
# dd <- as.dist(1-output_cor)
# pheatmap(output_cor,
#          clustering_distance_cols = dd,clustering_distance_rows = dd,
#          cluster_rows = T,cluster_cols = T,
#          show_colnames = F,fontsize = 6)
# write.table(output_cor,file=paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_cor_asubcluster_plot_ordered_20240214.txt"),quote = F,sep="\t")
# 
# 
# output_gsea_es <- readRDS(paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_gese_es_plot_20230606.rds"))
# rownames(output_gsea_es)[which(rownames(output_gsea_es)=="HCC_T_NK_cells: 14_TR_memory T")] <- "HCC_T_NK_cells: 14_CD8_TR_memory T"
# colnames(output_gsea_es)[which(colnames(output_gsea_es)=="HCC_T_NK_cells: 14_TR_memory T")] <- "HCC_T_NK_cells: 14_CD8_TR_memory T"
# rownames(output_gsea_es)[which(rownames(output_gsea_es)=="HCC_T_NK_cells: 35_central memory T")] <- "HCC_T_NK_cells: 35_CD4_central memory T"
# colnames(output_gsea_es)[which(colnames(output_gsea_es)=="HCC_T_NK_cells: 35_central memory T")] <- "HCC_T_NK_cells: 35_CD4_central memory T"
# output_gsea_es <- output_gsea_es[heatmap_order_row,heatmap_order_col]
# dd <- as.dist(1-output_gsea_es)
# pheatmap(output_gsea_es,
#          clustering_distance_cols = dd,clustering_distance_rows = dd,
#          cluster_rows = F,cluster_cols = F,
#          show_colnames = F,fontsize = 6)
# write.table(output_gsea_es,file=paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_gsea_es_plot_ordered_20240214.txt"),quote = F,sep="\t")



#########################################
#### Jing Yang
#### 2024.02.14 Figure 1 part 2. marker gene heatmap
#########################################
rm(list=ls())
library(purrr)
library(tidyverse)
library(scales)
library(Seurat)
library(infotheo)
library(pheatmap)
library(limma)
library(ggpubr)
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
                              "HCC: 28_T/NK",                           "HCC: 51_T/NK",                           
                              
                              
                              "STAD: LNK1",                            "STAD: LNK3",                             "STAD: LNK4",                             "STAD: LT2",
                              "STAD: LT3",                              "STAD: LNK2",                             "STAD: LT4",
                              "STAD: LT5",                              "STAD: LT1",                              "STAD: LPT")
workdir <- "/Users/jingyang/Dropbox/work/tumor immune cell atlas estimation/data/"
# TICA_markers <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells_markers_downsample.rds"))
# pancancer_T_cells_CD8_CD4 <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_markers.rds"))

TICA_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_TICA_T_NK_cells_markers_downsample.rds"))
TICA_T_NK_cells_markers$cluster <- paste0("TICA_T_NK_cells__",TICA_T_NK_cells_markers$cluster);dim(TICA_T_NK_cells_markers)

# pancancer_T_cells_CD8_CD4_markers <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_markers.rds"))
# pancancer_T_cells_CD8_CD4_markers <- readRDS(paste0(workdir,"atlas_pancancer_T_cells_CD8_CD4_markers_filled.rds"))
# pancancer_T_cells_CD8_CD4_markers$cluster <- paste0("pancancer_T_cells_CD8_CD4__",pancancer_T_cells_CD8_CD4_markers$cluster);dim(pancancer_T_cells_CD8_CD4_markers)

pancancer_blueprint_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_pancancer_blueprint_T_NK_cells_markers.rds"))
pancancer_blueprint_T_NK_cells_markers$cluster <- paste0("pancancer_blueprint_T_NK_cells__",pancancer_blueprint_T_NK_cells_markers$cluster);dim(pancancer_blueprint_T_NK_cells_markers)


ABTC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_ABTC_T_NK_cells_markers.rds"))
ABTC_T_NK_cells_markers$cluster <- paste0("ABTC_T_NK_cells__",ABTC_T_NK_cells_markers$cluster);dim(ABTC_T_NK_cells_markers)

HCC_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_HCC_T_NK_cells_markers.rds"))
HCC_T_NK_cells_markers$cluster <- paste0("HCC_T_NK_cells__",HCC_T_NK_cells_markers$cluster);dim(HCC_T_NK_cells_markers)

STAD_T_NK_cells_markers <- readRDS(paste0(workdir,"atlas_STAD_T_NK_cells_markers.rds"))
STAD_T_NK_cells_markers$cluster <- paste0("STAD_T_NK_cells__",STAD_T_NK_cells_markers$cluster);dim(STAD_T_NK_cells_markers)

markers <- rbind(TICA_T_NK_cells_markers,pancancer_blueprint_T_NK_cells_markers,ABTC_T_NK_cells_markers,HCC_T_NK_cells_markers,STAD_T_NK_cells_markers);dim(markers)



Mtpattern= "^MT-|^mt-"
rRNApattern="^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]"
rRNA.genes <- grep(pattern = rRNApattern,  markers$gene, value = TRUE)
Mt.genes<- grep (pattern= Mtpattern,markers$gene, value=TRUE ) 
keep.genes <- dplyr::setdiff(markers$gene, c(rRNA.genes,Mt.genes))
markers <- markers[markers$gene %in% keep.genes,];dim(markers)

# markers <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.25,];dim(markers)
# markers <- markers[order(-markers$avg_log2FC),];dim(markers)
# tmp1 <- lapply(markers %>% split(.$cluster),function(x) x[1:ceiling(nrow(x)/10),])
# tmp1 <- lapply(markers %>% split(.$cluster),function(x) x[1:5,])#[1:62,])
# output_tmp <- NULL
# for (i in 1:length(tmp1)) {
#   output_tmp <- rbind(output_tmp,tmp1[[i]])
# };dim(output_tmp);
# genes <- c(unique(output_tmp$gene),"CD4");length(genes)

genes <- c(
  # "STMN1","MKI67", #"TOP2A","CDK1",#proliferative
  # # "TYMS",#"CXCL13",
  # "FOXP3","TNFRSF4",#"LTB","BATF", #reg
  # "CCR7","IL7R",#"TCF7",#"SELL",#"LEF1",#"IL7R", #naive
  # "MALAT1",#"GZMA","GZMH", #STAD:CD8 effector T(mapping to effector and naïve CD8 T
  # # cells (CD8A, GZMH, GZMM, and NKG7; LT1–2),)/TICA: transitional memeory CD4 T/ABTC: DNT
  # "IL17A","IL17F",#Th17
  # "PDCD1","HAVCR2","IFIT3",#"CXCL13",#"CCL5",#,#"LAG3",#exhausted CD8 T
  # #"GZMA","CCL5",#"GZMB","GZMH","NKG7","CXCR4",
  # "GZMB","PRF1", # effector CD8 T
  # "GZMA","GZMH", #"GZMK",#pre-effector
  # # "IFIT3","HAVCR2",
  # "KLRB1",# marker of NKT "ZBTB16", "KLRD1",
  # "TRAV1-2", # marker of MAIT
  # "ZNF683", # memeory CD8 T
  #  #"TH17 
  # # "TBX21","IFNG","TNF",#Th1
  # # "IL4","IL5",#Th2
  # 
  # "KLRC1","NCR1","NCAM1",#"NCR2","NCR3","KLRB1","NCAM1",#"NCAM1",#"SPTSSB","XCL1", #NK
  # "FCGR3A",#NK
  "CD3D","CD3E","CD3G")#,
  # "KIT","RORC", #ILC markers
  # "TRDC","TRGC1",#delta gamma T cells
  #"NKG7",
  #"FGFBP2",#"CD69",#"FCGR3A","GNLY","GZMH","FGFBP2", #NK
  # "CXCR4","SLC7A5",
  #"NEAT1",#"RNF213",
  # "RHOA",
  # "PSAP",
  # "EIF3H",
  # "MX1",
  #"CCL4","HSPA1B",#"HSPA6",
  # "EEF1G",
  # "FXYD5",
  #"CD69",
  # "FOS","JUNB",
  #"IL17A",
  # "CD4","CD8A","CD8B")




# "GZMK","GZMB","PRF1","CD44","CD69","FAS","FASLG","PDCD1","LAG3","CTLA4","FGFBP2","GZMH","GNLY","NR4A1","BAG3","HSPA1A","IFIT1","MX1","DKK3",
#            "CCR4","EOMES","CNN2","LIMD2","CD27","LGR6","KLRC4","CD244","SEMA4A","ITGA1","KLRB1","PRDM1","CCR7","SELL","TCF7")
# markers <- output_tmp

# correlations
# charcters <- c("atlas_TICA_T_NK_cells.rds","atlas_pancancer_T_cells_CD8_CD4.rds","atlas_pancancer_blueprint_T_NK_cells",
#                "atlas_ABTC_T_NK_cells.rds","atlas_HCC_T_NK_cells.rds","atlas_STAD_T_NK_cells.rds")
charcters <- c("TICA_T_NK_cells","pancancer_blueprint_T_NK_cells",
               "ABTC_T_NK_cells","HCC_T_NK_cells","STAD_T_NK_cells")
add_index <- c("TICA: ","pan- blueprint: ", "ABTC: ", "HCC: ","STAD: ")
output_plot <- NULL
for (i in 1:length(charcters)) {
  print (i)
  atlas <- readRDS(paste0(workdir,"atlas_",charcters[i],".rds"))
  # if (i ==2 ) {
  #   Idents(atlas) <- atlas$cell_types_from_author
  #   ident <- as.character(unique(atlas@active.ident));ident
  #   output <- NULL
  #   for (m in 1:length(ident)) {
  #     # print (m)
  #     tmp <- as.matrix(GetAssayData(subset(atlas,idents=ident[m]),assay = "originalexp", slot="data"))
  #     for (n in 1:nrow(tmp)) {
  #       tmp[n,is.nan(tmp[n,])] <- mean(tmp[n,],na.rm =T)
  #     }
  #     output <- cbind(output,tmp)
  #   }
  #   output <- output[,colnames(atlas@assays$originalexp@data)]
  #   atlas@assays$originalexp@data <- output
  # } else {
  Idents(atlas) <- atlas$cell_types_from_author
  # }
  # if ( i == 2) {
  #   DefaultAssay(atlas) <- "originalexp"
  # } else {
  DefaultAssay(atlas) <- "RNA"
  # }
  
  output_pre <- DotPlot(atlas,features = genes,scale=FALSE)$data
  output_pre$id <- paste0(add_index[i],output_pre$id)
  output_pre <- data.frame(output_pre,source=add_index[i])
  output_plot <- rbind(output_plot,output_pre)
  rm("atlas")
  
}


# saveRDS(output_plot,file=paste0(workdir,"1_cluster_correlation_within_reference_atlases_clusterexpr_median_filled_with_pct_20240214.rds"))
# output_plot <- data.frame(gene=rownames(output_plot),output_plot)
# output_plot_sub <- output_plot[output_plot$features.plot %in% c("CD4","CD8A","CD8B"),];dim(output_plot_sub)
output_plot$id <- as.character(new_cluster_names[output_plot$id])

## manually ordered cell types
cluster_orders <- read.table(paste0(workdir,"1_cluster_correlation_within_reference_atlases_output_ove_plot_ordered_20240214.txt"),header=T,sep="\t",row.names = 1)
output_plot$id <- factor(output_plot$id,levels =rev(rownames(cluster_orders)))
output_plot$features.plot <- factor(output_plot$features.plot,levels=genes)
output_plot <- output_plot[!is.na(output_plot$id),]


pdf(paste0(workdir,"1_cluster_correlation_within_reference_atlases_clusterexpr_median_filled_with_pct_20250506.pdf"),height = 21, width = 15)
ggplot(data=output_plot,mapping=aes_string(y='id',x='features.plot',color='avg.exp.scaled')) + 
  geom_point(mapping=aes_string(size='pct.exp')) +
  scale_color_gradient(low="blue",high="red") +
  theme_classic() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(size = 28,color = "black",angle = 90,vjust=0.5,hjust=1),
        axis.text.y = element_text(size = 10,color = "black"),
        # axis.text.y = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_text(size=28)) + 
  # guides(color = guide_legend(override.aes = list(size=14), ncol=1) )+ 
  labs(y="")
dev.off()




