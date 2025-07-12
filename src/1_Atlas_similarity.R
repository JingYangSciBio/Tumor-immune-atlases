#########################################
#### Jing Yang
#### 2024.02.14 Figure 2A. #reorder plots on output_cor, output_ove and output_gsea_es 
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




