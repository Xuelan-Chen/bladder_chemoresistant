# Proteomics analysis

## Part1. we should load the packages we need for following analysis

~~~R
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(Matrix)
  library(proxy)
  library(gplots)
  library(Rtsne)
  library(densityClust)
  library(irlba)
  library(monocle)
  library(plyr)
  library(DOSE)
  library(clusterProfiler)
  library(topGO)
  library(pathview)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(cowplot)
  library(ggplot2)
  library(VGAM)
  library(gtools)
  library(RSCORE)
  library(GSVA)
library(GSEABase)
library(GSVAdata)
data(c2BroadSets)
c2BroadSets
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(edgeR)
library(BisqueRNA)
library(RColorBrewer)
library(trqwe)
library(ggpubr)
require(ComplexHeatmap)
require(BuenColors)
require(scales)
library(nichenetr)
library(tidyverse)
library(iTALK)
library(MuSiC)
library(tidyr)
library(DESeq2)
library("survival")
library("survminer")
library(pheatmap)
library(corrplot)
})
source("./MyBestFunction_scRNA.R")
library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
~~~

## Part2. we should load the data we need for following analysis

~~~R
Proteomics_data <- mcreadRDS("./Proteomics/RData_files/Proteomics_Sen_Ren_E64_all_result.rds",mc.cores=20)
ALL_GSEA_GMT <- read.gmt("./msigdb.v7.1.symbols.gmt")
common_93_gene <- read.csv(file="./Bulk_ATACseq/RData_files/1_common93_ATAC_RNA_sc_Up.csv")
RNAseq_result_E64 <- read.csv(file="./Bulk_RNAseq/csv_files/2_RNAseq_E64_count_tpm_all_result.csv")
~~~

## Part3.  the detail codes of visualization of Proteomics data

~~~R
ALL_GSEA_GMT$ont <- as.character(ALL_GSEA_GMT$ont)
EPITHELIAL <- unique(grep("*EPITHELIAL_CELL_DIFFERENTIATION",ALL_GSEA_GMT$ont,value=TRUE))
KERATIN <- unique(grep("*KERATIN",ALL_GSEA_GMT$ont,value=TRUE))
SQUAMOUS <- unique(grep("*SQUAM",ALL_GSEA_GMT$ont,value=TRUE))
all_pathway <- c(EPITHELIAL,KERATIN,SQUAMOUS)
KERATIN_GSEA <- XY_subset(ALL_GSEA_GMT,"ont",all_pathway)

Proteomics_data_new <- Proteomics_data %>% mutate(human_symbol = convert_mouse_to_human_symbols(as.character(Proteomics_data$Gene.name))) %>% drop_na()
Proteomics_data_new <- Proteomics_data_new[order(Proteomics_data_new$Ren_vs_Sen_log_FC,decreasing=TRUE),]
hsg_genelist <- Proteomics_data_new$Ren_vs_Sen_log_FC
names(hsg_genelist) <- Proteomics_data_new$human_symbol
hsg_egmt2 <- GSEA(hsg_genelist, TERM2GENE=KERATIN_GSEA, verbose=TRUE,minGSSize=2,pvalueCutoff = 0.05)
aa <- jdb_palette("brewer_heat",type = "continuous")[1:length(jdb_palette("brewer_heat",type = "continuous"))]
tmp_files <- XY_ridgeplot.gseaResult(hsg_egmt2,fill="NES", core_enrichment = TRUE)
tmp_files1 <- XY_ridgeplot.gseaResult(hsg_egmt2,fill="pvalue", core_enrichment = TRUE)
tmp_files$pvalue <- tmp_files1$pvalue
tmp_files <- subset(tmp_files,category!="GO_POSITIVE_REGULATION_OF_KERATINOCYTE_PROLIFERATION")
tmp_files <- subset(tmp_files,category!="GO_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION")

ggplot(tmp_files, aes_string(x="value", y="category", fill="NES")) + ggridges ::geom_density_ridges() +
    scale_fill_gradientn(name = "NES", colors=aa, guide=guide_colorbar(reverse=TRUE)) + xlim(0,2.5)+
    xlab(NULL) + ylab(NULL) +  theme_dose()+labs(title=c(" proteomics "))
~~~

![image-20211224214253070](Proteomics%20analysis.assets/image-20211224214253070.png)

~~~R
proteomics_data_common93 <- XY_subset(Proteomics_data,"Gene.name",common_93_gene$gene)
proteomics_data_common93 <- proteomics_data_common93[order(proteomics_data_common93$Ren_vs_Sen_pvalue,decreasing=FALSE),]
rownames(proteomics_data_common93) <- proteomics_data_common93$Gene.name
used_data <- proteomics_data_common93[,c(1:5)]
data_plot_zscore <- t(apply(used_data, 1, function(x) (x-mean(x))/sd(x)))
data_plot_zscore <- na.omit(data_plot_zscore)
break1=seq(min(data_plot_zscore), 0, length.out=25)
break2=seq(0,max(data_plot_zscore), length.out=26)
the_break<-c(break1,break2[-1])
pheatmap(data_plot_zscore,show_rownames=TRUE,border_color=NA,cluster_cols=F,cluster_row=FALSE,cellwidth=30,cellheight=6,fontsize=7,color = colorRampPalette(c("navy", "white","firebrick3"))(50),breaks=the_break)
~~~

![image-20211224214339069](Proteomics%20analysis.assets/image-20211224214339069.png)

~~~R
JJ_result <- XY_subset(Proteomics_data,"Gene.name",c("Epcam","Krt5","Krt14","Krt6a","Krt17","Krt16","Krt6b","Col17a1","Dsc3","Tgm1","Trp63","Sprr1a","Sprr1b","Dsc2","Dsc1","Tgm3"))
rownames(JJ_result) <- JJ_result$Gene.name
JJ_result_heatmap <- JJ_result[,c(1:8)]
heatmap_resdeal_all_1_zscore <- t(apply(JJ_result_heatmap, 1, function(x) (x-mean(x))/sd(x)))
heatmap_resdeal_all_1_zscore <- na.omit(heatmap_resdeal_all_1_zscore)
heatmap_resdeal_all_1_zscore <- as.data.frame(heatmap_resdeal_all_1_zscore)
heatmap_resdeal_all_1_zscore$mean_E <- apply(heatmap_resdeal_all_1_zscore[,c(6:8)],1,mean)
heatmap_resdeal_all_1_zscore$mean_P <- apply(heatmap_resdeal_all_1_zscore[,c(1:3)],1,mean)
heatmap_resdeal_all_1_zscore$Fc <- heatmap_resdeal_all_1_zscore$mean_E-heatmap_resdeal_all_1_zscore$mean_P
heatmap_resdeal_all_1_zscore <- heatmap_resdeal_all_1_zscore[order(heatmap_resdeal_all_1_zscore$Fc),]
heatmap_resdeal_all_1_zscore <- heatmap_resdeal_all_1_zscore[,c(1:8)]
heatmap_resdeal_all_1_zscore <- as.matrix(heatmap_resdeal_all_1_zscore)
data <- ifelse(heatmap_resdeal_all_1_zscore >1,1,heatmap_resdeal_all_1_zscore)
break1=seq(-1,-0.5, length.out=50)
break2=seq(-0.5, 1, length.out=51)
the_break<-c(break1,break2[-1])
pheatmap(data,border_color=NA,cluster_rows=F,cluster_cols=F,breaks=the_break,color = colorRampPalette(jdb_palette("brewer_celsius"))(100),show_rownames=TRUE,main="Proteomics semi-squamous")
~~~

![image-20211224214412549](Proteomics%20analysis.assets/image-20211224214412549.png)

~~~R
JJ_result <- XY_subset(Proteomics_data,"Gene.name",c("Krt5","Krt14","Krt6a","Krt17","Krt16","Krt6b","Col17a1","Dsc3","Tgm1","Dsc1","Dsc2","Trp63","Sprr1a","Sprr1b"))
rownames(JJ_result) <- JJ_result$Gene.name
JJ_result_heatmap <- JJ_result[,c(1:5)]
heatmap_resdeal_all_1_zscore <- t(apply(JJ_result_heatmap, 1, function(x) (x-mean(x))/sd(x)))
heatmap_resdeal_all_1_zscore <- na.omit(heatmap_resdeal_all_1_zscore)
library(pheatmap)
break1=seq(min(heatmap_resdeal_all_1_zscore), 0, length.out=25)
break2=seq(0, max(heatmap_resdeal_all_1_zscore), length.out=26)
the_break<-c(break1,break2[-1])
pheatmap(heatmap_resdeal_all_1_zscore,border_color=NA,cluster_rows=F,cluster_cols=F,cellwidth=15,fontsize=10,cellheight=20, legend_labels = c("-1.5", "0", "1.5"),color = colorRampPalette(c("navy", "white","firebrick3"))(50),show_rownames=TRUE)
~~~

![image-20211224214440386](Proteomics%20analysis.assets/image-20211224214440386.png)

~~~R
E64_vs_Vehicle_Up$entrez <- mapIds(org.Mm.eg.db,
                     keys=as.character(E64_vs_Vehicle_Up$Gene.name),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
ee  <-as.matrix(E64_vs_Vehicle_Up$entrez)
  dd <- as.vector(ee)
GOE64_vs_Vehicle_Up_all <- enrichGO(gene = dd, 
             OrgDb = org.Mm.eg.db,
        ont = "all", 
                 pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)

GO_up_RNA <- as.data.frame(GOE64_vs_Vehicle_Up_all)
GO_up_RNA_BP <- subset(GO_up_RNA,ONTOLOGY=="BP")
GO_up_RNA_BP$log10_Pvalue <- -log(GO_up_RNA_BP$pvalue,10)
GO_up_RNA_BP_5 <- head(GO_up_RNA_BP,5)
ggbarplot(GO_up_RNA_BP_5, 
  x = "Description", 
  y = "log10_Pvalue",
  color = "#FA6A63",            # Set bar border colors to white
  fill ="#FA6A63",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="Proteomics E64 Up GO BP")
~~~

![image-20211224214628806](Proteomics%20analysis.assets/image-20211224214628806.png)

~~~R
JJ_result <- XY_subset(Proteomics_data,"Gene.name",c("Gsdmc","Gsdmd","Gsdma","Casp8"))
rownames(JJ_result) <- JJ_result$Gene.name
JJ_result_heatmap <- JJ_result[,c(4:8)]
heatmap_resdeal_all_1_zscore <- t(apply(JJ_result_heatmap, 1, function(x) (x-mean(x))/sd(x)))
heatmap_resdeal_all_1_zscore <- na.omit(heatmap_resdeal_all_1_zscore)
library(pheatmap)
break1=seq(min(heatmap_resdeal_all_1_zscore), 0, length.out=25)
break2=seq(0, max(heatmap_resdeal_all_1_zscore), length.out=26)
the_break<-c(break1,break2[-1])
pheatmap(heatmap_resdeal_all_1_zscore,border_color=NA,cluster_rows=F,cluster_cols=F,cellwidth=15,fontsize=10,cellheight=20, legend_labels = c("-1.5", "0", "1.5"),color = colorRampPalette(c("navy", "white","firebrick3"))(50),show_rownames=TRUE)
~~~

![image-20211224214721957](Proteomics%20analysis.assets/image-20211224214721957.png)

~~~R
ALL_GSEA_GMT$ont <- as.character(ALL_GSEA_GMT$ont)
TNFA_genes  <- XY_subset(ALL_GSEA_GMT,"ont",c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"))
TNFA_genes$gene <- as.character(TNFA_genes$gene)
TNFA_genes_pathway = TNFA_genes %>% mutate(mouse_gene = convert_human_to_mouse_symbols(TNFA_genes$gene)) %>% drop_na()

Proteomics_TNFA  <- XY_subset(Proteomics_data,"Gene.name",as.character(TNFA_genes_pathway$mouse_gene))
Proteomics_TNFA_Up <- subset(Proteomics_TNFA,E64_vs_Vehicle_log_FC > 0)
RNA_TNFA  <- XY_subset(RNAseq_result_E64,"symbol",as.character(TNFA_genes_pathway$mouse_gene))
RNA_TNFA_Up <- subset(RNA_TNFA,E64_vs_Vehicle_log2FoldChange > 0)
all_common_gene <- intersect(Proteomics_TNFA_Up$Gene.name,RNA_TNFA_Up$symbol)
Proteomics_common  <- XY_subset(Proteomics_TNFA_Up,"Gene.name",all_common_gene)
RNA_common  <- XY_subset(RNA_TNFA_Up,"symbol",all_common_gene)

rownames(Proteomics_common) <- Proteomics_common$Gene.name
Proteomics_common1 <- Proteomics_common[,c(4:8)]
data_plot_zscore <- t(apply(Proteomics_common1, 1, function(x) (x-mean(x))/sd(x)))
data_plot_zscore <- na.omit(data_plot_zscore)
break1=seq(min(data_plot_zscore), 0, length.out=25)
break2=seq(0,max(data_plot_zscore), length.out=26)
the_break<-c(break1,break2[-1])
pheatmap(data_plot_zscore,show_rownames=TRUE,border_color=NA,cluster_cols=F,cluster_row=F,cellwidth=30,cellheight=14,fontsize=7,color = colorRampPalette(c("navy", "white","firebrick3"))(50),breaks=the_break)
~~~

![image-20211224214811233](Proteomics%20analysis.assets/image-20211224214811233.png)

~~~R
rownames(RNA_common) <- RNA_common$symbol
RNA_common1 <- RNA_common[,c(8:9,6:7)]
data_plot_zscore <- t(apply(RNA_common1, 1, function(x) (x-mean(x))/sd(x)))
data_plot_zscore <- na.omit(data_plot_zscore)
break1=seq(min(data_plot_zscore), 0, length.out=25)
break2=seq(0,max(data_plot_zscore), length.out=26)
the_break<-c(break1,break2[-1])
pheatmap(data_plot_zscore,show_rownames=TRUE,border_color=NA,cluster_cols=F,cluster_row=F,cellwidth=30,cellheight=14,fontsize=7,color = colorRampPalette(c("navy", "white","firebrick3"))(50),breaks=the_break)
~~~

![image-20211224214832934](Proteomics%20analysis.assets/image-20211224214832934.png)