# human sc-RNAseq data analysis

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
sc_human_merge <- mcreadRDS("./sc-RNAseq/RData_files/human_singlecell_RNAseq_all.rds",mc.cores=20)
sc_human_tumor_sdo <- mcreadRDS("./sc-RNAseq/RData_files/human_singlecell_RNAseq_tumor_pseudotime.rds",mc.cores=20)
sc_mouse_tumor_sdo <- mcreadRDS(file="./sc-RNAseq/RData_files/mouse_singlecell_RNAseq_tumor_pseudotime.rds",mc.cores=20)
~~~

## Part3. The detail codes of visualization of  human sc-RNAseq data

~~~R
Idents(object = sc_human_merge)  <- "xl_anno"
sc_human_merge$xl_anno <- factor(sc_human_merge$xl_anno,levels=c("Tumor1","Tumor2","Tcell","Macrophage","Stromal","endoderm","endothelia","Muscle","Bcell","mast"))
DimPlot(sc_human_merge, reduction = "tsne", group.by="xl_anno",label = FALSE, pt.size = 0.8,cols=as.character(jdb_palette("corona",10))) 
~~~

![image-20211225212319101](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212319101.png)

~~~R
sel_colors <- c("#2c80c5","#e10b2b")
names_clu <- c("Tumor1", "Tumor2")
col <- sel_colors
names(col) <- names_clu
sc_human_tumor_sdo$xl_anno <- factor(sc_human_tumor_sdo$xl_anno,levels=names_clu)
DimPlot(object = sc_human_tumor_sdo, reduction = "pca",
  label=FALSE,pt.size=1,
  group.by="xl_anno",
  cols=col[unique(sc_human_tumor_sdo$xl_anno)])
~~~

![image-20211225212442789](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212442789.png)

~~~R
sel_colors <- c("#2c80c5","#e10b2b")
names_clu <- c("human_Sen", "human_Ren")
col <- sel_colors
names(col) <- names_clu
sc_human_tumor_sdo$group <- factor(sc_human_tumor_sdo$group,levels=names_clu)
DimPlot(object = sc_human_tumor_sdo, reduction = "pca",
  label=FALSE,pt.size=1,
  group.by="group",
  cols=col[unique(sc_human_tumor_sdo$group)])
~~~

![image-20211225212504189](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212504189.png)

~~~R
table(sc_human_tumor_sdo$xl_anno)
table(sc_human_tumor_sdo$group)
sc_human_tumor_sdo$new_type <- paste(sc_human_tumor_sdo$group,sc_human_tumor_sdo$xl_anno,sep="_")
type <- c("Tumor2","Tumor1","Tumor2","Tumor1")
origin <- c("Sen","Sen","Ren","Ren")
percentage <- c(28,177,3084,434)
all_data <- data.frame(type,origin,percentage)
all_data$origin <- factor(all_data$origin,levels=c("Sen","Ren"))
all_data$type <- factor(all_data$type,levels=c("Tumor1","Tumor2"))
ggplot(all_data, aes(origin,percentage, fill=type))+geom_bar(stat='identity',position='fill') +labs(x = 'origin',y = 'tumor percentage') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'))
~~~

![image-20211225212527694](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212527694.png)

~~~R
XY_FeaturePlot(object = sc_human_tumor_sdo, features = c("CTSH","EPCAM"),pt.size=.6,reduction="pca",label=T,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000")) + NoLegend()
~~~

![image-20211225212552821](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212552821.png)

~~~R
E_64_targets <- c("CTSK","CTSL","CTSS","CTSB","CTSH","CAPN1","CAPN2","CAPN3","CAPN5","CAPN6","CAPN7","CAPN8","CAPN9","CAPN10","CAPN11","CAPN12","CAPN13","CAPN14","CAPN15")
E_64_targets <- as.data.frame(E_64_targets)
Lineage_marker <- unique(E_64_targets$E_64_targets)
 tmp <- sc_human_tumor_sdo
 Lineage_marker <- intersect(rownames(GetAssayData(object = tmp, slot = "data")),Lineage_marker)
 speci_raw <- FetchData(object = tmp, vars = Lineage_marker,slot="data")
sc_human_tumor_sdo[["E64_targets"]] <- (rowSums(speci_raw))/length(Lineage_marker)
XY_FeaturePlot(object = sc_human_tumor_sdo, features = c("E64_targets") ,pt.size=1,reduction="pca",label=T,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000"))+ NoLegend()
~~~

![image-20211225212612788](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212612788.png)

~~~R
ALL_GSEA_GMT <- read.gmt("./msigdb.v7.1.symbols.gmt")
ALL_GSEA_GMT$ont <- as.character(ALL_GSEA_GMT$ont)
used_pathway <- c("GO_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION","GO_EPITHELIAL_CELL_DIFFERENTIATION","GO_KERATINIZATION","GO_KERATIN_FILAMENT","GO_KERATINOCYTE_MIGRATION","GO_REGULATION_OF_KERATINOCYTE_MIGRATION","REACTOME_KERATINIZATION","GO_KERATINOCYTE_DIFFERENTIATION")
KERATIN_pathway  <- XY_subset(ALL_GSEA_GMT,"ont",used_pathway)
Lineage_marker <- unique(KERATIN_pathway$gene)
 tmp <- sc_human_tumor_sdo
 Lineage_marker <- intersect(rownames(GetAssayData(object = tmp, slot = "data")),Lineage_marker)
 speci_raw <- FetchData(object = tmp, vars = Lineage_marker,slot="data")
sc_human_tumor_sdo[["KERATIN_pathway"]] <- (rowSums(speci_raw))/length(Lineage_marker)
XY_FeaturePlot(object = sc_human_tumor_sdo, features = c("KERATIN_pathway") ,pt.size=1,reduction="pca",label=T,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000"))+ NoLegend()
~~~

![image-20211225212819830](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212819830.png)

~~~R
resProgressive_vs_Response <- read.csv(file="./sc-RNAseq/RData_files/TCGA_BLCA_resis_vs_respon_DSEq2_result.csv")
resProgressive_vs_Response <- as.data.frame(resProgressive_vs_Response)
resProgressive_vs_Response_p <- subset(resProgressive_vs_Response,pvalue < 0.05)
resProgressive_vs_Response_p$X <- as.character(resProgressive_vs_Response_p$X)
TCGA_resis = resProgressive_vs_Response_p %>% mutate(mouse_gene = convert_human_to_mouse_symbols(X)) %>% drop_na()
TCGA_resis <- TCGA_resis[order(TCGA_resis$log2FoldChange,decreasing=TRUE),]
tmp_pathway <- head(TCGA_resis,500)
Lineage_marker <- tmp_pathway$X
tmp <- sc_human_tumor_sdo
Lineage_marker <- intersect(rownames(GetAssayData(object = tmp, slot = "data")),Lineage_marker)
speci_raw <- FetchData(object = tmp, vars = Lineage_marker,slot="data")
sc_human_tumor_sdo[[c("TCGA_resis_Up")]] <- (rowSums(speci_raw))/length(Lineage_marker)
XY_FeaturePlot(object = sc_human_tumor_sdo, features = c("TCGA_resis_Up") ,pt.size=1,reduction="pca",label=T,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000"))+ NoLegend()
~~~



![image-20211225212853292](human%20sc-RNAseq%20data%20analysis.assets/image-20211225212853292.png)

~~~R
Idents(object = sc_mouse_tumor_sdo)  <- "xl_anno"
sc_mouse_tumor_sdo$xl_anno <- factor(sc_mouse_tumor_sdo$xl_anno,levels=c("Tumor1","Tumor2"))
DimPlot(object = sc_mouse_tumor_sdo, reduction = "dm")

All_counts <- as.data.frame(GetAssayData(object = sc_mouse_tumor_sdo, slot = "counts",assay="RNA"))
All_counts1 <- All_counts
All_counts1$mmu_gene <- rownames(All_counts1)
All_counts1 <- All_counts1 %>% mutate(hsa_gene = convert_mouse_to_human_symbols(mmu_gene)) %>% drop_na()
All_counts1 <- All_counts1[!duplicated(All_counts1$hsa_gene),]
rownames(All_counts1) <- All_counts1$hsa_gene
All_counts1 <- All_counts1[,colnames(All_counts)]
sc_mouse_tumor_sdo_hsa <- CreateSeuratObject(counts = All_counts1,assay = 'RNA',project = 'RNA',min.cells = 0,meta.data = sc_mouse_tumor_sdo@meta.data[colnames(All_counts1),])
sc_mouse_tumor_sdo_hsa <- sc_mouse_tumor_sdo_hsa %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
    ScaleData(verbose = TRUE) %>% 
    RunPCA(pc.genes = sc_mouse_tumor_sdo_hsa@var.genes, npcs = 30, verbose = FALSE)

DimPlot(object = sc_human_tumor_sdo, reduction = "pca",pt.size=0.3,group.by="xl_anno", label = TRUE, repel = TRUE)
DefaultAssay(sc_human_tumor_sdo) <- "RNA"
sc_human_tumor_sdo <- sc_human_tumor_sdo %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
    ScaleData(verbose = TRUE) %>% 
    RunPCA(pc.genes = sc_human_tumor_sdo@var.genes, npcs = 30, verbose = FALSE)
sc_human_tumor_sdo$xl_anno <- factor(sc_human_tumor_sdo$xl_anno,levels=c("Tumor1","Tumor2"))
sc_human_tumor_sdo.anchors <- FindTransferAnchors(reference = sc_mouse_tumor_sdo_hsa, query = sc_human_tumor_sdo, dims = 1:30)
predictions <- TransferData(anchorset = sc_human_tumor_sdo.anchors, refdata = sc_mouse_tumor_sdo_hsa$xl_anno, dims = 1:30,k.weight=5)
sc_human_tumor_sdo$predicted.id <- predictions$predicted.id
sc_human_tumor_sdo$prediction.score.max <- predictions$prediction.score.max
DimPlot(object = sc_human_tumor_sdo, reduction = "pca",label=T,group.by="predicted.id")

b = table(sc_human_tumor_sdo$predicted.id,sc_human_tumor_sdo$xl_anno)
colnames(b) <- c("H.Tumor1","H.Tumor2")
rownames(b) <- c("M.Tumor1","M.Tumor2")
break1=seq(0, 0.5, length.out=50)
break2=seq(0.5, 1, length.out=51)
the_break<-c(break1,break2[-1])
pheatmap(b/rowSums(b),border_color=NA,cluster_cols=F,breaks=the_break,
  color = colorRampPalette(c("navy", "white","firebrick3"))(100),show_rownames=TRUE,cluster_rows=F,main="Cor between models and Patient")
~~~

![image-20211225213008860](human%20sc-RNAseq%20data%20analysis.assets/image-20211225213008860.png)

