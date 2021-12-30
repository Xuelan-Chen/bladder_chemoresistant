# Bulk RNA-seq analysis

## Part1. we should load the packages we need for following analysis

~~~R
require(Rsamtools)
require(GenomicFeatures)
require(GenomicAlignments)
require(BiocParallel)
require(pheatmap)
require(RColorBrewer)
require(PoiClaClu)
require(org.Mm.eg.db)
require(AnnotationDbi)
require(DOSE)
require(clusterProfiler)
require(topGO)
require(pathview)
require(org.Hs.eg.db)
require(AnnotationDbi)
require(DOSE)
require(clusterProfiler)
require(topGO)
require(ggplot2)
require("DESeq2")
require(future)
require(future.apply)
require(scales)
require(nichenetr)
require(tidyverse)
require(iTALK)
require(Seurat)
require(BuenColors)
options(future.globals.maxSize = 200 * 1024^3)
plan("multiprocess", workers = 20)
plan()
~~~

## Part2. we should load the data we need for following analysis

~~~R
count_and_tpm <- read.csv(file="./Bulk_RNAseq/csv_files/1_RNAseq_count_tpm_data.csv")
res_Sen_vs_normal_result <- read.csv(file="./Bulk_RNAseq/csv_files/1_RNAseq_Sen_vs_normal_count_tpm_all_result.csv")
res_Resis_vs_Sen_result <- read.csv(file="./Bulk_RNAseq/csv_files/1_RNAseq_Resis_vs_Sen_count_tpm_all_result.csv")
all_result_DDP <- read.csv(file="./Bulk_RNAseq/csv_files/1_RNAseq_DDP_vs_Sen_count_tpm_all_result.csv")
all_result_E64 <- read.csv(file="./Bulk_RNAseq/csv_files/2_RNAseq_E64_count_tpm_all_result.csv")
all_human_DDP_data <- read.csv(file="./Bulk_RNAseq/csv_files/3_RNAseq_human_DDP_count_tpm_data.csv")
~~~

## Part3. RNA-seq1 (sensitive tumor vs normal baldder) analysis

~~~R
load("./Bulk_RNAseq/RData_files/RNAseq_data1_se.RData")
load("./Bulk_RNAseq/RData_files/txdb_mm10.RData")
load("./Bulk_RNAseq/RData_files/ebg_mm10.RData")
sampleTable <- read.csv("./csv_files/RNAseq_sampletable1.csv")
countdata <- assay(se)
colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~ new)
exons <- as.data.frame(ebg)
exon_length_tmp <- future_lapply(1:nrow(assay(se)),function(i){
	exonsfromgene <- exons[which(exons$group==i),]
	length <- sum(exonsfromgene$width)
	gene_name <- exonsfromgene[1,2]
	exon_length <- data.frame(gene_name=gene_name,length=length)
	return(exon_length)
	})
exon_length <- as.data.frame(data.table::rbindlist(exon_length_tmp))
#求出每个样本的总read数（百万）
whole_length <- as.data.frame(round(colSums(assay(dds))/1e6,1))
#计算FPKM矩阵
countM <- as.data.frame(assay(se))
FPKM <- countM
FPKM_tmp <- future_lapply(1:nrow(FPKM),function(x){
	sel_data <- FPKM[x,]
	gene <- rownames(sel_data)
	length <- as.numeric(exon_length[which(exon_length$gene_name==gene),]$length)
	data_tmp <- future_lapply(1:ncol(sel_data),function(yy){
		tmp_data <- sel_data[,yy]/(whole_length[yy,]*length*0.001)
		return(tmp_data)
	})
	data_tmp <- as.data.frame(data_tmp)
	colnames(data_tmp) <- colnames(sel_data)
	data_tmp$id <- rownames(data_tmp)
	return(data_tmp)
	})
FPKM <- as.data.frame(data.table::rbindlist(FPKM_tmp))
rownames(FPKM) <- rownames(countM)
FPKM <- FPKM[,-ncol(FPKM)]
fpkmToTpm <- function(FPKM){
  exp(log(FPKM) - log(sum(FPKM)) + log(1e6))
}
tpm <- fpkmToTpm(FPKM)
count <- assay(dds)
colnames(tpm) <- paste0("tpm_",sampleTable$deal)
colnames(count) <- paste0("count_",sampleTable$deal)
count_and_tpm <- cbind(count,tpm)[,c(1:3,7:21,25:36)]
count_and_tpm$ENSEMBL <- mapIds(x = org.Mm.eg.db,
                      keys = rownames(count_and_tpm),
          keytype ="SYMBOL",
          column ="ENSEMBL",
          multiVals="first")
write.csv(count_and_tpm,file="./Bulk_RNAseq/csv_files/1_RNAseq_count_tpm_data.csv")

ddsMat <- DESeqDataSetFromMatrix(countData = countdata[,c(1:3,7:9)],
                                 colData = coldata[c(1:3,7:9),],
                                 design = ~ new)
ddsMat <- DESeq(ddsMat)
res_Sen_vs_normal <- results(ddsMat, contrast=c("new","mouse_Sensitive","normal_bladder"))
colnames(res_Sen_vs_normal) = paste("Sen_vs_normal",colnames(res_Sen_vs_normal),sep="_")
count_and_tpm <- read.csv(file="./Bulk_RNAseq/csv_files/1_RNAseq_count_tpm_data.csv")
res_Sen_vs_normal_result <- cbind(count_and_tpm[,c(1,32,2:7,17:22)],res_Sen_vs_normal)
res_Sen_vs_normal_result$Sen_vs_normal_threshold = as.factor(ifelse(res_Sen_vs_normal_result$Sen_vs_normal_pvalue < 0.05 & abs(res_Sen_vs_normal_result$Sen_vs_normal_log2FoldChange) > 1, ifelse(res_Sen_vs_normal_result$Sen_vs_normal_log2FoldChange < -1 ,'Down','Up'),'NoSignifi'))
res_Sen_vs_normal_result$entrez <- mapIds(x = org.Mm.eg.db,
                        keys = as.character(res_Sen_vs_normal_result$X),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
res_Sen_vs_normal_result <- as.data.frame(res_Sen_vs_normal_result)
res_Sen_vs_normal_result <- res_Sen_vs_normal_result[,c(3:21,1:2,22)]
colnames(res_Sen_vs_normal_result)[20] <- c("symbol")
write.csv(res_Sen_vs_normal_result,file="./Bulk_RNAseq/csv_files/1_RNAseq_Sen_vs_normal_count_tpm_all_result.csv")
~~~

## Part4. RNA-seq2 (resistant tumor vs sensitive tumor) analysis

~~~R
ddsMat <- DESeqDataSetFromMatrix(countData = countdata[,c(16:18,7:9)],
                                 colData = coldata[c(16:18,7:9),],
                                 design = ~ new)
ddsMat <- DESeq(ddsMat)
res_Resis_vs_Sen <- results(ddsMat, contrast=c("new","mouse_Resis","mouse_Sensitive"))
colnames(res_Resis_vs_Sen) = paste("Resis_vs_Sen",colnames(res_Resis_vs_Sen),sep="_")
count_and_tpm <- read.csv(file="./Bulk_RNAseq/csv_files/1_RNAseq_count_tpm_data.csv")
res_Resis_vs_Sen_result <- cbind(count_and_tpm[,c(1,32,5:7,14:16,20:22,29:31)],res_Resis_vs_Sen)
res_Resis_vs_Sen_result$Resis_vs_Sen_threshold = as.factor(ifelse(res_Resis_vs_Sen_result$Resis_vs_Sen_pvalue < 0.01 & abs(res_Resis_vs_Sen_result$Resis_vs_Sen_log2FoldChange) > 1, ifelse(res_Resis_vs_Sen_result$Resis_vs_Sen_log2FoldChange < -1 ,'Down','Up'),'NoSignifi'))
res_Resis_vs_Sen_result$entrez <- mapIds(x = org.Mm.eg.db,
                        keys = as.character(res_Resis_vs_Sen_result$X),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
res_Resis_vs_Sen_result <- as.data.frame(res_Resis_vs_Sen_result)
res_Resis_vs_Sen_result <- res_Resis_vs_Sen_result[,c(3:21,1:2,22)]
colnames(res_Resis_vs_Sen_result)[20] <- c("symbol")
write.csv(res_Resis_vs_Sen_result,file="./Bulk_RNAseq/csv_files/1_RNAseq_Resis_vs_Sen_count_tpm_all_result.csv")
~~~



## Part5. RNA-seq3 (DDP tumor vs sensitive tumor) analysis

~~~R
ddsMat <- DESeqDataSetFromMatrix(countData = countdata[,c(7:18)],
                                 colData = coldata[c(7:18),],
                                 design = ~ new)
ddsMat <- DESeq(ddsMat)
res_DDP2_vs_Sen <- results(ddsMat, contrast=c("new","mouse_DDP2","mouse_Sensitive"))
res_DDP4_vs_Sen <- results(ddsMat, contrast=c("new","mouse_DDP4","mouse_Sensitive"))
res_DDP4_vs_DDP2 <- results(ddsMat, contrast=c("new","mouse_DDP4","mouse_DDP2"))
colnames(res_DDP2_vs_Sen) = paste("DDP2_vs_Sen",colnames(res_DDP2_vs_Sen),sep="_")
colnames(res_DDP4_vs_Sen) = paste("DDP4_vs_Sen",colnames(res_DDP4_vs_Sen),sep="_")
colnames(res_DDP4_vs_DDP2) = paste("DDP4_vs_DDP2",colnames(res_DDP4_vs_DDP2),sep="_")
count_and_tpm <- read.csv(file="./Bulk_RNAseq/csv_files/1_RNAseq_count_tpm_data.csv")
all_result <- cbind(count_and_tpm[,c(1,32,5:16,20:31)],res_DDP2_vs_Sen)
all_result <- cbind(all_result,res_DDP4_vs_Sen)
all_result <- cbind(all_result,res_DDP4_vs_DDP2)

all_result$entrez <- mapIds(x = org.Mm.eg.db,
                        keys = as.character(all_result$X),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
all_result <- as.data.frame(all_result)
all_result <- all_result[,c(3:44,1:2,45)]
colnames(all_result)[43] <- c("symbol")
write.csv(all_result,file="./Bulk_RNAseq/csv_files/1_RNAseq_DDP_vs_Sen_count_tpm_all_result.csv")

squamouse_gene <- c("CD44","KRT6A","KRT14","TGM1","PI3","DSC3","GSDMC","TRP63","COL17A1","KRT5","KRT17","KRT1","CDH3","KRT16","KRT6B","KRT6C") 
squamouse_gene <- as.data.frame(squamouse_gene)
squamouse_gene$squamouse_gene <- as.character(squamouse_gene$squamouse_gene)
squamouse_gene1 = squamouse_gene %>% mutate(mouse_gene = convert_human_to_mouse_symbols(squamouse_gene)) %>% drop_na()

baldder_SCC_like_VS_other <- read.csv( file = "./csv_files/RNAseq_GSE32894_baldder_SCC_like_VS_other_res_pj_mouse.csv")
baldder_SCC_like_VS_other <- as.data.frame(baldder_SCC_like_VS_other)
baldder_SCC_like_VS_other <- baldder_SCC_like_VS_other[order(baldder_SCC_like_VS_other$t,decreasing=T),]

rownames(all_result) <- all_result$symbol
plot_data <- all_result[,c(13:24)]
plot_data <- log(plot_data+1,2)
bulk_data <- CreateSeuratObject(counts = plot_data, project = "MLL3")
bulk_data$cluster <- c("Sen","Sen","Sen","DDP2","DDP2","DDP2","DDP4","DDP4","DDP4","Res","Res","Res")

Lineage_marker <- head(baldder_SCC_like_VS_other$gene,100)
    tmp <- bulk_data
    Lineage_marker <- intersect(rownames(GetAssayData(object = tmp, slot = "data")),Lineage_marker)
    speci_raw <- FetchData(object = tmp, vars = Lineage_marker,slot="data")
bulk_data[["SCC_100"]] <- (rowSums(speci_raw))/length(Lineage_marker)

Lineage_marker <- squamouse_gene1$mouse_gene
tmp <- bulk_data
Lineage_marker <- intersect(rownames(GetAssayData(object = tmp, slot = "data")),Lineage_marker)
speci_raw <- FetchData(object = tmp, vars = Lineage_marker,slot="data")
bulk_data[["squamous_sig"]] <- (rowSums(speci_raw))/length(Lineage_marker)
bulk_data[["Ctsh"]] <- FetchData(object = bulk_data, vars = "Ctsh",slot="data")

my_comparisons <- list( c("Sen", "DDP2"),c("Sen","DDP4"),c("DDP2","DDP4"),c("DDP4","Res"))
ggviolin(bulk_data[[]], x = "cluster", y = "squamous_sig", fill = "cluster",
         palette = c("#00AFBB", "#E7B800", "#FC4E07", "#CA2927"),
         add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons =my_comparisons,
    label = "p.format", method = "t.test",label.y=c(1.7, 2.0, 2.3,2.6),paired=FALSE) 
~~~

![image-20211223175319307](Bulk%20RNA-seq%20analysis.assets/image-20211223175319307.png)

~~~R
ggviolin(bulk_data[[]], x = "cluster", y = "SCC_100", fill = "cluster",
         palette = c("#00AFBB", "#E7B800", "#FC4E07", "#CA2927"),
         add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons =my_comparisons,
    label = "p.format", method = "t.test",label.y=c(1.6, 1.9, 2.2,2.5),paired=FALSE) 
~~~

![image-20211223175403664](Bulk%20RNA-seq%20analysis.assets/image-20211223175403664.png)

~~~R
my_comparisons <- list( c("Sen", "DDP2"),c("Sen","DDP4"),c("Sen","Res"),c("DDP2","Res"))
ggbarplot(bulk_data[[]], x = "cluster", y = "Ctsh", fill="cluster",
       add = c("mean_se"),ylim=c(0,4), legend = "none",title="Normalised value")+
       rotate_x_text(angle = 45)+ stat_compare_means(comparisons =my_comparisons,
    label = "p.format", method = "t.test",label.y=c( 2.0, 2.5,3.0,3.5),paired=FALSE) +
    scale_fill_manual(values=as.character(jdb_palette("corona",4)))
~~~

![image-20211223175433737](Bulk%20RNA-seq%20analysis.assets/image-20211223175433737.png)

## Part6. RNA-seq4(resistant E64 vs resistant Vehicle) analysis

~~~R
setwd("/mnt/data/user_data/xuelan/project/3_singlecell/7_wangmanli/PDF_result/8_2020new/15_figure_share/Bulk_RNAseq")
load("./Bulk_RNAseq/RData_files/RNAseq_dataE64_se.RData")
load("./Bulk_RNAseq/RData_files/txdb_mm10.RData")
load("./Bulk_RNAseq/RData_files/ebg_mm10.RData")
sampleTable <- read.csv("./Bulk_RNAseq/csv_files/RNAseq_sampletable_64.csv")
countdata <- assay(se)
colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~ new)
exons <- as.data.frame(ebg)
exon_length_tmp <- future_lapply(1:nrow(assay(se)),function(i){
	exonsfromgene <- exons[which(exons$group==i),]
	length <- sum(exonsfromgene$width)
	gene_name <- exonsfromgene[1,2]
	exon_length <- data.frame(gene_name=gene_name,length=length)
	return(exon_length)
	})
exon_length <- as.data.frame(data.table::rbindlist(exon_length_tmp))
#求出每个样本的总read数（百万）
whole_length <- as.data.frame(round(colSums(assay(dds))/1e6,1))
#计算FPKM矩阵
countM <- as.data.frame(assay(se))
FPKM <- countM
FPKM_tmp <- future_lapply(1:nrow(FPKM),function(x){
	sel_data <- FPKM[x,]
	gene <- rownames(sel_data)
	length <- as.numeric(exon_length[which(exon_length$gene_name==gene),]$length)
	data_tmp <- future_lapply(1:ncol(sel_data),function(yy){
		tmp_data <- sel_data[,yy]/(whole_length[yy,]*length*0.001)
		return(tmp_data)
	})
	data_tmp <- as.data.frame(data_tmp)
	colnames(data_tmp) <- colnames(sel_data)
	data_tmp$id <- rownames(data_tmp)
	return(data_tmp)
	})
FPKM <- as.data.frame(data.table::rbindlist(FPKM_tmp))
rownames(FPKM) <- rownames(countM)
FPKM <- FPKM[,-ncol(FPKM)]
fpkmToTpm <- function(FPKM){
  exp(log(FPKM) - log(sum(FPKM)) + log(1e6))
}
tpm <- fpkmToTpm(FPKM)
count <- assay(dds)
colnames(tpm) <- paste0("tpm_",sampleTable$deal)
colnames(count) <- paste0("count_",sampleTable$deal)
count_and_tpm <- cbind(count[,c(2:5)],tpm[,c(2:5)])
count_and_tpm$ENSEMBL <- mapIds(x = org.Mm.eg.db,
                      keys = rownames(count_and_tpm),
          keytype ="SYMBOL",
          column ="ENSEMBL",
          multiVals="first")
write.csv(count_and_tpm,file="./Bulk_RNAseq/csv_files/2_RNAseq_E64_count_tpm_data.csv")

ddsMat <- DESeqDataSetFromMatrix(countData = countdata[,c(2:5)],
                                 colData = coldata[c(2:5),],
                                 design = ~ new)
ddsMat <- DESeq(ddsMat)
res_E64_vs_Vehicle <- results(ddsMat, contrast=c("new","Resis_E64","Resis_Vehicle"))
summary(res_E64_vs_Vehicle)
colnames(res_E64_vs_Vehicle) = paste("E64_vs_Vehicle",colnames(res_E64_vs_Vehicle),sep="_")
all_result <- cbind(count_and_tpm,res_E64_vs_Vehicle)
all_result$symbol <- rownames(all_result)
all_result <- as.data.frame(all_result)
all_result <- all_result[,c(1:8,10:16,9)]
all_result$E64_vs_Vehicle_threshold = as.factor(ifelse(all_result$E64_vs_Vehicle_pvalue < 0.05 & abs(all_result$E64_vs_Vehicle_log2FoldChange) > 1, ifelse(all_result$E64_vs_Vehicle_log2FoldChange < -1 ,'Down','Up'),'NoSignifi'))
all_result$entrez <- mapIds(x = org.Mm.eg.db,
                        keys = as.character(all_result$symbol),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
write.csv(all_result,file="./Bulk_RNAseq/csv_files/2_RNAseq_E64_count_tpm_all_result.csv")

all_E64_p_up <- subset(all_result,E64_vs_Vehicle_threshold =="Up")
all_E64_p_down <- subset(all_result,E64_vs_Vehicle_threshold =="Down")
dim(all_E64_p_up)
dim(all_E64_p_down)
ee  <-as.matrix(all_E64_p_up$entrez)
dd <- as.vector(ee)
GO_all_E64_p_up <- enrichGO(gene = dd, 
           OrgDb = org.Mm.eg.db,
      ont = "all", 
               pvalueCutoff = 0.01, 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05,
                   minGSSize = 10, 
                   maxGSSize = 500, 
                   readable = T, 
                   pool = FALSE)
GO_up_RNA <- as.data.frame(GO_all_E64_p_up)

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
  title="E64 Up GO BP")
~~~

![image-20211223211129279](Bulk%20RNA-seq%20analysis.assets/image-20211223211129279.png)

~~~R
GO_up_RNA_MF <- subset(GO_up_RNA,ONTOLOGY=="MF")
GO_up_RNA_MF$log10_Pvalue <- -log(GO_up_RNA_MF$pvalue,10)
GO_up_RNA_MF_5 <- head(GO_up_RNA_MF,5)
ggbarplot(GO_up_RNA_MF_5, 
  x = "Description", 
  y = "log10_Pvalue",
  color = "#FA6A63",            # Set bar border colors to white
  fill ="#FA6A63",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="E64 Up GO MF")
~~~

![image-20211223211157291](Bulk%20RNA-seq%20analysis.assets/image-20211223211157291.png)

## Part7. RNA-seq4 (human DDP) analysis

~~~R
load("./Bulk_RNAseq/RData_files/RNAseq_data_human_DDP_se.RData")
load("./Bulk_RNAseq/RData_files/txdb_hg19.RData")
load("./Bulk_RNAseq/RData_files/ebg_hg19.RData")
sampleTable <- read.csv("./Bulk_RNAseq/csv_files/RNAseq_sampletable_human_DDP.csv")
countdata <- assay(se)
colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~ new)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
colnames(normalized_counts) <- c("DESeq_BC20190927_T","DESeq_BC20190927_DDP12","DESeq_BC20190927_DDP19","DESeq_BC20200528_2_T","DESeq_BC20200528_2_DDP5","DESeq_BC20200528_T","DESeq_BC20200528_DDP6","DESeq_BC20200617_T","DESeq_BC20200617_DDP3")
normalized_counts <- as.data.frame(normalized_counts)
count_data <- as.data.frame(assay(se))
colnames(count_data) <- c("count_BC20190927_T","count_BC20190927_DDP12","count_BC20190927_DDP19","count_BC20200528_2_T","count_BC20200528_2_DDP5","count_BC20200528_T","count_BC20200528_DDP6","count_BC20200617_T","count_BC20200617_DDP3")
all_human_DDP_data <- cbind(count_data,normalized_counts)
all_human_DDP_data$SYMBOL <- rownames(all_human_DDP_data)
all_human_DDP_data$ENSEMBL <- mapIds(x = org.Hs.eg.db,
                      keys = rownames(all_human_DDP_data),
          keytype ="SYMBOL",
          column ="ENSEMBL",
          multiVals="first")
write.csv(all_human_DDP_data,file="./Bulk_RNAseq/csv_files/3_RNAseq_human_DDP_count_tpm_data.csv")
~~~



