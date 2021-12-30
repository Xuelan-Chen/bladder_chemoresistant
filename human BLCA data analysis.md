# human BLCA data analysis

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
cornell_2016_sample <- read.table("./human_BLCA/BLCA_data1/1_blca_cornell_2016_clinical_sample.txt",fill = TRUE, header = TRUE,sep="\t")
cornell_2016_CNA <- read.table("./human_BLCA/BLCA_data1/1_blca_cornell_2016_CNA.txt",fill = TRUE, header = TRUE,sep="\t")
TCGA_BLCA1 <- mcreadRDS("./human_BLCA/BLCA_data2/TCGA_BLCA_normalised.rds",mc.cores=20)
TCGA_BLCA <- mcreadRDS("./human_BLCA/BLCA_data2/TCGA_BLCA_raw.rds",mc.cores=20)
TCGA_BLCA_clinical <- mcreadRDS("./human_BLCA/BLCA_data2/TCGA_BLCA_clinical.rds",mc.cores=20)
sdo_timegenes <- mcreadRDS(file="./sc-RNAseq/RData_files/1_mouse_sc_RNA_pseudotime_gene_module_gene.rds",mc.cores=20)
~~~

## Part3. the detail codes of visualization of  BLCA data

~~~R
common_sample <- intersect(colnames(cornell_2016_CNA)[3:55],cornell_2016_sample$SAMPLE_ID)
cornell_2016_sample_used <- XY_subset(cornell_2016_sample,"SAMPLE_ID",common_sample)
cornell_2016_CNA_used  <- cbind(cornell_2016_CNA[,c(1,2)],cornell_2016_CNA[,common_sample])
rownames(cornell_2016_CNA_used) <- cornell_2016_CNA_used$Hugo_Symbol
cornell_2016_CNA_used <- cornell_2016_CNA_used[,-c(1,2)]
XL_CNV <- cornell_2016_CNA_used[c("CTSK","CTSL","CTSS","CTSB","CTSH","CAPN1","CAPN2","CAPN3","CAPN5","CAPN7","CAPN8","CAPN9","CAPN10","CAPN11","CAPN12","CAPN13","CAPN14","CAPN15"),]

group_prechemo <- subset(cornell_2016_sample_used,SPECIMEN_COLLECTION_PRE_OR_POST_CHEMO=="pre-chemotherapy")
group_postchemo <- subset(cornell_2016_sample_used,SPECIMEN_COLLECTION_PRE_OR_POST_CHEMO=="post-chemotherapy")
XL_CNV_prechemo <- XL_CNV[,as.character(group_prechemo$SAMPLE_ID)]
XL_CNV_postchemo <- XL_CNV[,as.character(group_postchemo$SAMPLE_ID)]
XL_CNV_prechemo$per_number_Amp <- apply(XL_CNV_prechemo,1,function (x) length(which(x==1 |x==2 )))
XL_CNV_prechemo$per_number_del <- apply(XL_CNV_prechemo,1,function (x) length(which(x==-1|x==-2 )))
XL_CNV_prechemo$per_Amp <- ((XL_CNV_prechemo$per_number_Amp)/length(group_prechemo$SAMPLE_ID))*100
XL_CNV_prechemo$per_del <- ((XL_CNV_prechemo$per_number_del)/length(group_prechemo$SAMPLE_ID))*100
XL_CNV_postchemo$post_number_Amp <- apply(XL_CNV_postchemo,1,function (x) length(which(x==1 |x==2 )))
XL_CNV_postchemo$post_number_del <- apply(XL_CNV_postchemo,1,function (x) length(which(x==-1|x==-2 )))
XL_CNV_postchemo$post_Amp <- ((XL_CNV_postchemo$post_number_Amp)/length(group_postchemo$SAMPLE_ID))*100
XL_CNV_postchemo$post_del <- ((XL_CNV_postchemo$post_number_del)/length(group_postchemo$SAMPLE_ID))*100
XL_CNV_data_used <- cbind(XL_CNV_prechemo[,c("per_number_Amp","per_number_del","per_Amp","per_del")],XL_CNV_postchemo[,c("post_number_Amp","post_number_del","post_Amp","post_del")])
XL_CNV_data_used$post_pre_Amp <- XL_CNV_data_used$post_Amp - XL_CNV_data_used$per_Amp
XL_CNV_data_used$post_pre_del <- XL_CNV_data_used$post_del - XL_CNV_data_used$per_del
XL_CNV_data_used <- XL_CNV_data_used[order(XL_CNV_data_used$post_pre_Amp,decreasing=TRUE),]

XL_CNV_data_used_new <- XL_CNV_data_used[,c("per_Amp","post_Amp")]
XL_CNV_data_used_new$gene <- rownames(XL_CNV_data_used_new)
data1 <- XL_CNV_data_used_new[,c("per_Amp","gene")]
data1$group <- c("per_chem")
data2 <- XL_CNV_data_used_new[,c("post_Amp","gene")]
data2$group <- c("post_chem")
colnames(data1)[1] <- c("Amp_ratio")
colnames(data2)[1] <- c("Amp_ratio")
all_data_plot <- rbind(data1,data2)
all_data_plot$group <- factor(all_data_plot$group,levels=c("per_chem","post_chem"))

CTSH_data <- subset(all_data_plot,gene=="CTSH")
chidata<- matrix(c(3,16,13,21),ncol=2,byrow=F)
crowname <- c("Amp","no_Amp")
ccolname <- c("Pre_chem","post_chem")
dimnames(chidata) <- list(group=ccolname,result=crowname)
chidata
chisq.test(chidata,correct=F)

ggbarplot(CTSH_data, x = "group", y = "Amp_ratio",color="group",
 position = position_dodge())+ labs(title="CTSH per_Amp_n=3 post_Amp_n=16 per_n=16 post_n=37 p=0.08783")
~~~

![image-20211225185908825](J:%5C%E9%99%88%E9%9B%AA%E5%85%B0%5C%E8%84%9A%E6%9C%AC%E8%AE%B0%E5%BD%95_WIN%5C%E8%84%9A%E6%9C%AC%E8%AE%B0%E5%BD%95_WIN%5Cproject%5C14_bladder_wangmanli%5C2_bladder_resistance%5CMIBC_Bulk_RNAseq%5Chuman_BLCA%20analysis.assets%5Cimage-20211225185908825.png)

~~~R
ggpaired(all_data_plot, x = "group", y = "Amp_ratio",
         color = "group", line.color = "gray", line.size = 0.4,
         palette = "jco")+
stat_compare_means(paired = TRUE)+ labs(title="E64 cluster gene ")
~~~

![image-20211225190012422](J:%5C%E9%99%88%E9%9B%AA%E5%85%B0%5C%E8%84%9A%E6%9C%AC%E8%AE%B0%E5%BD%95_WIN%5C%E8%84%9A%E6%9C%AC%E8%AE%B0%E5%BD%95_WIN%5Cproject%5C14_bladder_wangmanli%5C2_bladder_resistance%5CMIBC_Bulk_RNAseq%5Chuman_BLCA%20analysis.assets%5Cimage-20211225190012422.png)

~~~R
sdo_timegenes$symbol <- rownames(sdo_timegenes)
sdo_timegenes$symbol <- as.character(sdo_timegenes$symbol)
sdo_timegenes$cluster <- as.character(sdo_timegenes$Module)
table(sdo_timegenes$cluster)
TCGA_BLCA_Sel_ <- future_lapply(1:length(unique(sdo_timegenes$cluster)),function(x){
  tmp_clu <- subset(sdo_timegenes,cluster==unique(sdo_timegenes$cluster)[x])
  tmp_clu = tmp_clu %>% mutate(from = convert_mouse_to_human_symbols(symbol), to = convert_mouse_to_human_symbols(symbol)) %>% drop_na()
  blca_tcga_pub_tmp <- data.frame(tmp=as.character(apply(TCGA_BLCA1[,intersect(colnames(TCGA_BLCA1),unique(tmp_clu$from))],1,mean)),
    tmp2=as.character(apply(TCGA_BLCA1[,intersect(colnames(TCGA_BLCA1),unique(tmp_clu$from))],1,mean)),
    row.names=rownames(TCGA_BLCA1))
  blca_tcga_pub_tmp$tmp <- as.numeric(as.character(scale((as.numeric(as.character(blca_tcga_pub_tmp$tmp))))))
# blca_tcga_pub_tmp$tmp <- as.numeric(as.character(blca_tcga_pub_tmp$tmp))
  colnames(blca_tcga_pub_tmp) <- c(unique(sdo_timegenes$cluster)[x],"tmp2")
  return(blca_tcga_pub_tmp[,1])
  })
TCGA_BLCA_Sel <- do.call(cbind,TCGA_BLCA_Sel_)
TCGA_BLCA_Sel <- as.data.frame(TCGA_BLCA_Sel)
colnames(TCGA_BLCA_Sel) <- unique(sdo_timegenes$cluster)
rownames(TCGA_BLCA_Sel) <- rownames(TCGA_BLCA)

TCGA_BLCA_clinical_sel <- TCGA_BLCA_clinical[rownames(TCGA_BLCA_Sel),]
Escc_clinical_ <- cbind(TCGA_BLCA_clinical_sel,TCGA_BLCA_Sel)
meta <- Escc_clinical_
meta$days_to_last_follow_up[is.na(meta$days_to_last_follow_up)] <- "HHH"
tmp <- subset(meta,days_to_last_follow_up=="HHH")
tmp$days_to_last_follow_up <- tmp$days_to_death
no_na <- meta[setdiff(rownames(meta),rownames(tmp)),]
all_merge <- rbind(tmp,no_na)
all_merge <- subset(all_merge,days_to_last_follow_up != "HHH")
all_merge$vital_status <- as.character(all_merge$vital_status)
all_merge$status <- ifelse(all_merge$vital_status=="Alive",0,1)
all_merge$days_to_last_follow_up <- as.numeric(all_merge$days_to_last_follow_up)

all_merge.cut <- surv_cutpoint(
   all_merge,
   time = "days_to_last_follow_up",
   event = "status",
   variables = c("Module_4"),
   progressbar=TRUE,
   minprop=0.1
)
summary(all_merge.cut)
plot(all_merge.cut, "Module_4")
all_merge.cut.cat <- surv_categorize(all_merge.cut) 
library(survival)
fit <- survfit(Surv(days_to_last_follow_up, status) ~ Module_4, data = all_merge.cut.cat)
ggsurvplot(fit, data = all_merge.cut.cat,
surv.median.line = "hv",
pval = TRUE,
ggtheme = theme_bw(),
risk.table=TRUE)
~~~

![image-20211225190058999](C:%5CUsers%5CAdministrator%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20211225190058999.png)