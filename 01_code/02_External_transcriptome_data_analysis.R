## =========================================
## 1.Differential expression analysis
if(F){
  ## for bulk sample
  library(data.table)
  library(tximport)
  library(DESeq2)
  dir<-dir("/media/data1/YJ/TNBC_LX/2_kallisto/3_T_pt_cdna/",patter="_T")
  files<-file.path("/media/data1/YJ/TNBC_LX/2_kallisto/3_T_pt_cdna/",dir,"abundance.h5")
  names(files)<-dir
  
  dir<-dir("/media/data1/YJ/TNBC_LX/2_kallisto/4_N_pt_cdna/",patter="_N")
  tmp<-file.path("/media/data1/YJ/TNBC_LX/2_kallisto/4_N_pt_cdna/",dir,"abundance.h5")
  names(tmp)<-dir
  files<-c(files,tmp)
  
  kallisto.allsample<-tximport(files,type="kallisto",tx2gene = tx2gene,ignoreTxVersion = T)
 
  ### Perform differential expression analysis to find DEGs
  sampleTable <- data.frame(condition = gsub("P.*_","",colnames(kallisto.allsample$counts)))
  sampleTable$condition<-factor(sampleTable$condition,levels = c("N","T"))
  rownames(sampleTable) <- colnames(kallisto.allsample$counts)
  head(sampleTable)
  
  dds <- DESeqDataSetFromTximport(kallisto.allsample, sampleTable, ~condition)
  keep2<-rowSums(counts(dds)>0)>=3
  dds<-dds[keep2,]
  dds<-DESeq(dds)
  bulk_DEG.res<-results(dds)
  head(bulk_DEG.res)
  res.sig<-as.data.frame(bulk_DEG.res)
  res.sig<-res.sig[res.sig$log2FoldChange > 0.5 & res.sig$padj <0.01,]
  res.sig<-res.sig[order(res.sig$log2FoldChange,decreasing = T),]
  res.sig<-res.sig[!is.na(res.sig$baseMean),]
  dim(res.sig[!is.na(res.sig$baseMean),])
  
  ## for scRNA-seq data
  library(Seurat)
  library(Matrix)
  library(scibet)
  dir<-dir("./data/scRNA/",patter="GSM")
  files<-file.path("./data/scRNA",dir)
  names(files)<-dir
  TNBC_scRNA.list<-list()
  for(sample in names(files)){
    data <- Read10X(data.dir = files[sample])
    # Initialize the Seurat object with the raw (non-normalized data).
    data <- CreateSeuratObject(counts = data, project = sample, min.cells = 3, min.features = 500)
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    data <- subset(data, subset = percent.mt < 20)
    data$sample<-sample
    data$sample<-gsub(".*_","",data$sample)
    TNBC_scRNA.list[[sample]]<-data
  }
  TNBC_scRNA.seurat<-merge(TNBC_scRNA.list[["GSM4909255_N-N280-Epi"]], y = list(TNBC_scRNA.list[["GSM4909256_N-PM0095-Epi"]],
                                                                                TNBC_scRNA.list[["GSM4909258_N-NF-Epi"]],
                                                                                TNBC_scRNA.list[["GSM4909259_N-NE-Epi"]],
                                                                                TNBC_scRNA.list[["GSM4909281_TN-MH0126"]],
                                                                                TNBC_scRNA.list[["GSM4909282_TN-MH0135"]],
                                                                                TNBC_scRNA.list[["GSM4909283_TN-SH0106"]],
                                                                                TNBC_scRNA.list[["GSM4909284_TN-MH0114-T2"]]))
  TNBC_scRNA.seurat <- NormalizeData(TNBC_scRNA.seurat)
  TNBC_scRNA.seurat <- FindVariableFeatures(TNBC_scRNA.seurat, selection.method = "vst", nfeatures = 2000)
  TNBC_scRNA.seurat <- ScaleData(TNBC_scRNA.seurat, vars.to.regress = "percent.mt")
  TNBC_scRNA.seurat <- RunPCA(TNBC_scRNA.seurat, features = VariableFeatures(object = TNBC_scRNA.seurat))
  TNBC_scRNA.seurat <- FindNeighbors(TNBC_scRNA.seurat, dims = 1:20)
  TNBC_scRNA.seurat <- FindClusters(TNBC_scRNA.seurat, resolution = c(0.5,1))
  TNBC_scRNA.seurat <- RunUMAP(TNBC_scRNA.seurat, dims = 1:20)
  DimPlot(TNBC_scRNA.seurat, reduction = "umap",group.by = "sample")
  FeaturePlot(TNBC_scRNA.seurat, features = c("EPCAM"))
  VlnPlot(TNBC_scRNA.seurat, features = c("EPCAM"))
  ### Perform scibet analysis to determine cell type
  model<-readr::read_csv("/media/data1/YJ/HCC_singlecell/3_Ranalysis/major_human_cell_types.csv")
  model<-pro.core(model)
  query<-t(TNBC_scRNA.seurat@assays$RNA@data)
  prd <- LoadModel(model)
  label <- prd(query)
  label<-data.frame(scibat=label)
  rownames(label)<-rownames(TNBC_scRNA.seurat@meta.data)
  label$CellID<-rownames(label)
  
  TNBC_scRNA.seurat$scibet_label<-label$scibat
  
  scibet_res<-as.data.frame(rbind(table(TNBC_scRNA.seurat$RNA_snn_res.1,TNBC_scRNA.seurat$scibet_label)))
  
  ### Epithelial cells were selected according to scibet and the expression of EPCAM
  TNBC_scRNA.seurat<-RenameIdents(TNBC_scRNA.seurat, 
                                  "0"="epithelial",
                                  "1"="epithelial",
                                  "2"="epithelial",
                                  "3"="epithelial",
                                  "4"="epithelial",
                                  "5"="epithelial",
                                  "6"="epithelial",
                                  "8"="epithelial",
                                  "9"="epithelial",
                                  "11"="epithelial",
                                  "13"="epithelial",
                                  "14"="epithelial",
                                  "15"="epithelial",
                                  "16"="epithelial",
                                  "19"="epithelial",
                                  "22"="epithelial",
                                  "23"="epithelial",
                                  "24"="epithelial",
                                  "25"="epithelial",
                                  
                                  "7"="other",
                                  "10"="other",
                                  "12"="other",
                                  "17"="other",
                                  "18"="other",
                                  "20"="other",
                                  "21"="other"
  )
  TNBC_scRNA.seurat$celltype<-Idents(TNBC_scRNA.seurat)
  table(TNBC_scRNA.seurat$celltype,TNBC_scRNA.seurat$sample)
  TNBC_scRNA.seurat$tissue<-gsub("-.*","",TNBC_scRNA.seurat$sample)
  head(TNBC_scRNA.seurat@meta.data)
  table(TNBC_scRNA.seurat$tissue)
  
  TNBC_scRNA_epi.seurat<-subset(TNBC_scRNA.seurat,subset=celltype=="epithelial")
  ## Screening of differentially expressed genes between normal and tumor epithelial cells
  Idents(TNBC_scRNA_epi.seurat)<-"tissue"
  Tepi.de.markers <- FindMarkers(TNBC_scRNA_epi.seurat, ident.1 = "TN", ident.2 = "N",logfc.threshold = 0)
  Tepi.de.markers["PSMG2",]
  Tepi.de.markers.sig<-Tepi.de.markers[Tepi.de.markers$p_val_adj<0.01 & Tepi.de.markers$avg_logFC>0.2,]
  dim(Tepi.de.markers.sig)
  head(Tepi.de.markers.sig)
}

## ======================================================
## 2. Pathway Correlation Analysis
## For scRNA-seq data GSE75688
if(F){
  sc_TPM<-read.table("./data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",header = T) 
  sc_sampleinfo<-read.table("./data/GSE75688_final_sample_information.txt",header = T)  
  sc_sampleinfo$Pt<-gsub("_.*","",sc_sampleinfo$sample)
  ## Screen out the tumor cells
  sc_sampleinfo<-sc_sampleinfo[(sc_sampleinfo$Pt %in% c("BC07","BC08","BC09","BC10","BC11")) & sc_sampleinfo$type=="SC",]
  sc_sampleinfo<-sc_sampleinfo[sc_sampleinfo$index=="Tumor",]
  head(sc_sampleinfo)
  
  sc_TPM<-sc_TPM[,c("gene_name",sc_sampleinfo$sample)]
  sc_TPM$mean<-rowMeans(sc_TPM[,2:ncol(sc_TPM)])
  sc_TPM<-sc_TPM[order(sc_TPM$gene_name,-sc_TPM$mean),]
  sc_TPM<-sc_TPM[!duplicated(sc_TPM$gene_name),]
  rownames(sc_TPM)<-sc_TPM$gene_name
  sc_TPM<-sc_TPM[,2:(ncol(sc_TPM)-1)]
  sc_TPM[1:3,1:6]
  
  ## Perform GSVA analysis
  data<-as.data.frame(sc_TPM)
  data<-data[rowSums(data)>0,]
  data<-log2(data+1)
  h<-getGmt("/media/data1/YJ1/TNBC_LX_Nano/R_analysis_nano/data/h.all.v7.2.symbols.gmt")
  ## Screening out the proteasome pathway
  c5<-getGmt("/media/data1/YJ1/TNBC_LX_Nano/R_analysis_nano/data/c5.go.v7.2.symbols.gmt")
  c2<-getGmt("/media/data1/YJ1/TNBC_LX_Nano/R_analysis_nano/data/c2.all.v7.2.symbols.gmt")
  c2_c5_geneset<-GeneSetCollection(c(c2,c5))
  tmp<-c2_c5_geneset[c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY")]
  h<-GeneSetCollection(c(h,tmp))
  sc_H_gsva<-gsva(as.matrix(data),h,method="gsva",kcdf="Gaussian",parallel.sz=32)
  
  data_cor <- rcorr(as.matrix(t(sc_H_gsva)))
  CorMatrix <- function(cor,p) {
    ut <- upper.tri(cor) 
    data.frame(row = rownames(cor)[row(cor)[ut]] ,
               column = rownames(cor)[col(cor)[ut]], 
               cor =(cor)[ut], 
               p = p[ut] )
  }
  data_cor<-CorMatrix(data_cor$r,data_cor$P)
  data_cor<-data_cor[data_cor$column %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  data_cor<-data_cor[!data_cor$row %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  GSE75688_data_cor<-data_cor[data_cor$column=="KEGG_PROTEASOME",]
  GSE75688_data_cor<-GSE75688_data_cor[order(GSE75688_data_cor$row),]
}

## For scRNA-seq data GSE11838
if(F){
  GSE11838_TPM<-read.table("./data/GSE118389_tpm_rsem.txt",header = T)
  GSE11838_TPM[1:3,1:3]
  GSE11838_info<-read.table("./data/cell_types_tab_S9.txt")
  GSE11838_info<-GSE11838_info[GSE11838_info$V2=="epithelial",]
  head(GSE11838_info)
  dim(GSE11838_info)
  ## screening out tumor cells
  GSE11838_TPM<-GSE11838_TPM[,GSE11838_info$V1]
  ## Perform GSVA analysis
  data<-as.data.frame(GSE11838_TPM)
  data<-data[rowSums(data)>0,]
  data<-log2(data+1)
  
  GSE11838_H_gsva<-gsva(as.matrix(data),h,method="gsva",kcdf="Gaussian",parallel.sz=32)
  
  data_cor <- rcorr(as.matrix(t(GSE11838_H_gsva)))
  CorMatrix <- function(cor,p) {
    ut <- upper.tri(cor) 
    data.frame(row = rownames(cor)[row(cor)[ut]] ,
               column = rownames(cor)[col(cor)[ut]], 
               cor =(cor)[ut], 
               p = p[ut] )
  }
  data_cor<-CorMatrix(data_cor$r,data_cor$P)
  data_cor<-data_cor[data_cor$column %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  data_cor<-data_cor[!data_cor$row %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  GSE11838_data_cor<-data_cor[data_cor$column=="KEGG_PROTEASOME",]
  GSE11838_data_cor<-GSE11838_data_cor[order(GSE11838_data_cor$row),]
  save(GSE11838_data_cor,GSE75688_data_cor,file="./data/Hallmark_gene_set_cor.Rdata")
}

## Bulk sequence data PRJNA553096
if(F){
  tx2gene<-read.table("/data/reference/GRCh38.97/tx2gene2symbol.txt")
  tx2gene<-tx2gene[,c(2,3,1)]
  
  library(data.table)
  library(tximport)
  dir<-dir("/media/data1/YJ/TNBC_LX/2_kallisto/3_T_pt_cdna/",patter="_T")
  files<-file.path("/media/data1/YJ/TNBC_LX/2_kallisto/3_T_pt_cdna/",dir,"abundance.h5")
  names(files)<-dir
  kallisto<-tximport(files,type="kallisto",tx2gene = tx2gene,ignoreTxVersion = T)
  own_TNBC_tmp<-as.data.frame(kallisto$abundance)
  
  ## Perform GSVA analysis
  data<-as.data.frame(own_TNBC_tmp)
  data<-data[rowSums(data)>0,]
  data<-log2(data+1)

  own_H_gsva<-gsva(as.matrix(data),h,method="gsva",kcdf="Gaussian",parallel.sz=32)
  data_cor <- rcorr(as.matrix(t(own_H_gsva)))
  CorMatrix <- function(cor,p) {
    ut <- upper.tri(cor) 
    data.frame(row = rownames(cor)[row(cor)[ut]] ,
               column = rownames(cor)[col(cor)[ut]], 
               cor =(cor)[ut], 
               p = p[ut] )
  }
  data_cor<-CorMatrix(data_cor$r,data_cor$P)
  data_cor<-data_cor[data_cor$column %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  data_cor<-data_cor[!data_cor$row %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  own_data_cor<-data_cor[data_cor$column=="KEGG_PROTEASOME",]
  own_data_cor<-own_data_cor[order(own_data_cor$row),]
}

## Bulk sequence data TCGA-TNBC
if(F){
  ### Screening TNBC samples in TCGA-BRCA data
  B_cli<-read.table("/media/data1/YJ/TCGA/TCGA_BRCA_clinical.txt",header = T,sep="\t",quote = "")
  B_cli[1:3,1:3]
  B_cli$A0_Samples<-gsub("-",".",B_cli$A0_Samples)
  molecule<-B_cli[,c(1,19,20,42)]
  head(molecule)
  colnames(molecule)<-c("ID","ER","PR","her2")
  TNBC<-molecule %>% filter(ER=="Negative",PR=="Negative",her2=="Negative")
  TNBC_clin<-B_cli[B_cli$A0_Samples %in% TNBC$ID,]
  head(TNBC_clin)
  dim(TNBC_clin)
  
  TCGA_fpkm<-read.table("/media/data1/YJ/TCGA/TCGA-BRCA.htseq_fpkm.tsv.gz",header = T)
  ID2symbol<-read.table("/media/data1/YJ/TCGA/gencode.v22.annotation.gene.probeMap",header = T)
  head(ID2symbol)
  ID2symbol<-ID2symbol[,c(1,2)]
  colnames(ID2symbol)[1]<-"Ensembl_ID"
  TCGA_fpkm<-merge(TCGA_fpkm,ID2symbol)
  rownames(TCGA_fpkm)<-TCGA_fpkm$Ensembl_ID
  TCGA_fpkm$symbol<-TCGA_fpkm$gene
  TCGA_fpkm<-TCGA_fpkm[,c(colnames(TCGA_fpkm)[grepl("TCGA",colnames(TCGA_fpkm))],"symbol")]
  
  TCGA_fpkm<-TCGA_fpkm[TCGA_fpkm$symbol!="",]
  TCGA_fpkm<-TCGA_fpkm[!is.na(TCGA_fpkm$symbol),]
  TCGA_fpkm$median<-rowMedians(as.matrix(TCGA_fpkm[,1:1217]))
  TCGA_fpkm<-TCGA_fpkm[order(TCGA_fpkm$symbol,TCGA_fpkm$median,decreasing = T),]
  TCGA_fpkm<-TCGA_fpkm[!duplicated(TCGA_fpkm$symbol),]
  rownames(TCGA_fpkm)<-TCGA_fpkm$symbol
  TCGA_fpkm<-TCGA_fpkm[,1:1217]
  
  ### Screening out TNBC samples
  TCGA_fpkm_tnbc<-TCGA_fpkm[,colnames(TCGA_fpkm)[substring(colnames(TCGA_fpkm),1,12) %in% TNBC_clin$A0_Samples]]
  ### Screening out tumor samples
  TCGA_fpkm_tnbc_T<-TCGA_fpkm_tnbc[,substring(colnames(TCGA_fpkm_tnbc),14,15)<10]
  ## Perform GSVA analysis
  data<-as.data.frame(TCGA_fpkm_tnbc_T)
  data<-data[rowSums(data)>0,]
  
  TCGA_H_gsva<-gsva(as.matrix(data),h,method="gsva",kcdf="Gaussian",parallel.sz=32)
  data_cor <- rcorr(as.matrix(t(TCGA_H_gsva)))
  CorMatrix <- function(cor,p) {
    ut <- upper.tri(cor) 
    data.frame(row = rownames(cor)[row(cor)[ut]] ,
               column = rownames(cor)[col(cor)[ut]], 
               cor =(cor)[ut], 
               p = p[ut] )
  }
  data_cor<-CorMatrix(data_cor$r,data_cor$P)
  data_cor<-data_cor[data_cor$column %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  data_cor<-data_cor[!data_cor$row %in% c("KEGG_PROTEASOME","GO_PROTEASOME_ASSEMBLY"),]
  TCGA_data_cor<-data_cor[data_cor$column=="KEGG_PROTEASOME",]
  TCGA_data_cor<-TCGA_data_cor[order(TCGA_data_cor$row),]
  
  save(own_data_cor,TCGA_data_cor,file="./data/Hallmark_gene_set_cor_bulk.Rdata")
}

