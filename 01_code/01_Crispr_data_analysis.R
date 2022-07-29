##=============================================================
## 1、Filter essential gene
if(F){
  lib_gene_summary<-read.table("/media/data1/YJ/WXY_TNBC/mageck_out/BT549_Human_GeCKo_25th.gene_summary.txt",header = T)
  
  ## Calculate Fitness score using BAGEL
  ### In preparation for the BAGEL analysis
  load("/media/data1/YJ1/softwave/BAGELR-master/data/BAGEL_v2_ESSENTIAL_GENES.rdata")
  load("/media/data1/YJ1/softwave/BAGELR-master/data/BAGEL_nonEssential.rdata")
  lib_sgRNA_summary<-read.table("//media/data1/YJ/WXY_TNBC/mageck_out/BT549_Human_GeCKo_25th.sgrna_summary.txt",header=T)
  head(lib_sgRNA_summary)
  df<-lib_sgRNA_summary[,c(1,2,7)]
  head(df)
  colnames(df)<-c("SEQID","GENE","BT549")
  
  refGuidesLibrary<-read.csv("/data/YJ/LYC/Crispr/Reference/human_geckov2_library_ab_09032015.csv",header = F)
  head(refGuidesLibrary)
  refGuidesLibrary<-refGuidesLibrary[,c(1,3)]
  write.table(df,file="../BAGEL_analysis/BT549_input.txt",quote = F,row.names = F,sep = "\t")
  write.table(BAGEL_essential,file = "../BAGEL_analysis/BAGEL_essential.txt",quote = F,row.names = F,sep = "\t")
  write.table(BAGEL_nonEssential,file = "../BAGEL_analysis/BAGEL_nonEssential.txt",quote = F,row.names = F,sep = "\t")
  
  ## Load the results of the BAGEl analysis
  BAGEL.bf<-read.table("/media/data1/YJ/WXY_TNBC/BAGEL_analysis/BT549_bagel_res.bf",header = T)
  BAGEL.bf$Fitness_Score<- -1*BAGEL.bf$BF
  colnames(BAGEL.bf)[1]<-"id"
  head(BAGEL.bf)
  
  head(lib_gene_summary)
  BAGEL.bf<-merge(BAGEL.bf,lib_gene_summary)
  
  ## Screen out the essential gene
  BAGEL.bf$Group<-ifelse((BAGEL.bf$neg.p.value<0.05) & (BAGEL.bf$Fitness_Score< -0),"Essential"," ")
  table(BAGEL.bf$Group)
  non_essential_gene<-BAGEL.bf[BAGEL.bf$Group !="Essential",]
  essential_gene<-BAGEL.bf[BAGEL.bf$Group=="Essential",]
}

##===============================================================
## 2、Screening RPKM>10, it is necessary to screen out genes with a certain amount of RNA expression
if(F){
  BT549_rpkm<-read.table("./data/GSE112365_all_samples_rpkm_and_RC.txt",header = T,sep = "\t")
  head(BT549_rpkm)
  BT549_WT_rpkm<-BT549_rpkm[,c(1,15,17,19)]
  BT549_WT_rpkm$contrlmeanFPKM<-rowMeans(BT549_WT_rpkm[,2:4])
  BT549_WT_rpkm<-BT549_WT_rpkm[BT549_WT_rpkm$contrlmeanFPKM>10,]
  rownames(BT549_WT_rpkm)<-BT549_WT_rpkm$Gene.name
  head(BT549_WT_rpkm)
  
  BT549_WT_rpkm<-BT549_WT_rpkm[,c(1,5)]
  colnames(BT549_WT_rpkm)[1]<-"id"
  head(BT549_WT_rpkm)
}

##=============================================================
## 3、Screen GeCKO AZD vs. DMSO group
if(F){
  ## Load the results of the analysis using MAGeCK
  G_gene_sum<-read.table("/media/data1/YJ/WXY_TNBC/mageck_out/BT549_GeCKo_25th.gene_summary.txt",header = T)
  head(G_gene_sum)
  G_gene_sum[G_gene_sum$id=="PSMG2",]
  G_sgrna_sum<-read.table("/media/data1/YJ/WXY_TNBC/mageck_out/BT549_GeCKo_25th.sgrna_summary.txt",header = T)
  head(G_sgrna_sum)
  G_sgrna_sum[G_sgrna_sum$Gene=="PSMG2",]
  ## Any sgRNA with expression value greater than 25 in control or treatment was screened out
  G_sgrna_sum<-filter_at(G_sgrna_sum, vars(control_count,treatment_count), any_vars(.>25)) 
  G_sgrna_sum<-G_sgrna_sum[order(G_sgrna_sum$control_count,G_sgrna_sum$treatment_count),]
  
  ## The results of sgRNA were mapped to the corresponding gene
  df<-matrix(nrow = length(unique(G_sgrna_sum$Gene)),ncol = 2)
  rownames(df)<-unique(G_sgrna_sum$Gene)
  colnames(df)<-c("Up","Down")
  head(df)
  for(gene in rownames(df)){
    tmp<-G_sgrna_sum[G_sgrna_sum$Gene==gene,,drop=F]
    up=0
    down=0
    for(i in 1:nrow(tmp)){
      if(tmp[i,7] > 0){
        up=up+1
      }else{
        if(tmp[i,7] < 0){
          down=down+1
        }
      }
    }
    df[gene,1]<-up
    df[gene,2]<-down
  }
  
  df<-as.data.frame(df)
  df$total<-df$Up+df$Down
  df$"Down/total"<-df$Down/df$total
  df$id<-rownames(df)
  head(df)
  Gecko_sgRNA2gene_zero<-df
  
  G_gene_sum_zero<-merge(G_gene_sum,df)
  G_gene_sum_zero$`Up/total`<-G_gene_sum_zero$Up/G_gene_sum_zero$total
}
