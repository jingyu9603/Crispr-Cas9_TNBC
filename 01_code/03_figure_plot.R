## Fig1A
if(F){
  df<-G_gene_sum[G_gene_sum$`Down/total`!="NaN",]
  df$neg.p.value<- -log10(df$neg.p.value)
  df<-df[order(df$`Down/total`),]
  df$pos.p.value<- -log10(df$pos.p.value)
  
  tmp<-df[,c("neg.lfc","neg.p.value","total","Down/total")]
  tmp$neg.p.value<- -1*tmp$neg.p.value
  colnames(tmp)<-c("lfc","pvalue","total","recurrence")
  
  tmp1<-df[,c("neg.lfc","pos.p.value","total","Up/total")]
  colnames(tmp1)<-c("lfc","pvalue","total","recurrence")   
  tmp<-rbind(tmp,tmp1)
  
  p<-ggplot(tmp, aes(pvalue,lfc)) +
    geom_point(aes(size = recurrence),alpha=0.5,color="gray",)+ 
    geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
    geom_vline(aes(xintercept=1.3), colour="black", linetype="dashed")+
    geom_vline(aes(xintercept=-1.3), colour="black", linetype="dashed")+
    ylab("Log Fold change")+
    xlab("log10(Pos.pvalue)/log10(Neg.pvalue)")+
    geom_point(data=G_gene_sum_filtered_up,aes( -log10(pos.p.value),neg.lfc,size = Up/total,color=as.factor(total)),alpha=0.4)+
    geom_point(data=G_gene_sum_filtered,aes( log10(neg.p.value),neg.lfc,size = Down/total,color=as.factor(total)),alpha=0.4)+
    theme_bw()+
    # scale_color_gradient2(low = '#6B58EE7F', high = '#FF707F7F')+
    scale_color_brewer(palette="Set2")+
    scale_radius()+
    # ylim(c(1,13))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          legend.title = element_text(color = "black",size=18),
          legend.text  = element_text(color = "black",size=16),
          legend.position="bottom",
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(color="Detected sgRNA per gene",size="Recurrence")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  pdf(paste0("Fig1A_BT549_GeCKO_point_plot_",t,".pdf"),height = 7,width = 6.5)
  print(p)
  dev.off()
}

## Fig1B
if(F){
  G_gene_sum_filtered<-G_gene_sum[(G_gene_sum$neg.p.value<0.05) & (G_gene_sum$neg.lfc<0) & (G_gene_sum$total >=3 ) & (G_gene_sum$`Down/total`>0.5),]
  G_gene_sum_filtered<-G_gene_sum_filtered[order(G_gene_sum_filtered$`Down/total`),]
  G_gene_sum_filtered<-G_gene_sum_filtered[order(G_gene_sum_filtered$total),]
  
  df<-G_gene_sum[G_gene_sum$`Down/total`!="NaN",]
  df$neg.p.value<- -log10(df$neg.p.value)
  df<-df[order(df$total),]
  
  
  tmp<-df[,c("neg.lfc","neg.p.value","total","Down/total")]
  colnames(tmp)<-c("lfc","pvalue","total","Down_total")
  
  
  p<-ggplot(tmp, aes(pvalue,lfc)) +
    geom_point(aes(size=Down_total),alpha=0.4,color="gray",)+ 
    geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
    geom_vline(aes(xintercept=1.3), colour="black", linetype="dashed")+
    geom_point(data=G_gene_sum_filtered,aes( -log10(neg.p.value),neg.lfc,color = as.character(total),size=Down/total),alpha=0.5)+
    ylab("Log2 Fold change")+
    xlab("-Log10(p value)")+
    theme_bw()+
    scale_color_manual(values = pal_npg()(4)[c(2,3,4,1)])+
    # scale_color_gradient2(low = 'blue', high = 'red')+
    scale_radius()+
    # ylim(c(1,13))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          legend.title = element_text(color = "black",size=18),
          legend.text  = element_text(color = "black",size=16),
          legend.position="bottom",
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs( size="Down recurrence",color="Detected sgRNA per gene")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  pdf(paste0("./Fig1B_BT549_GeCKO_sig_down_gene_point_plot_",t,".pdf"),height = 7,width = 7)
  print(p)
  dev.off()
}

## Fig1C
if(F){
  G_gene_sum_zero_filtered<-G_gene_sum_zero[(G_gene_sum_zero$neg.p.value<0.05) & (G_gene_sum_zero$neg.lfc<0) & (G_gene_sum_zero$total >=3 ) & (G_gene_sum_zero$`Down/total`>0.5),]
  G_gene_sum_zero_filtered<-G_gene_sum_zero_filtered[order(G_gene_sum_zero_filtered$`Down/total`),]
  x = list(
    "Significant down regulated in AZD-Day7" = G_gene_sum_zero_filtered$id,
    "Efficiently expressed gene in BT549" = BT549_WT_rpkm$id,
    "Up-regulated in tumor tissue(bulk)"=rownames(res.sig),
    "Up-regulated in tumor cell(scRNA)"=rownames(Tepi.de.markers.sig),
    "Non-housekeeping gene" = non_essential_gene$id
  )
  library(UpSetR)
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  col_fun = colorRamp2(c(0,49), c("blue","red"))
  col<-col_fun(seq(1:49))
  names(col)<-sort(unique(comb_size(m1)))
  
  m1 = make_comb_mat(x,mode = "intersect",min_set_size = 36)
  # m1<-m1[comb_degree(m1)> 1]
  m1<-m1[comb_name(m1) %in% c(
    # "10000","01000","00100","00010","00001",
    "11000","11110","11111")]
  cs = comb_size(m1)
  pdf(paste0("./10_modified_figure/Fig1D_gene_upset_plot_high_",t,".pdf"),width = 6,height = 4.5)
  ht<-UpSet(m1,
            comb_order = order(-comb_size(m1)),
            # comb_col = col[as.character(comb_size(m1))],
            set_order = names(x),
            # right_annotation=NULL,
            right_annotation = upset_right_annotation(m1,gp = gpar(fill = pal_npg()(4),col="white")),
            top_annotation = HeatmapAnnotation(
              "Gene Intersections" = anno_barplot(cs, 
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE,
                                                  gp = gpar(fill = pal_nejm()(3),
                                                            col=pal_nejm()(3)), 
                                                  height = unit(5, "cm")
              ), 
              annotation_name_side = "left", 
              annotation_name_rot = 90))
  ht = draw(ht)
  od = column_order(ht)
  decorate_annotation("Gene Intersections", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(10, "pt"), 
              default.units = "native", just = c("bottom"), 
              gp = gpar(fontsize = 16, col = "#404040"))
  })
  
  dev.off()
}

## Fig1E
if(F){
  tmp<-G_gene_sum_zero_filtered
  tmp<-tmp[tmp$id %in% BT549_WT_rpkm$id,]
  dim(tmp)
  tmp<-tmp[tmp$id %in% non_essential_gene$id,]
  tmp<-tmp[order(-tmp$`Down/total`),]
  
  tmp<-tmp[tmp$id %in% rownames(res.sig),]
  # tmp<-tmp[tmp$id %in% rownames(SRP157974_res.sig),]
  tmp<-tmp[tmp$id %in% rownames(Tepi.de.markers.sig),]
  
  df<-tmp[order(tmp$neg.lfc),]
  head(df)
  dim(df)
  final52<-df
  
  df$neg.p.value<- -log10(df$neg.p.value)
  p<-ggplot(df, aes(neg.p.value,neg.lfc)) +
    geom_point(alpha=1,size=5,color="gray")+ 
    geom_hline(aes(yintercept=-1), colour="black", linetype="dashed")+
    # geom_vline(aes(xintercept=-1), colour="black", linetype="dashed")+
    ylab("Log2 Fold change")+
    xlab("-Log10(p value)")+
    geom_point(data=df[1:10,],aes(color=id),alpha=1,size=6)+ 
    geom_text_repel(data=df[1:10,],aes(label=id,color=id), size=7)+
    theme_bw()+
    scale_color_npg()+
    scale_radius()+
    # ylim(c(1,13))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  pdf(paste0("./Fig1D_Filtered_90gene_scatterplot_",t,".pdf"),family = "Arial",width = 6.5,height = 6)
  print(p)
  dev.off()
}

## supFig1B
if(F){
  p<-ggplot(BAGEL.bf, aes(log10_pvalue, Fitness_Score)) +
    geom_point(aes(color=Group),alpha=0.6,shape=19)+ 
    geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
    geom_vline(aes(xintercept=1.3), colour="black", linetype="dashed")+
    # geom_point(data=df_sub,color=pal_npg()(10))+
    # geom_text_repel(aes(label=id), size=4, data=df_sub)+
    # geom_vline(aes(xintercept=0.5), colour="black", linetype="dashed")+
    scale_color_manual(values = c("Essential"="red","Non-essential"="gray"))+
    ylab("Fitness score")+
    xlab("-Log10(p value)")+
    theme_bw()+
    # ylim(c(1,13))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(color = "black",size=18),
          legend.text  = element_text(color = "black",size=16))+
    labs(color="")
  p
  pdf(paste0("./10_modified_figure/Fig1C_BT549_essential_gene_point_plot_BAGEL_",t,".pdf"),height = 6,width = 6)
  print(p)
  dev.off()
}

## Fig3 A and B
if(F){
  ## Fig3A
  library(RobustRankAggreg)
  load("./data/Hallmark_gene_set_cor.Rdata")
  GSE11838_data_cor<-GSE11838_data_cor[order(GSE11838_data_cor$cor,decreasing = T),]
  GSE75688_data_cor<-GSE75688_data_cor[order(GSE75688_data_cor$cor,decreasing = T),]
  
  des_order.list<-list("GSE11838"=GSE11838_data_cor$row,
                       "GSE75688"=GSE75688_data_cor$row
  )
  des_order_matrix<-rankMatrix(des_order.list)
  des_orderAR=aggregateRanks(rmat=des_order_matrix)
  # colnames(des_orderAR)<-c("Name","Pvalue")
  des_orderAR$rank<- round(-log10(des_orderAR$Score),digits = 2)
  des_orderAR$Name<-factor(des_orderAR$Name,levels = rev(des_orderAR$Name))
  head(des_orderAR)
  
  GSE11838_data_cor<-GSE11838_data_cor[order(GSE11838_data_cor$row,decreasing = T),]
  GSE75688_data_cor<-GSE75688_data_cor[order(GSE75688_data_cor$row,decreasing = T),]
  data<-cbind(GSE11838_data_cor,GSE75688_data_cor)
  # data$GSE11838_padjust<-p.adjust(data$GSE11838_p,method = "BH")
  all(data[,1]==data[,5])
  colnames(data)<-c(paste0("GSE11838_",colnames(GSE11838_data_cor)),paste0("GSE75688_",colnames(GSE75688_data_cor)))
  data$min_p<-min(data$GSE11838_p,data$GSE75688_p)
  data$col<- -log10(data$GSE75688_p + 1E-10)
  data$szie<- -log10(data$GSE11838_p + 1E-10)
  data<-data[order(-data$GSE11838_cor),]
  head(data)
  sub_data<-data[data$GSE11838_row %in% des_orderAR$Name[1:5],]
  data$color<-ifelse(data$GSE11838_row %in% sub_data$GSE75688_row,data$GSE11838_row,"No")
  data$stroke<-ifelse(data$GSE11838_row %in% sub_data$GSE11838_row,1.2,0.2)
  sub_data$GSE11838_row<-factor(sub_data$GSE11838_row,levels = des_orderAR$Name[1:5])
  p<-ggplot(data, aes(x=GSE11838_cor, y=GSE75688_cor)) + 
    geom_point(aes(fill=col, size=szie,colour=color,stroke=stroke),shape=21) +
    theme_bw()+
    labs(
      y="Correlation coefficient(GSE75688)", 
      x="Correlation coefficient(GSE11838)"
    )+
    scale_fill_gradientn(name="-log10(p value)(GSE75688)",
                         colours= c("darkgreen","gray",pal_material(palette = "red",alpha = 1)(5)),
                         values=c(0,0.13,0.131,0.2,0.201,0.5,0.75,1),
                         limits=c(0,10),
                         breaks = c(0,1.3,2,5,7.5,10))+
    scale_colour_manual(values = c("black",pal_nejm()(10)[c(2,4,5)],pal_jama()(9)[c(3)],"white"),
                        breaks = c(levels(sub_data$GSE11838_row),"No"))+
    scale_size(name = "-log10(p value)(GSE11838)",
               range = c(0.8, 5),
               breaks = c(1.3, 2, 5, 7.5))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          legend.text = element_text(color = "black",size=16),
          legend.title = element_text(color = "black",size=18))
  
  pdf(paste0("./Fig3A_scRNAseq_proteasom_hallomark_cor_scatterplot_",t,".pdf"),width = 12,height = 6)
  print(p)
  dev.off()
  
  ## Fig3B
  # save(GSE11838_H_gsva,sc_H_gsva,file = "./data/scRNA_haoomark_gsva_res.Rdata")
  data<-cbind(GSE11838_H_gsva,sc_H_gsva)
  data<-data[c("HALLMARK_MTORC1_SIGNALING","KEGG_PROTEASOME"),]  
  data<-as.data.frame(t(data))
  data$group<-c(rep("GSE11838",ncol(GSE11838_H_gsva)),rep("GSE75688",ncol(sc_H_gsva)))
  table(data$group)
  
  display.brewer.pal(n = 10, name = 'Paired')
  p<-ggplot(data, aes(x=HALLMARK_MTORC1_SIGNALING, y=KEGG_PROTEASOME,col=group)) + 
    geom_point(aes(col=group, shape=group),alpha=1,size=2) +
    geom_smooth(aes(col=group), method="lm", se=F)+
    geom_rug()+
    theme_bw()+
    labs(
      x="mTORc1 signaling", 
      y="Proteasome"
    )+
    stat_cor(aes(color = group))+
    scale_color_manual(values = c("#386CB0","#F0027F"))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          legend.text = element_text(color = "black",size=16),
          legend.title = element_text(color = "black",size=18),
          legend.position="bottom",
          legend.box="vertical")+
    labs(color="Dataset")
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  pdf(paste0("./Fig3B_sc_RNAseq_proteasom_PI3Kpathway_cor_scatterplot_",t,".pdf"),width = 7,height = 7)
  print(p)
  dev.off()
}

## FigS4 A and B 
if(F){
  library(RobustRankAggreg)
  load("./data/Hallmark_gene_set_cor_bulk.Rdata")
  own_data_cor<-own_data_cor[order(own_data_cor$cor,decreasing = T),]
  TCGA_data_cor<-TCGA_data_cor[order(TCGA_data_cor$cor,decreasing = T),]
  
  des_order.list<-list("PRJNA553095"=own_data_cor$row,
                       "TCGA_TNBC"=TCGA_data_cor$row
  )
  des_order_matrix<-rankMatrix(des_order.list)
  des_orderAR=aggregateRanks(rmat=des_order_matrix)
  # colnames(des_orderAR)<-c("Name","Pvalue")
  des_orderAR$rank<- round(-log10(des_orderAR$Score),digits = 2)
  des_orderAR$Name<-factor(des_orderAR$Name,levels = rev(des_orderAR$Name))
  head(des_orderAR)

  own_data_cor<-own_data_cor[order(own_data_cor$row,decreasing = T),]
  TCGA_data_cor<-TCGA_data_cor[order(TCGA_data_cor$row,decreasing = T),]
  data<-cbind(own_data_cor,TCGA_data_cor)
  all(data[,1]==data[,5])
  colnames(data)<-c(paste0("PRJNA553095_",colnames(own_data_cor)),paste0("TCGA_",colnames(TCGA_data_cor)))
  
  data$col<- -log10(data$PRJNA553095_p + 1E-10)
  data$szie<- -log10(data$TCGA_p + 1E-10)
  data<-data[order(-data$PRJNA553095_cor),]
  head(data)
  sub_data<-data[data$PRJNA553095_row %in% des_orderAR$Name[1:5],]
  data$color<-ifelse(data$PRJNA553095_row %in% sub_data$PRJNA553095_row,data$PRJNA553095_row,"No")
  data$stroke<-ifelse(data$PRJNA553095_row %in% sub_data$PRJNA553095_row,1.2,0.2)
  sub_data$PRJNA553095_row<-factor(sub_data$PRJNA553095_row,levels = des_orderAR$Name[1:5])
  p<-ggplot(data, aes(x=PRJNA553095_cor, y=TCGA_cor)) + 
    geom_point(aes(fill=col, size=szie,colour=color,stroke=stroke),shape=21) +
    theme_bw()+
    labs(
      y="Correlation coefficient(TCGA-TNBC)", 
      x="Correlation coefficient(PRJNA553095)"
    )+
    scale_fill_gradientn(name="-log10(p value)(PRJNA553095)",
                         colours=c("darkgreen","gray",pal_material(palette = "red",alpha = 1)(5)),
                         values=c(0,0.13,0.131,0.2,0.201,0.5,0.75,1),
                         limits=c(0,10),
                         breaks = c(0,1.3,2,5,7.5,10))+
    scale_colour_manual(values = c("black",pal_nejm()(10)[c(2,4,5)],pal_jama()(9)[c(3)],"white"),
                        breaks = c(levels(sub_data$PRJNA553095_row),"No"))+
    scale_size(name = "-log10(p value)(TCGA-TNBC)",
               range = c(0.8, 5),
               breaks = c(1.3, 2, 5, 7.5))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          legend.text = element_text(color = "black",size=16),
          legend.title = element_text(color = "black",size=18))
  
  pdf(paste0("./FigS4A_bulk_RNAseq_proteasom_hallomark_cor_scatterplot_",t,".pdf"),family = "Arial",width = 12,height = 6)
  print(p)
  dev.off()
  
  ## bulk RNAseq level
  load("./data/bulk_hallmark_gsva_res.Rdata")
  data<-cbind(own_H_gsva,TCGA_H_gsva)
  data<-data[c("HALLMARK_MTORC1_SIGNALING","KEGG_PROTEASOME"),]  
  data<-as.data.frame(t(data))
  data$group<-c(rep("PRJNA553095",ncol(own_H_gsva)),rep("TCGA-TNBC",ncol(TCGA_H_gsva)))
  table(data$group)
  
  display.brewer.pal(n = 10, name = 'Paired')
  p<-ggplot(data, aes(x=HALLMARK_MTORC1_SIGNALING, y=KEGG_PROTEASOME,col=group)) + 
    geom_point(aes(col=group, shape=group),alpha=1,size=2) +
    geom_smooth(aes(col=group), method="lm", se=F)+
    geom_rug()+
    theme_bw()+
    labs(
      x="mTORc1 signaling", 
      y="Proteasome"
    )+
    stat_cor(aes(color = group))+
    scale_color_manual(values = c("#FB9A99","#33A02C"))+
    theme(axis.text = element_text(color = "black",size=20),
          axis.title = element_text(color = "black",size=25,face = "plain"),
          legend.text = element_text(color = "black",size=16),
          legend.title = element_text(color = "black",size=18),
          legend.position="bottom",
          legend.box="vertical")+
    labs(color="Dataset")
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  pdf(paste0("/Fig4SB_bulk_RNAseq_proteasom_PI3Kpathway_cor_scatterplot_",t,".pdf"),family = "Arial",width = 7,height = 7)
  print(p)
  dev.off()
}