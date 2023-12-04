# LOAD PACKAGES
library(magrittr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(EnvStats)
library(ggpubr)
library(gridExtra)
library(ggthemes)
library(reshape2)
library(ComplexHeatmap)

`%notin%`=Negate(`%in%`)
options("scipen"=100, "digits"=4)


# DEFINE ANNOTATIONS COLORS
pam50_cols = RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,5,3,4)]
names(pam50_cols)=c("LumA","LumB","Basal","Her2","Normal")

annot_col = list(Histology = c("Ductal"="blue","Lobular"="red"),
                 ER = circlize::colorRamp2(c(-1,0,1), c("#1B5E20", "black", "#f44336")),
                 CDH1 = circlize::colorRamp2(c(-1,0,1), c("#1B5E20", "black", "#f44336")),
                 PAM50 = pam50_cols,
                 Patient = c("MDLC-2"="#8FAADC","MDLC-3"="#F4B183","MDLC-1"="#C19100"))


cols = c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")
colFun = circlize::colorRamp2(c(-2, -1, 0, 1, 2), cols, space = 'RGB')


# LOAD FUNCTIONS
lfc_2_geneList = function(lfc_mat, geneNames, convert = F, fromType = "SYMBOL", toType = "ENTREZID"){
  geneList = lfc_mat
  names(geneList)=geneNames
  geneList = geneList[order(geneList,decreasing = T)]
  #geneList = geneList[abs(geneList) > 1] # lfc cutoff
  geneList_symbols = geneList
  
  if (convert){
    # Convert genes to entrez
    map = gene_name_conversion(geneNames, fromType = fromType, toType = toType)
    map = map[!duplicated(map[,1]),]
    rownames(map)=map[,1]
    names(geneList)=map[names(geneList),2]
    geneList = geneList[!is.na(names(geneList))]
    
    geneList_df = list(geneList = unlist(geneList), geneList_symbols = unlist(geneList_symbols))
    names(geneList_df)=c(toType, fromType)
    return(geneList_df)
  }
  return(geneList_symbols)
}

gene_name_conversion = function(geneNames, fromType = "SYMBOL", toType = "ENTREZID"){
  load("../../inputs/bioMart_geneMap.Rdata")
  library(org.Hs.eg.db)
  library(clusterProfiler)
  map = bitr(geneID = geneNames, fromType = fromType, toType = toType, OrgDb = org.Hs.eg.db)
  map = map[!duplicated(map[,1]),]
  rownames(map)=map[,1]
  geneNames = intersect(map[,1],geneNames)
  return(map)
}

mkplotDF <- function(dim_res, samples){
  plt_df <- data.frame(DIM_1 = dim_res[,1], DIM_2 = dim_res[,2], 
                       labels = colnames(mDLC_DSP_log2[,samples]),
                       Histology = factor(ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "SegmentLabel"]),
                       ROILabel = ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "ROILabel"],
                       Patient = patient_map[ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "SampleID"]],
                       AOIarea = log10(ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "AOISurfaceArea"]),
                       LOQ_01 = ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "LOQ_01"],
                       ID = ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "ROILabel"],
                       NucleiCnt = ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "AOINucleiCount"],
                       AlignedReads = log10(ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "AlignedReads"]+1),
                       NormFactorQ3_postDrop = ROI_annot[mDLC_DSP_log2[,samples] %>% colnames(), "NormFactorQ3_postDrop"])
  return(plt_df)
}

mk_gene_plot = function(gene, sampleIDs){
  df = data.frame(exp = as.numeric(mDLC_DSP_log2[gene, sampleIDs]),
                  hist = ROI_annot[mDLC_DSP_log2[,sampleIDs] %>% colnames(), "SegmentLabel"],
                  ID = ROI_annot[mDLC_DSP_log2[, sampleIDs] %>% colnames(), "ROILabel"],
                  Patient = patient_map[ROI_annot[mDLC_DSP_log2[,sampleIDs] %>% colnames(), "SampleID"]])
  
  
  plt=ggplot(df, aes(hist, exp, color = hist)) + geom_boxplot(alpha = 0.5, outlier.color = NA) + 
    ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 4, shape = 20) + 
    myTheme + scale_color_manual("Histology", values = annot_col$Histology) + 
    stat_compare_means(comparisons = list(c("Lobular","Ductal")), label = "p.signif", 
                       color = 'red', size = 6) +  #, vjust = 1.6
    ylab(paste0(gene, " expression")) + xlab("") 
  return(plt)
}

FeatureModulePlot=function(obj, geneList, geneListName, reduction="tsne", pt.size = 2, backgroundCtrl=200, cols, label=F,removeAxis=T){
  obj <- AddModuleScore(obj, features = list(geneList), ctrl = backgroundCtrl, name = geneListName,seed = 123)
  obj[[geneListName]]=obj[[paste0(geneListName,"1")]]
  plt=FeaturePlot(obj, features = geneListName,reduction = reduction, pt.size = pt.size, slot = "scale.data",
                  label = label, cols = cols,combine = !removeAxis)
  if (removeAxis){
    for(i in 1:length(plt)) {
      plt[[i]] <- plt[[i]] + NoAxes()
    }
    plt=cowplot::plot_grid(plotlist = plt)
  }
  return(plt)
}

library(limma)

run_Limma = function(fit, contrast, mat, adjust.method="BH", fileName=NULL){
  fit1 <- contrasts.fit(fit, contrast)
  fit1 <- eBayes(fit1, robust = T)
  res <- topTable(fit1, number = nrow(mat), coef=colnames(fit1$coefficients), adjust.method=adjust.method)
  res$absLogFC=abs(res$logFC)
  
  res$label = rownames(res)
  
  res$log10Pvalue = -log10(res$P.Value)
  res$log10Padj = -log10(res$adj.P.Val)
  
  return(res)
}

volcanoPlt = function(res, sig_res, pcut = 0.01, fc_cut = 1, topN = 15, annotsize = 4){
  topUP = subset(sig_res, logFC >= fc_cut & adj.P.Val <= pcut)
  dim(topUP)
  topUP = topUP[order(topUP$log10Padj,topUP$absLogFC,  decreasing = T),] %>% head(topN)
  
  topDn = subset(sig_res, logFC <= -fc_cut & adj.P.Val <= pcut)
  dim(topDn)
  topDn = topDn[order(topDn$log10Padj,topDn$absLogFC, decreasing = T),] %>% head(topN)
  
  plt = ggplot(res, aes(x = logFC, y = log10Padj)) + myTheme + 
    geom_point(data = subset(res, adj.P.Val <= pcut & abs(logFC) >= fc_cut), aes(color = factor(logFC > 0), size = log10Padj), alpha=0.8) + 
    geom_point(data = subset(res, adj.P.Val > pcut | abs(logFC) < fc_cut), color = "gray", alpha=0.8, show.legend = F) + 
    scale_color_manual("Log2 FC",values = c("blue", "red"), labels = c("Dn", "Up")) + 
    guides(color = guide_legend(override.aes = list(size=6)), size="none") + xlab("Log2 FC") + 
    geom_hline(yintercept = -log10(pcut), linetype="dashed", col = "darkgray") +
    geom_vline(xintercept = fc_cut, linetype="dashed", col = "darkgray") + 
    geom_vline(xintercept = -fc_cut, linetype="dashed", col = "darkgray") + 
    ggrepel::geom_label_repel(data = rbind(topUP, topDn), max.overlaps = 100, fontface = "bold", color = "black",
                              aes(label=label), show.legend = F, box.padding = 0.5) + 
    annotate(geom = "text", x = 0.8, y = 0, size=annotsize, fontface="bold",color="darkgray", 
             label = paste0("Lobular (",as.numeric(table(sig_res$logFC>0)[2]),")")) + 
    annotate(geom = "text", x = -0.8, y = 0, size=annotsize, fontface="bold",color="darkgray",
             label= paste0("Ductal (",as.numeric(table(sig_res$logFC>0)[1]),")")) + 
    theme(legend.position = "Dn")
  return(plt)
}

myOverlapPlt=function(mat,size=6){
  mat=sapply(mat, function(x){
    sapply(mat, function(y) length(intersect(x, y))/(length(c(x,y))-length(intersect(x, y))), simplify = F)})
  mode(mat)="numeric"
  df=reshape2::melt(mat)
  df$log_value=log2(df$value+1)
  df$value_simple=df$value
  df[round(df$value_simple,1) == 1,"value_simple"]=NA
  plt=ggplot(df, aes(x=Var1, y=Var2, fill=value_simple))+
    geom_tile(alpha=0.5, show.legend = T, color="black") + scale_fill_gradient("Jaccard\nSimilarity", low = "white", high = "#543005") + 
    myTheme+xlab("")+ylab("")+
    geom_text(show.legend = F, size=size, aes(label=round(value_simple,2)))+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    coord_flip()+blankTheme + 
    guides(color = FALSE, size = FALSE)
  return(plt)
}


enrichmentAnalysis = function(sig_score_mat, estimateCut = 0.15, groupA_ids, groupB_ids){
  #### Pvalue
  selectPathways = rownames(sig_score_mat) #intersect(rownames(sig_score), names(Hallmark$genesets))
  ductal = sig_score_mat[selectPathways, groupA_ids]
  lobular = sig_score_mat[selectPathways, groupB_ids]
  
  pvalues = sapply(selectPathways, function(x){
    t.test(ductal[x,], lobular[x,])$p.value
  })
  
  fdr = p.adjust(p = pvalues)
  
  estimate = rowMeans(lobular) - rowMeans(ductal)
  
  plot(density(estimate))
  mean(estimate) + sd(estimate)
  mean(estimate) - sd(estimate)
  
  group = ifelse(estimate > estimateCut, "Lobular", 
                 ifelse(estimate < -estimateCut, "Ductal", "None"))
  
  table(group)
  
  pvalues_simple = pvalues
  pvalues_simple[pvalues > 0.05] = "ns"
  pvalues_simple[pvalues <= 0.05] = "*"
  pvalues_simple[pvalues <= 0.01] = "**"
  pvalues_simple[pvalues <= 0.001] = "***"
  pvalues_simple[pvalues <= 0.0001] = "****"
  
  pvalues_cap = pvalues
  pvalues_cap[pvalues <= 0.0001] = 0.0001
  pvalues_cap[pvalues > 0.05] = 1
  
  fdr_cap = fdr
  fdr_cap[fdr <= 0.0001] = 0.0001
  fdr_cap[fdr > 0.05] = 1
  
  fdr_simple = fdr
  fdr_simple[fdr > 0.05] = "ns"
  fdr_simple[fdr <= 0.05] = "*"
  fdr_simple[fdr <= 0.01] = "**"
  fdr_simple[fdr <= 0.001] = "***"
  fdr_simple[fdr <= 0.0001] = "****"
  
  return(list(fdr_simple = fdr_simple, fdr = fdr, fdr_cap = fdr_cap, pvalues_cap = pvalues_cap, pvalues = pvalues, 
              pvalues_simple = pvalues_simple, estimate = estimate, ductal = rowMeans(ductal), lobular = rowMeans(lobular), group = group))
}

myWordCloud_terms = function(terms, wordstoignore, minFreq = 2){
  wordCloud_df = as.data.frame(table(unlist(stringr::str_split(terms,pattern = " "))))
  wordCloud_df = subset(wordCloud_df, Var1 %notin% wordstoignore)
  set.seed(123); wordcloud::wordcloud(words = as.character(wordCloud_df$Var1), wordCloud_df$Freq, min.freq = minFreq,
                                      colors = c(RColorBrewer::brewer.pal(10, "Dark2"), RColorBrewer::brewer.pal(10, "Set1")))
  plt <- recordPlot()
  return(plt)
}


enrichedInAnnot = function(ids, which = "column",label_pval = T, pvals = pvalues_cap, pchs = pvalues_simple){
  HeatmapAnnotation(#`Enriched In` = group[ids], 
                    `Score` = scale(estimate)[ids,],
                    `Pval` = anno_simple(x = floor(-log10(pvals[ids])),  
                                         pch = if(label_pval == T){
                                           pchs[ids]
                                         }else{
                                           rep(" ", length(ids))
                                         }, 
                                         col = c("0" = "#F9FBE7","1"= "#fffac6", "2" = "#FFF176", "3" = "#cddc39", "4" = "#388e3c")),
                    which = which,
                    col = list(`Enriched In` = c("Lobular" = "red", "Ductal" = "blue", "None" = "gray"),
                               `Score` = circlize::colorRamp2(c(-2, 0,  2), c("blue","white","red"), space = 'RGB')))}


colFun3 = circlize::colorRamp2(c(0,0.5,1), colors = c("white","#ca9446", "#543005"), space = 'RGB')

bot_annot_left = function(sampleIDs){
  HeatmapAnnotation(annotation_name_side = "left", 
                    Histology = ROI_annot[sampleIDs, "SegmentLabel"], 
                    Patient = patient_map[ROI_annot[sampleIDs, "SampleID"]],
                    annotation_name_gp = gpar(fontsize = 10, fontface="bold", col="#616161"),
                    col = annot_col)}

makeHeatmap = function(goi = goi, sig_res=sig_res, goi_group = goi_group, sampleIDs = sampleIDs, 
                       show_genes = F, show_column_names = T){
  mat = t(scale(t(mDLC_DSP_log2[goi, sampleIDs])))
  dim(mat)
  
  set.seed(123);  top_DEG_ht = Heatmap(mat, col = colFun, name = "Log2 Exp",# column_km = 3,
                                       show_column_names = show_column_names, show_row_names = show_genes, 
                                       column_labels = ROI_annot[sampleIDs, "ROILabel"], 
                                       show_column_dend = F, show_row_dend = F,border = T, 
                                       column_names_gp = gpar(col = ifelse(ROI_annot[sampleIDs, "ROILabel"] == 17, "red", "black"), fontsize = 12),
                                       row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                       top_annotation = bot_annot_left(sampleIDs))
  
  log2MAT = as.data.frame(sig_res[rownames(mat), "logFC"])
  rownames(log2MAT) = sig_res[rownames(mat),]$label
  colnames(log2MAT) = "log2FC"
  colFun1 = circlize::colorRamp2(c(min(log2MAT), 0, max(log2MAT)), 
                                 colors = c("blue",'white',"red"), space = 'RGB')
  log2HT = Heatmap(as.matrix(log2MAT), name = "Log2FC", col = colFun1, show_row_names = F,border = T,
                   width = unit(0.5,"cm"),
                   show_column_dend = F, show_row_dend = F)
  
  pvalMAT = as.data.frame(sig_res[rownames(mat), "log10Padj"])
  rownames(pvalMAT) = sig_res[rownames(mat),]$label
  colnames(pvalMAT) = "log-pval"
  colFun2 = circlize::colorRamp2(c(min(pvalMAT), max(pvalMAT)), 
                                 colors = c('white','red'), space = 'RGB')
  pvalHT = Heatmap(as.matrix(pvalMAT), name = "log-pval", col = colFun2, row_split = goi_group, show_row_names = F, 
                   width = unit(0.5,"cm"),
                   show_column_dend = F, border = T, show_row_dend = F)
  
  return(list(pvalHT, log2HT, top_DEG_ht))
}


# THEMES
myTheme=theme_bw()+theme(# TEXT
  text = element_text(size=15),
  axis.text = element_text(size=12, face = "bold", color = "#1A237E"),
  legend.text = element_text(size=12, face = "bold"),
  plot.subtitle = element_text(size=8,face="bold"),
  # BACKGROUND
  legend.background = element_rect(fill = "white", size = 4, colour = "white"),
  panel.grid.major = element_blank(),#(color = "#9E9E9E",linetype = "dotted", size = 0.1),
  panel.grid.minor = element_line(color="#BDBDBD",linetype = "dotted"),
  # TICKS
  #LEGEND POS
)


myTheme_barplot=theme_bw()+theme(# TEXT
  text = element_text(size=15, family = "Helvetica",face = "bold"),
  legend.text = element_text(size = 10,face="bold"),
  legend.title = element_text(size=12,face="bold"),
  rect=element_blank(),
  #rect = element_rect(linetype = "dotted",colour = "grey"),
  # axis.title = element_text(size = 20),
  axis.line.y = element_line(colour = "#424242",size = 0.5),
  #axis.line.x = element_blank(),
  plot.title = element_text(hjust = 0.5),
  # line = element_blank(),
  panel.background = element_rect(fill =  "white"),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank()
  # TICKS
  #LEGEND POS
)

blankTheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                   panel.background = element_blank(), axis.ticks.x = element_blank(),plot.margin=grid::unit(c(0,0,0,0), "mm"))


myTheme_umap=theme_bw()+theme(# TEXT
  text = element_text(size=15),
  legend.text = element_text(size = 10,face="bold"),
  legend.title = element_text(size=12,face="bold"),
  rect=element_blank(),
  # rect = element_rect(linetype = "dotted",colour = "grey"),
  # axis.title = element_text(size = 20),
  axis.line = element_line(linetype = "dashed",colour = "gray",#"#0091EA",
                           size = 0.5,
                           arrow=ggplot2::arrow(angle=30,length = unit(0.24, "cm"), type = "closed")),
  plot.title = element_text(hjust = 0.5),
  # line = element_blank(),
  # panel.background = element_rect(colour = "white",color = "white"),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank()
  # TICKS
  #LEGEND POS
)
