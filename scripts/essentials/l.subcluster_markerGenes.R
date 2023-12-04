# MDLC1_scRNAseq_epi@active.ident = MDLC1_scRNAseq_epi$RNA_snn_res.1_new
# table(MDLC1_scRNAseq_epi@active.ident)
# 
# subcluster_markers <- FindAllMarkers(MDLC1_scRNAseq_epi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# subcluster_markers$absLogFC = abs(subcluster_markers$avg_log2FC)
# 
# 
# write.table(subcluster_markers, file = "../outputs/MDLC1_scRNAseq_subcluster marker genes.tsv",
#             row.names = F,sep = "\t", quote = F)



subcluster_markers = read.delim("../outputs/MDLC1_scRNAseq_subcluster marker genes.tsv")
subcluster_markers_sig = subset(subcluster_markers, p_val_adj <= 0.05 & absLogFC > 0.25)

table(subcluster_markers_sig$cluster)



##  all unique makers
subcluster_markers_sig -> topN
topN = topN[!duplicated(topN$gene),]

topN$histology = ifelse(topN$cluster %in% c("D1", "D2", "D3"), "Ductal", "Lobular")

MDLC1_scRNAseq_epi$labels = ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1_new %in% c("D1", "D2", "D3"), "Ductal", "Lobular")

MDLC1_scRNAseq_epi$order = ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 0, 1, # D3
                                  ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 1, 6, # L3
                                         ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 2, 2, # D2
                                                ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 3, 5, # L2
                                                       ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 4, 4, # L1
                                                              ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 5, 3, NA)))))) # D1

topN$order = ifelse(topN$cluster %in% "D3", 1,
                    ifelse(topN$cluster %in% "L3", 6,
                           ifelse(topN$cluster %in% "D2", 2,
                                  ifelse(topN$cluster %in% "L2", 5, 
                                         ifelse(topN$cluster %in% "L1", 4,
                                                ifelse(topN$cluster %in% "D1", 3, NA))))))


library(ComplexHeatmap)
features = topN$gene
mat = as.matrix(MDLC1_scRNAseq_epi@assays$RNA@data[features, ])

cols = c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")
colFun = circlize::colorRamp2(c(-2, -1, 0, 1, 2), cols, space = 'RGB')

set.seed(123); Markers_Ht = Heatmap(t(scale(t(mat))), show_column_names = F,show_row_names = F, 
                                    cluster_row_slices = F, name = "mRNA", cluster_column_slices = F, 
                                    row_names_gp = gpar(fontsize = 8,fontface="bold"), 
                                    use_raster = F, row_title = " ", column_title = " ",
                                    column_split = MDLC1_scRNAseq_epi$order, row_split = topN$order, 
                                    border = F, row_gap = unit(0,"cm"), column_gap = unit(0,'cm'),
                                    show_row_dend = F, show_column_dend = F, col = colFun, 
                                    top_annotation = HeatmapAnnotation(
                                      Cluster = MDLC1_scRNAseq_epi$RNA_snn_res.1_new,
                                      col = list(Histology = annot_col$Histology, 
                                                 Cluster = cluster_cols)))







## TOP 5 markers

subcluster_markers_sig %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = -p_val_adj) %>% as.data.frame() -> topN
split(topN$gene, topN$cluster)


topN$histology = ifelse(topN$cluster %in% c("D1", "D2", "D3"), "Ductal", "Lobular")

MDLC1_scRNAseq_epi$labels = ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1_new %in% c("D1", "D2", "D3"), "Ductal", "Lobular")

MDLC1_scRNAseq_epi$order = ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 0, 5, # D3
                                  ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 1, 3, # L3
                                         ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 2, 6, # D2
                                                ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 3, 2, # L2
                                                       ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 4, 1, # L1
                                                              ifelse(MDLC1_scRNAseq_epi$RNA_snn_res.1 %in% 5, 4,NA)))))) # D1

table(MDLC1_scRNAseq_epi$RNA_snn_res.1, MDLC1_scRNAseq_epi$RNA_snn_res.1_new)
topN$order = ifelse(topN$cluster %in% "D3", 5,
                    ifelse(topN$cluster %in% "L3", 3,
                           ifelse(topN$cluster %in% "D2", 6,
                                  ifelse(topN$cluster %in% "L2", 2, 
                                         ifelse(topN$cluster %in% "L1", 1,
                                                ifelse(topN$cluster %in% "D1", 4,NA))))))


library(ComplexHeatmap)
features = topN$gene
mat = as.matrix(MDLC1_scRNAseq_epi@assays$RNA@data[features, ])

cols = c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")
colFun = circlize::colorRamp2(c(-2, -1, 0, 1, 2), cols, space = 'RGB')

set.seed(123); topMarkers_Ht = Heatmap(t(scale(t(mat))), show_column_names = F,show_row_names = T, 
                                       cluster_row_slices = F, name = "mRNA", cluster_column_slices = F, 
                                       row_names_gp = gpar(fontsize = 8,fontface="bold"), 
                                       use_raster = F, row_title = " ", column_title = " ",
                                       column_split = MDLC1_scRNAseq_epi$order, row_split = topN$order, 
                                       border = F, row_gap = unit(0,"cm"), column_gap = unit(0,'cm'),
                                       show_row_dend = F, show_column_dend = F, col = colFun, 
                                       top_annotation = HeatmapAnnotation(Histology = MDLC1_scRNAseq_epi$CellTypes_new,
                                                                          Cluster = MDLC1_scRNAseq_epi$RNA_snn_res.1_new,
                                                                          col = list(Histology = annot_col$Histology, 
                                                                                     Cluster = cluster_cols)))



