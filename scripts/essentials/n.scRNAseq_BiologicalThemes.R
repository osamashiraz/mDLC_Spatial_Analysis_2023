scRNAseqMat = as.matrix(MDLC1_scRNAseq_epi@assays$RNA@data)


# MSIG DB
gset_ = c(GO_BP$genesets, GO_CC$genesets, GO_MF$genesets, Hallmark$genesets, KEGG$genesets, REACTOME$genesets)
gset_ = gset_[grep(x = names(gset_), pattern = "Collagen|ECM|Extracellular Matrix|Ncam1|Immunoglobulin|Mediated Immunity|Leptin|Glucose Metaboli|Gluconeogenesis|Mitotic|Cell Cycle|Senescence|Proteoglycan|Tp53|Proteosome|Mhc Class I|Antigen Processing|Oxidative Phosphorylation|Electron Transport Chain|Respiration|Sumoylation|Acetyltransferase|Oxidative Stress|Double Strand|DNA Repair|Mtor|Phosphatidylinositol|Protein Maturation|Unfolded Protein Response|Mapk Targets")]


gset_ = gset_[biological_themes_df$pathway] # pathways part of biological themes


DSP_Significant_BiologicalThemes = read.delim("../outputs/biologicalThemes_DSP.tsv")
rownames(DSP_Significant_BiologicalThemes) = DSP_Significant_BiologicalThemes$pathway

lobular_DSP_Significant_BT = subset(DSP_Significant_BiologicalThemes, group %in% "Lobular" & pvalues <= 0.05 & bg_pct < 5)$pathway
ductal_DSP_Significant_BT = subset(DSP_Significant_BiologicalThemes, group %in% "Ductal" & pvalues <= 0.05 & bg_pct < 5)$pathway


gset_scRNA = sapply(gset_, function(x){
  x = setdiff(x, "CDH1")
  return(intersect(x, rownames(scRNAseqMat)))
})

degs_scRNA = subset(MDLC1_scRNAseq_epi_markers, absLogFC > 0.25 & p_val_adj < 0.1) %>% rownames()

gset_sub_scRNA = sapply(gset_scRNA, function(x){
  return(intersect(x, degs_scRNA))
})

sort(sapply(gset_sub_scRNA, length), decreasing = T)
names(sapply(gset_sub_scRNA, length))[sapply(gset_sub_scRNA, length) >= 3] -> selectPathways

selectPathways = intersect(selectPathways, unique(c(lobular_DSP_Significant_BT,
                                                    ductal_DSP_Significant_BT)))
intersect(selectPathways, names(gset_))

length(gset_sub_scRNA)
gset_sub_scRNA = gset_sub_scRNA[selectPathways]
length(gset_sub_scRNA)


library(GSVA)
set.seed(123); sig_score_scRNAseq = gsva(expr = as.matrix(scRNAseqMat), 
                                         gset.idx.list = gset_sub_scRNA,
                                         method="ssgsea", kcdf = "Poisson", parallel.sz = 1)


sig_score_df = reshape2::melt(t(scale(t(sig_score_scRNAseq))))
colnames(sig_score_df) = c("Pathway", "Cell", "Score")
head(sig_score_df)
sig_score_df$Cluster = as.character(MDLC1_scRNAseq_epi@active.ident[sig_score_df$Cell])
MDLC1_scRNAseq_epi@active.ident = MDLC1_scRNAseq_epi$labels %>% as.factor()
sig_score_df$Histology = as.character(MDLC1_scRNAseq_epi@active.ident[sig_score_df$Cell])
MDLC1_scRNAseq_epi@active.ident = as.factor(MDLC1_scRNAseq_epi$RNA_snn_res.1_new)


selectPathways = rownames(sig_score_scRNAseq)

# D vs L - pvalues
ductal = sig_score_scRNAseq[selectPathways, unique(as.character(subset(sig_score_df, Histology %in% "Ductal")$Cell))]
lobular = sig_score_scRNAseq[selectPathways, unique(as.character(subset(sig_score_df, Histology %in% "Lobular")$Cell))]
pvalues = sapply(selectPathways, function(x){
  #print(x)
  t.test(ductal[x,], lobular[x,])$p.value
})

table(pvalues < 0.05)

fdr = p.adjust(p = pvalues)
estimate = rowMedians(lobular) - rowMedians(ductal)
names(estimate) = names(pvalues)

mean(estimate) + 2*sd(estimate)
mean(estimate) - 2*sd(estimate)

table(abs(estimate) > 0.05)
estimate[abs(estimate) > 0.05]

group = ifelse(estimate > 0.05, "Lobular", 
               ifelse(estimate < -0.05, "Ductal", "None"))

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


group_sub = group[group %notin% "None"]
pvalues_sub = pvalues[names(group_sub)]

names(pvalues_sub)[pvalues_sub <= 0.05] -> SigPathways
length(SigPathways)

group_sig = group[SigPathways]

myDSP = c(intersect(lobular_DSP_Significant_BT, names(group_sig[group_sig == "Lobular"])),
          intersect(ductal_DSP_Significant_BT, names(group_sig[group_sig == "Ductal"])))



cols = c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")
colFun4 = circlize::colorRamp2(c(-1, -0.25, 0, 0.25, 1), cols, space = 'RGB')

cluster_cols = c("D1" = "#90CAF9", "D2" = "#2196F3", "D3" = "#0D47A1", "L1" = "#EF9A9A", "L2" = "#F44336", "L3" = "#B71C1C")

colFun3 = circlize::colorRamp2(c(0,0.5,1), colors = c("white","#ca9446", "#543005"), space = 'RGB')

myDSP_sig = na.omit(names(pvalues_sub[myDSP][pvalues_sub[myDSP] <= 0.05]))
set.seed(123);  biologicalTheme_heatmap_scRNA = Heatmap(t(scale(t(sig_score_scRNAseq[myDSP_sig, colnames(MDLC1_scRNAseq_epi)]))), 
                                                      col = colFun, name = "GSVA", #width = unit(10,'cm'),
                                                      show_column_names = F, show_row_names = F,  border = T, row_gap = unit(0,"cm"), 
                                                      column_gap = unit(0,'cm'), row_title_rot = 0,
                                                      row_labels = gsub(x = myDSP_sig, pattern = "Hallmark - |KEGG - ", 
                                                                        replacement = ""),
                                                      left_annotation = enrichedInAnnot(myDSP_sig, label_pval = F, which = "row"),
                                                      column_split = MDLC1_scRNAseq_epi@meta.data[, c("labels")],
                                                      row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                                      row_split = DSP_Significant_BiologicalThemes[myDSP_sig, "Clusters"], 
                                                      
                                                      cluster_column_slices = F,
                                                      top_annotation = HeatmapAnnotation(#Cluster = MDLC1_scRNAseq_epi$RNA_snn_res.1_new,
                                                        Histology = MDLC1_scRNAseq_epi$CellTypes_new,
                                                        col = list(Histology = annot_col$Histology, 
                                                                   Cluster = cluster_cols)))








# estimate
Cluster_Cells = split(colnames(MDLC1_scRNAseq_epi), MDLC1_scRNAseq_epi@meta.data$RNA_snn_res.1_new)
table(MDLC1_scRNAseq_epi@meta.data$RNA_snn_res.1_new, MDLC1_scRNAseq_epi@meta.data$labels)

library(matrixStats)
D1 = rowMedians(sig_score_scRNAseq[,  Cluster_Cells$D1])
D2 = rowMedians(sig_score_scRNAseq[,  Cluster_Cells$D2])
D3 = rowMedians(sig_score_scRNAseq[,  Cluster_Cells$D3])
L1 = rowMedians(sig_score_scRNAseq[,  Cluster_Cells$L1])
L2 = rowMedians(sig_score_scRNAseq[,  Cluster_Cells$L2])
L3 = rowMedians(sig_score_scRNAseq[,  Cluster_Cells$L3])

sig_score_means = data.frame(D1 = D1, D2 = D2, D3 = D3, L1 = L1, L2 = L2, L3 = L3)
rownames(sig_score_means) = make.unique(rownames(sig_score_scRNAseq))


bio_themes_df = biological_themes_df[match(myDSP, biological_themes_df$pathway),]


set.seed(123);  biologicalTheme_heatmap_scRNA_mean = Heatmap(t(scale(t(sig_score_means[myDSP_sig,]))), 
                                                  col = colFun, #width = unit(5,'cm'),
                                                  name = "GSVA", row_title_rot = 0,
                                                  row_split = bio_themes_df[match(myDSP_sig, 
                                                                                  bio_themes_df$pathway), "signature_clusters"],
                                                  row_labels = gsub(x = myDSP_sig, pattern = "Hallmark - |KEGG - ", 
                                                                    replacement = ""),
                                                  show_column_names = T, show_row_names = F,  border = T, row_gap = unit(0,"cm"), 
                                                  column_gap = unit(0,'cm'),
                                                  row_names_gp = gpar(fontsize = 10, fontface = "bold"))

